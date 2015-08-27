#!/usr/bin/python
# encoding: utf-8
'''
cluster -- find similar biosynthetic clusters in multiple genomes.

First, cluster.py processes all input genbank files to associate genome sequence
with ID, species and strain information. If, during the analysis, a duplicate ID
is found, warning is issued and a numeric dot-suffix (e.g. ".1") is appended to
the ID. If all components are identical (ID, species, strain, and genome
length/contig count), then an exception is raised.
For each input file, nucleotide sequence is extracted as a [multi-]fasta file,
and annotated using antismash (--cpus 1, but up to as many parallel instances as
there are cores). This step ensures that there are no differences in the methods
of gene finding. Antismash annotation, depending on the options, will or will not
extend identified clusters. Translations of antismash-annotated CDS are then
extracted to a protein multi-fasta file. At the end of the 1st stage,
program maintains a relation between genome attributes (ID, etc) and 2 files:
antismash-generated genbank, and .faa.

Second, if orthology links information is requested, and after all the input files
have been processed as described above, an all-vs-all pblast is performed between
all possible genome pairs, and also each genome is pblasted against itself.
InParanoid is applied to all pairs of genomes to find orthologs. Finally,
quickparanoid is applied to generate a single table of orthologous groups.
These orthology-related steps can be skipped with --skip-orthology option.

Finally, cluster processes genomes and calculates similarity
between biosynthetic clusters in those genomes.

cluster.py outputs CSV files with collected information for further analysis.
'''


from __future__ import print_function
import sys
import logging
import os
import glob
import csv
import cPickle as pickle

from pprint import pprint
from os import mkdir, symlink, getcwd, chdir, remove
from os.path import exists, join, dirname, realpath, basename
from shutil import move
from argparse import ArgumentParser
from multiprocessing import Process, Queue, cpu_count
from itertools import permutations, combinations, product, combinations_with_replacement
from collections import defaultdict

from lib import utils
from lib.ClusterPair import ClusterPair, Cluster
from lib.Genome import Genome
from lib.MultiParanoid import MultiParanoid


#__all__ = []
__version__ = 0.7
__date__ = '2013-07-10'
__updated__ = '2014-05-11'


def preprocess_input_files(inputs, args):
    '''
    inputs: dict[genome.id] = Genome, is populated by this function.
    args.paths: list of genbank files to process.
    Read IDs and SeqRecords from GenBank files, make sure IDs are unique.
    Write 'all_clusters.csv' file with 3 columns: genomeid_cluster, Genome_ID, cluster.
    '''
    # Process genome IDs.
    for infile in args.paths:
        g = Genome(infile, args.project)
        # Check for duplicate ID.
        if g.id in inputs:
            # Check for an exact duplicate.
            if (inputs[g.id].contigs == g.contigs and
                inputs[g.id].genome_size == g.genome_size and
                inputs[g.id].species == g.species):
                logging.error('Found identical genomes of "%s" with ID %s, %s bp long (in %s fragments).',
                              g.species, g.id, g.genome_size, g.contigs)
                raise Exception('IdenticalGenomesError')
            else:
                for i in range(1, 999):
                    new_id = g.id + '_' + str(i) # TODO: zero-padded would look better
                    logging.debug('Appending %s to %s to make it unique: %s.', i, g.id, new_id)
                    if new_id in inputs:
                        continue
                    else:
                        g.id = new_id
                        del new_id
                        break
        inputs[g.id] = g

    total_genomes = len(inputs)
    workers = min(cpu_count(), total_genomes)
    task_queue = Queue()
    done_queue = Queue()
    # 1. Define worker.
    def geneparser(tasks, done, args):
        while True:
            try:
                # g = Genome object
                g = tasks.get() # By default, there is no timeout.
            except: # Queue.Empty
                logging.warning("geneparser(tasks, done) encountered an exception trying tasks.get()")
                raise
            # Exit if 'STOP' element is found.
            if g == 'STOP':
                logging.debug("STOP found, exiting.")
                break
            g.gb2fna()
            g.run_antismash(args.force, args.antismash_warning_shown,
                            args.no_extensions, cores = 1)
            # FIXME: cannot directly assign to args.antismash_warning_shown!!! shared object!!!
            #args.antismash_warning_shown = g.antismash_warning_shown
            g.parse_gene_cluster_relations(args)
            # Do not process zero-cluster genomes.
            if g.num_clusters() == 0:
                logging.info('%s (%s) has zero clusters, removing from further analysis',
                             g.species, g.id)
                done.put(None)
                continue
            g.as2faa(args.force)
            done.put(g)
    # 2. Start workers.
    workers_list = []
    logging.info("Starting %s Genome workers.", workers)
    for _ in range(workers):
        p = Process(target=geneparser, args=(task_queue, done_queue, args))
        p.start()
        workers_list.append(p)
    del p
    # 3. Populate tasks.
    logging.debug("Populating task_queue.")
    for g in inputs.itervalues():
        task_queue.put(g)
    # 4. Add STOP messages.
    logging.debug("Adding %s STOP messages to task_queue.", workers)
    for _ in range(workers):
        task_queue.put('STOP')
    del workers
    # 5. Start collecting results.
    # First, empty 'inputs'; this preserves reference to the external 'inputs'.
    for k in inputs.keys():
        del inputs[k]
    # Setup handle to a CSV file.
    # TODO: show a warning if overwriting an existing file.
    csvall = open(join(args.project, 'all_clusters.csv'), 'w')
    writer = csv.writer(csvall, delimiter = '\t', quoting = csv.QUOTE_NONE)
    # Output header.
    header = ['id', 'genome_ID', 'cluster', 'type']
    writer.writerow(header)
    for _ in range(total_genomes):
        g = done_queue.get()
        if g != None:
            inputs[g.id] = g
            # Output all clusters of this genome to the CSV file.
            for cluster in g.clusters:
                writer.writerow([''.join([g.id, str(cluster)]), g.id, cluster, g.number2products[cluster]])
        # Populate all_clusters.
        #for c in g.clusters:
        #    all_clusters.append(cluster(g.id, c))
    # Close CSV.
    csvall.close()
    # 6. Wait for processes to finish, close queues.
    logging.debug("Joining processes.")
    for p in workers_list:
        p.join()
    del workers_list, total_genomes
    task_queue.close()
    done_queue.close()


def prepare_inparanoid(inputs, args):
    '''
    Set up directory and .faa-files symlink for running InParanoid and
    quickparanoid. Prepare list of .faa-files w/o paths, and return it.
    '''
    # Collect all faafile names into a list without paths.
    faafiles = []
    for _ in inputs.itervalues():
        faafiles.append(basename(_.faafile))
    faafiles.sort(key = lambda s: s.lower())

    # Put everything inparanoid-related into a subdir.
    inparanoidir = realpath(join(args.project, 'inparanoid'))
    if not (args.force and exists(inparanoidir)):
        mkdir(inparanoidir)

    # Symlink all faa files there.
    for _ in faafiles:
        linkname = join(inparanoidir, _)
        if not (args.force and exists(linkname)):
            symlink(join('..', _), linkname)
    del linkname, _

    return inparanoidir, faafiles


def run_inparanoid(inparanoidir, faafiles, emulate_inparanoid):
    # TODO: make this into Inparanoid class, join with prepare_inparanoid.
    # Prepend custom_inparanoid path to PATH, so that it is used first.
    # alternative path finding method: dirname(sys.argv[0])
    custom_inparanoid = join(dirname(realpath(__file__)), 'custom_inparanoid')
    os.environ['PATH'] = custom_inparanoid + ':' + os.environ['PATH']
    logging.debug('PATH after prepending custom_inparanoid: %s', os.environ['PATH'])
    del custom_inparanoid

    # CD to the inparanoid directory.
    curr_path = getcwd()
    chdir(inparanoidir)
    logging.debug("Changed directory from %s to %s.", curr_path, inparanoidir)

    if emulate_inparanoid:
        num_workers = 1 # to avoid output mangling
    else:
        num_workers = cpu_count()
    qsize = 100 * num_workers

    total_genomes = len(faafiles)
    if not emulate_inparanoid:
        print("BLASTing %s single genomes." % total_genomes)
    tasks = Queue(maxsize = qsize)
    # 1. Define worker.
    def single_blaster(tasks, total_genomes):
        '''
        parallel worker for single blasts
        '''
        while True:
            try:
                # serial = counter, faaname = path to .faa file
                (serial, faaname) = tasks.get() # By default, there is no timeout.
            except: # Queue.Empty
                logging.warning("single_blaster(tasks) encountered an exception trying tasks.get()")
                raise
            # Exit if 'STOP' element is found.
            if faaname == 'STOP':
                logging.debug("STOP found, exiting.")
                break
            blast_single = ['inparanoid.pl', '--blast-only', faaname]
            if emulate_inparanoid:
                print(' '.join(blast_single))
            else:
                logging.info('Start blast %s / %s on a single genome: %s',
                             serial, total_genomes, ' '.join(blast_single))
                out, err, retcode = utils.execute(blast_single)
                logging.info(' Done blast %s / %s on a single genome: %s',
                             serial, total_genomes, ' '.join(blast_single))
                if retcode != 0:
                    logging.debug('inparanoid returned %d: %r while blasting %r, full output follows:\n%s',
                                  retcode, err, faaname, out)
                del blast_single, out, err, retcode
    # 2. Start workers.
    workers = []
    for _ in range(num_workers):
        p = Process(target=single_blaster, args=(tasks, total_genomes))
        workers.append(p)
        p.start()
    del p
    # 3. Populate the tasks queue.
    logging.debug("Populating the queue.")
    counter = 0
    for _ in faafiles:
        counter += 1
        tasks.put((counter, _)) # will block until all-qsize items are consumed
    del counter
    # 4. Add STOP messages.
    logging.debug("Adding %s STOP messages to task_queue.", num_workers)
    for _ in range(num_workers):
        tasks.put((0, 'STOP'))
    del _
    # 5. Wait for all processes to finish, close the queue.
    for p in workers:
        p.join()
    tasks.close()
    del total_genomes, workers

    total_permutations = len(list(permutations(faafiles, 2)))
    tasks = Queue(maxsize = qsize)
    if not emulate_inparanoid:
        print("BLASTing %s pairwise permutations." % total_permutations)
    # 1. Define worker.
    def pair_blaster(tasks, total_permutations):
        'parallel worker for paired blasts'
        while True:
            try:
                # serial = counter, faa1/2 = paths to .faa files
                (serial, faa1, faa2) = tasks.get() # By default, there is no timeout.
            except: # Queue.Empty
                logging.warning("pair_blaster(tasks) encountered an exception trying tasks.get()")
                raise
            # Exit if 'STOP' element is found.
            if faa1 == 'STOP':
                logging.debug("STOP found, exiting.")
                break
            blast_pair = ['inparanoid.pl', '--blast-only', faa1, faa2]
            if emulate_inparanoid:
                print(' '.join(blast_pair))
            else:
                logging.info('Start blast %s / %s on genome pair: %s', serial,
                             total_permutations, ' '.join(blast_pair))
                out, err, retcode = utils.execute(blast_pair)
                logging.info(' Done blast %s / %s on genome pair: %s', serial,
                             total_permutations, ' '.join(blast_pair))
                if retcode != 0:
                    logging.debug('inparanoid returned %d: %r while blasting %r and %r, full output follows:\n%s',
                                  retcode, err, faa1, faa2, out)
                del blast_pair, out, err, retcode
    # 2. Start workers.
    workers = []
    for _ in range(num_workers):
        p = Process(target=pair_blaster, args=(tasks, total_permutations))
        workers.append(p)
        p.start()
    del p
    # 3. Populate the tasks queue.
    logging.debug("Populating the queue.")
    counter = 0
    for pair in permutations(faafiles, 2):
        # First, formatdb input files.
        for f in pair:
            # Check if the formatdb file exists.
            if not exists(f + '.psq'):
                formatdb = ['formatdb', '-i', f]
                out, err, retcode = utils.execute(formatdb)
                if retcode != 0:
                    logging.warning('formatdb returned %d: %r while formatting %r, full output follows:\n%s',
                                    retcode, err, f, out)
                else:
                    logging.info('formatted %s: %s', f, ' '.join(formatdb) )
                del formatdb, out, err, retcode
        counter += 1
        tasks.put((counter, pair[0], pair[1])) # will block until all-qsize items are consumed
    del counter
    # 4. Add STOP messages.
    logging.debug("Adding %s STOP messages to task_queue.", num_workers)
    for _ in range(num_workers):
        tasks.put((0, 'STOP', ''))
    del _
    # 5. Wait for all processes to finish, close the queue.
    for p in workers:
        p.join()
    tasks.close()
    del total_permutations, workers

    if not emulate_inparanoid:
        total_combinations = len(list(combinations(faafiles, 2)))
        print("Analyzing with inparanoid %s pairwise combinations." % total_combinations)
        tasks = Queue(maxsize = qsize)
        # 1. Define worker.
        def inparanoider(tasks, total_combinations):
            '''
            parallel worker for paired blasts
            '''
            while True:
                try:
                    # serial = counter, faa1/2 = paths to .faa files
                    (serial, faa1, faa2) = tasks.get() # By default, there is no timeout.
                except: # Queue.Empty
                    logging.warning("inparanoider(tasks) encountered an exception trying tasks.get()")
                    raise
                # Exit if 'STOP' element is found.
                if faa1 == 'STOP':
                    logging.debug("STOP found, exiting.")
                    break
                inparanoid_pair = ['inparanoid.pl', faa1, faa2]
                logging.info('Start inparanoid %s / %s: %s', serial,
                             total_combinations, ' '.join(inparanoid_pair))
                out, err, retcode = utils.execute(inparanoid_pair)
                logging.info(' Done inparanoid %s / %s: %s', serial,
                             total_combinations, ' '.join(inparanoid_pair))
                if retcode != 0:
                    logging.debug('inparanoid returned %d: %r while analyzing %r and %r, full output follows:\n%s',
                                  retcode, err, faa1, faa2, out)
                del inparanoid_pair, out, err, retcode
        # 2. Start workers.
        workers = []
        for _ in range(num_workers):
            p = Process(target=inparanoider, args=(tasks, total_combinations))
            workers.append(p)
            p.start()
        del p
        # 3. Populate the tasks queue.
        logging.debug("Populating the queue.")
        counter = 0
        for pair in combinations(faafiles, 2):
            counter += 1
            tasks.put((counter, pair[0], pair[1])) # will block until all-qsize items are consumed
        del pair, counter
        # 4. Add STOP messages.
        logging.debug("Adding %s STOP messages to task_queue.", num_workers)
        for _ in range(num_workers):
            tasks.put((0, 'STOP', ''))
        # 5. Wait for all processes to finish, close the queue.
        for p in workers:
            p.join()
        tasks.close()
        del total_combinations, workers, num_workers
    else:
        logging.info('Skipped inparanoid analysis because of --emulate-inparanoid.')

    # Cleanup .phr, .pin, .psq formatdb output files after all inparanoid runs,
    # including analysis; parallel inparanoid does not do this.
    if not emulate_inparanoid:
        formatdb_extensions = ['*.phr', '*.pin', '*.psq']
        for ext in formatdb_extensions:
            for _ in glob.glob(join(inparanoidir, ext)):
                logging.debug("Deleting %s.", _)
                remove(_)
        del formatdb_extensions

    # CD back to the initial directory.
    chdir(curr_path)
    logging.debug("Changed directory from %s to %s.", inparanoidir, curr_path)
    del inparanoidir, curr_path


def run_quickparanoid(inparanoidir, faafiles, project):
    '''
    Requires a config file, which is simply a list of all .faa files, 1 per line.
    quickparanoid MUST BE RUN from quickparanoid dir!
    quickparanoid will generate a Makefile.in in its own directory.
    When all is configured, run 'qp inparanoidir configfile_path execfile_prefix'
    to generate EXEC_FILE and EXEC_FILES in the quickparanoid directory.
    Then
    ./EXEC_FILE > quickparanoid_result.txt
    and
    ./EXEC_FILEs for some stats
    Return absolute result_name path.
    '''
    configfile = realpath(join(project, 'quickparanoid.config'))
    logging.debug("Generating %s.", configfile)
    # TODO: possibly add a check for existing configfile and quickparanoid result,
    #       and skip this entirely? (can use args.force and exists() for checking).
    logging.debug("(configfile is re-generated even if it exists)")
    with open(configfile, 'w') as conf_handle:
        conf_handle.write("\n".join(faafiles))

    # chdir() to where quickparanoid will be run from.
    # Alternative path finding method: dirname(sys.argv[0])
    quickparanoid = realpath(join(dirname(realpath(__file__)), 'quickparanoid'))
    curr_path = getcwd()
    logging.debug("Remembering current directory %s, changing to %s.",
                  curr_path, quickparanoid)
    chdir(quickparanoid)

    # quickparanoid requires a trailing slash.
    qp = ['./qp', inparanoidir + os.sep, configfile, project]
    logging.info('Running quickparanoid analysis: %s', ' '.join(qp))
    out, err, retcode = utils.execute(qp)
    if retcode != 0 or not exists(join(quickparanoid, project)):
        logging.error('quickparanoid returned %d: %r while analyzing %r in %r',
                      retcode, err, configfile, project)
        logging.error('Full output:\n%s\n', out)
        raise Exception('QuickParanoidError')
    del out, err, retcode

    # Move generated 'project' and 'projects' executables to the {project} directory.
    move(join(quickparanoid, project), join(curr_path, project, project))
    move(join(quickparanoid, project + 's'),
         join(curr_path, project, project + 's'))

    # Delete leftover garbage from quickparanoid.
    for _ in ['dump', 'gen_header', 'hashtable_itr.o', 'ortholog.o', 'qp.h',
              '__ortholog.h']:
        remove(join(quickparanoid, _))

    # Get back from quickparanoid.
    project_dir = join(curr_path, project)
    logging.debug("Going to %s from %s.", project_dir, quickparanoid)
    chdir(project_dir)

    # Run generated executable to get results file.
    result_name = 'quickparanoid-' + project + '.txt'
    exe = ['./' + project]
    logging.info('Running generated executable: %s', ' '.join(exe))
    out, err, retcode = utils.execute(exe)
    if retcode != 0:
        logging.debug('executable returned %d: %r', retcode, err)
    else:
        with open(result_name, 'w') as fh:
            fh.write(out)
    del out, err, retcode

    logging.debug("Returning to %s from %s.", curr_path, project_dir)
    chdir(curr_path)
    del quickparanoid, curr_path
    return join(project_dir, result_name)


def main():
    '''Command line options.'''
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)

    parser = ArgumentParser()
    parser.add_argument("-d", "--debug", dest="debug", action="store_true", default=False, help="set verbosity level to debug [default: %(default)s]")
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False, help="report only warnings and errors [default: %(default)s]")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument("--trim", dest="trim", action="store_true", default=False, help="trim away antismash2 cluster extensions [default: %(default)s]")
    parser.add_argument("--fulldp", dest="fulldp", action="store_true", default=False, help="use full dynamic programming solution in usearch alignment (much slower!) [default: %(default)s]")
    parser.add_argument("--highmem", dest="highmem", action="store_true", default=False, help="assume huge RAM: all GenBanks are loaded early and kept in RAM (faster processing) [default: %(default)s]")
    parser.add_argument("--cutoff", dest="cutoff", type = float, default = 0.6, help="protein identity cut-off when aligning with usearch [default: %(default)s]")
    parser.add_argument("--skip-putative", dest="skipp", action="store_true", default=False, help="exclude putative clusters from the analysis [default: %(default)s]")
    parser.add_argument("--skip-orthology", action="store_true", default=False, help="do not run any orthology analysis [default: %(default)s]")
    parser.add_argument("--strict", dest="strict", action="store_true", default=False, help="orthology link weight between clusters with 5 and 10 genes will never exceed 0.5 [default: %(default)s]")
    parser.add_argument("--no-name-problems", dest="no_name_problems", action="store_true", default=False, help="only use ortho-clusters which do not have diff.names tree_conflict problems [default: %(default)s]")
    parser.add_argument("--no-tree-problems", dest="no_tree_problems", action="store_true", default=False, help="only use ortho-clusters which do not have [diff.names/diff.numbers] tree_conflict problems [default: %(default)s]")
    parser.add_argument('--emulate-inparanoid', action = 'store_true', default = False, help='only print generated inparanoid commands, do not run; exit after inparanoid blasting; suppress some normal output [default: %(default)s]')
    parser.add_argument("--prefix", default='out', help="output CSV files prefix [default: %(default)s]")
    parser.add_argument("--project", default='cluster_project', help="put all the project files into this directory [default: %(default)s]")
    parser.add_argument('--force', action = 'store_true', default = False, help='insist on re-using existing project directory (this will re-use existing intermediate files) [default: %(default)s]')
    parser.add_argument('--no-extensions', action = 'store_true', default = False, help='pass --no-extensions option to the modified antismash2 (see README for details) [default: %(default)s]')
    parser.add_argument('--threshold', action = 'store', type=float, default = 0.0, help='cluster links with weight below this one will be discarded [default: %(default)s]')
    parser.add_argument('--from-file', action = 'store', help='read paths to GenBank files (one per line) from the provided file')
    parser.add_argument(dest="paths", help="paths to GenBank files with genomes to analyze", metavar="path", nargs='*')
    args = parser.parse_args()

    # Special flag to skip detailed warnings about existing antismash files
    # after the first is shown.
    args.antismash_warning_shown = False

    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    else:
        level = logging.INFO
    logging.basicConfig(level=level,
                        format='%(asctime)s::%(levelname)s::%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    # Where do we take the paths from?
    if len(args.paths) > 0:
        logging.debug("GenBank paths were provided on the command line, using these.")
    if args.from_file != None:
        logging.debug("GenBank paths were provided in '%s', using these.",
                      args.from_file)
        if exists(args.from_file):
            with open(args.from_file) as from_handle:
                for _ in from_handle:
                    args.paths.append(_.strip())
        else:
            logging.exception("File %s does not exist.", args.from_file)
            raise Exception("FileDoesNotExistError")

    # Remove any possible duplicate input filenames.
    args.paths = list(set(args.paths))

    if len(args.paths) == 0:
        logging.error("No GenBank files were given to the program.")
        sys.exit(4)
    logging.info("Will process %s input files.", len(args.paths))

    # Useful for generating inparanoid commands to run on different computers.
    if args.emulate_inparanoid:
        print('Used arguments and options:', file = sys.stderr)
        pprint(vars(args), stream = sys.stderr)
    else:
        print('Used arguments and options:')
        pprint(vars(args))

    if exists(args.project):
        if args.force:
            logging.warning('Re-using existing project directory "%s" as requested.', args.project)
        else:
            logging.error('Specified project directory "%s" already exists! Use --force to continue anyway.', args.project)
            sys.exit(1)
    else: # create
        mkdir(args.project)

    # "Catalog" all genomes.
    genomes = {} # Map genome ID to Genome object.

    if len(args.paths) == 1:
        logging.warning("Single input file specified, program will exit after preprocessing.")
    preprocess_input_files(genomes, args)
    if len(args.paths) == 1: # single input - exit
        sys.exit(3)
    if not args.skip_orthology:
        inparanoidir, faafiles = prepare_inparanoid(genomes, args)
        run_inparanoid(inparanoidir, faafiles, args.emulate_inparanoid)
    if args.emulate_inparanoid:
        logging.info('Exiting prematurely because of --emulate-inparanoid.')
        return 0

    if not args.skip_orthology:
        result_path = run_quickparanoid(inparanoidir, faafiles, args.project)
        del inparanoidir, faafiles
        # Parse quickparanoid results.
        mp = MultiParanoid(result_path, args.no_tree_problems, args.no_name_problems)

    # Setup cache directories.
    ortholinks = realpath(join(args.project, 'ortholinks'))
    if not exists(ortholinks):
        mkdir(ortholinks)
    usearch = realpath(join(args.project, 'usearch_' + str(int(args.cutoff * 100))))
    if args.fulldp:
        usearch += '_fulldp'
    if not exists(usearch):
        mkdir(usearch)

    # Total (theoretical) number of cluster pairs considered.
    total_pairs = 0
    # Counter of cluster pairs submitted to usearch.
    submitted_tasks = 0
    # Counter of cluster pairs with at least 1 usearch result.
    cluster_pairs_counter = 0

    # Open the output CSV file for writing, prepare CSV writer.
    res_fname = 'results'
    if args.fulldp:
        res_fname += '_fulldp'
    if args.skip_orthology:
        res_fname += '_skip_orthology'
    if args.trim:
        res_fname += '_trim'
    if args.no_extensions:
        res_fname += '_no_extensions'
    res_fname += '.csv'
    csvout = open(join(args.project, res_fname), 'w')
    writer = csv.writer(csvout, delimiter = '\t', quoting = csv.QUOTE_NONE)
    # Output header.
    header = ['is_intra', 'genome1_ID', 'genome2_ID', 'species1', 'species2',
              'cluster1', 'cluster2', 'type1', 'type2', 'genes1_count',
              'genes2_count', 'size1_kb', 'size2_kb', 'ortholinks_count1',
              'ortholinks_count2', 'similar_genes_count',
              'avg_protein_identity', 'P', 'K', 'S']
    writer.writerow(header)
    del header

    # Every loaded genome (g1) is expected to be used many times before unloading.
    combinations_counter = 0 # to track progress, outer genomes loop only
    total_combinations = len(list(combinations_with_replacement(genomes.keys(), r = 2)))
    prev_g1 = None # tracker of the previous g1 genome, to know when to unload it

    # Preparing for parallel usearch.
    # 0. Define queues.
    tasks = Queue(maxsize = cpu_count() * 4)
    done  = Queue()
    # 1. Define worker.
    def usearcher():
        '''
        Run_usearch on a single NamedTemporaryFile seqfile from 'tasks' queue,
        put (CPid, (avg_identity, gene1_to_gene2, protein_identities)) into 'done'.
        '''
        while True:
            try:
                task = tasks.get()
            except:
                logging.exception('tasks queue empty in usearcher?')
                raise
            if task == 'STOP':
                #logging.debug('STOP found, exiting.')
                break
            # Calculate gene-level and clusterpair-average protein identities in clusters.
            CPid, seqfilename = task
            result = utils.run_usearch(seqfilename, args.cutoff, args.fulldp)
            remove(seqfilename)
            logging.debug('Deleted %s.', seqfilename)
            if result != None:
                done.put((CPid, result))
                continue
            done.put((CPid, None))
    # 2. Start workers.
    workers_list = []
    logging.debug('Starting %s cp_processors.', cpu_count())
    for _ in range(cpu_count()):
        p = Process(target=usearcher, args=())
        p.start()
        workers_list.append(p)
    del p

    for g1, g2 in combinations_with_replacement(genomes.keys(), r = 2):
        combinations_counter += 1
        print('%s / %s\t' % (combinations_counter, total_combinations), genomes[g1].id, genomes[g2].id)
        # Unload previous "second genome", if it is not equal to the "first genome".
        if not args.highmem and prev_g1 and g1 != prev_g1:
            genomes[prev_g1].unload()
            prev_g1 = g1
        # Dict of ClusterPairs, for this pair of genomes.
        cluster_pairs = {}
        genomes[g1].load()
        genomes[g2].load() # if g1 == g2, this will do nothing (g1 already loaded)
        genome1 = genomes[g1]
        genome2 = genomes[g2]

        # Load caches, if any.
        pair_id = '_'.join([g1, str(genome1.total_genes_in_clusters),
                            g2, str(genome2.total_genes_in_clusters)])
        ortho_file = join(ortholinks, pair_id)
        # Flag, if True - write cache back to the file.
        orthocache_updated = False
        if exists(ortho_file):
            with open(ortho_file) as ortho_handle:
                orthocache = pickle.load(ortho_handle)
        else:
            orthocache = {}
            orthocache_updated = True
        usearch_file = join(usearch, pair_id)
        # Similar flag for the usearch cache.
        usearchcache_updated = False
        if exists(usearch_file):
            with open(usearch_file) as usearch_handle:
                usearchcache = pickle.load(usearch_handle)
        else:
            usearchcache = {}
            usearchcache_updated = True
        del pair_id

        # Iterate all possible cluster pairs between these 2 genomes, generate cluster pairs.
        # When g1 == g2, multiple internal cluster comparisons (c1 to c2, then c2 to c1, etc) would happen with 'product'.
        if g1 == g2:
            # Generate only unique intra-species cluster-cluster pairs (combinations).
            pairslist = [(Cluster(g1, c1), Cluster(g1, c2))
                         for c1, c2 in combinations(genome1.clusters, r = 2)]
        else:
            # Pair each cluster from g1 with each cluster from g2 (product).
            pairslist = [(Cluster(g1, c1), Cluster(g2, c2))
                         for c1, c2 in product(genome1.clusters,
                                               genome2.clusters)]
        # 3. Populate tasks queue.
        logging.debug('Populating usearcher tasks queue.')
        for cl1, cl2 in pairslist:
            cp = ClusterPair(cl1, cl2)
            CPid = (cp.g1, cp.c1, cp.g2, cp.c2)
            if not args.skip_orthology:
                # Check cache.
                if CPid in orthocache:
                    cp.link1, cp.link2 = orthocache[CPid]
                else:
                    cp.assign_orthologous_link(mp, genomes, args)
                    orthocache[CPid] = (cp.link1, cp.link2)
                    usearchcache_updated = True
            if args.skip_orthology or cp.link1 > 0 or cp.link2 > 0:
                # Check cache.
                if CPid in usearchcache:
                    # usearch[CPid] = (cp.avg_identity, cp.gene1_to_gene2, cp.protein_identities)
                    done.put((CPid, usearchcache[CPid]))
                else:
                    # Prepare for running usearch.
                    seqfilename = cp.pre_CDS_identities(genome1, genome2)
                    tasks.put((CPid, seqfilename))
                    logging.debug('Submitted %s with %s.', CPid, seqfilename)
                    orthocache_updated = True
                # Save the cp object for re-use.
                cluster_pairs[CPid] = cp
        submitted_tasks += len(cluster_pairs)
        total_pairs += len(pairslist)
        # 4. Collect results and write them to file.
        logging.debug('Collecting results.')
        for _ in range(len(cluster_pairs)):
            result = done.get()
            logging.debug('Collected task %s of %s local.', _ + 1, len(cluster_pairs))
            if result[1] == None: # no result for this cluster pair
                # Still save this to the cache
                if result[0] not in usearchcache:
                    usearchcache[result[0]] = None
                continue
            # Get cp back from cluster_pairs dict.
            cp = cluster_pairs[result[0]]
            cp.avg_identity, cp.gene1_to_gene2, cp.protein_identities = result[1]
            if result[0] not in usearchcache: # CPid
                usearchcache[result[0]] = result[1]
            cluster_pairs_counter += 1
            #logging.debug('Average identity of %s and %s is %s ',
            #              cp.gc1, cp.gc2, cp.avg_identity)
            #print(cp.protein_identities)
            if cp.avg_identity[0] > 2:
                # Calculate gene order preservation (for similar genes).
                cp.gene_order(genome1, genome2)
            # Calculate predicted domains order preservation within similar genes.
            #cp.domains(genomes)
            # Calculate cluster-level nucleotide identity.
            #cp.nucleotide_similarity(genomes)
            writer.writerow([int(cp.intra), genome1.id, genome2.id,
                   genome1.species, genome2.species, cp.c1, cp.c2,
                   genome1.number2products[cp.c1],
                   genome2.number2products[cp.c2],
                   cp.num_c1_genes(genome1), cp.num_c2_genes(genome2),
                   genome1.clustersizes[cp.c1],
                   genome2.clustersizes[cp.c2], cp.link1, cp.link2,
                   cp.avg_identity[0], round(cp.avg_identity[1], 1),
                   round(cp.pearson, 2), round(cp.kendall, 2),
                   round(cp.spearman, 2)])
        logging.debug('Done collecting results.')

        if g1 != g2 and not args.highmem:
            genomes[g2].unload()
        # Save caches.
        if orthocache_updated:
            with open(ortho_file, 'w') as ortho_handle:
                pickle.dump(orthocache, ortho_handle, pickle.HIGHEST_PROTOCOL)
                #logging.debug('Dumped orthocache: %s', orthocache)
        if usearchcache_updated:
            with open(usearch_file, 'w') as usearch_handle:
                pickle.dump(usearchcache, usearch_handle, pickle.HIGHEST_PROTOCOL)
        del genome1, genome2, cluster_pairs

    # 5. Add STOP messages.
    logging.debug('Adding %s STOP messages to tasks.', cpu_count())
    for _ in range(cpu_count()):
        tasks.put('STOP')
    # 6. Wait for processes to finish, close queues.
    logging.debug('Joining processes.')
    for p in workers_list:
        p.join(timeout = 10)
    del workers_list
    tasks.close()
    done.close()

    print('Usearch non-empty for %s of %s cluster pairs (%s total).' %
          (cluster_pairs_counter, submitted_tasks, total_pairs))
    del submitted_tasks, cluster_pairs_counter
    csvout.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
