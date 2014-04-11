'''
Represents a single genome.
Has methods to annotate it, and generate derivative files.
Contains gene-cluster-product-coordinate mappings.
'''


from __future__ import print_function
import logging
import glob
from Bio import SeqIO
from Bio.Alphabet import generic_dna#, generic_protein
from shutil import rmtree, move
from os.path import join, exists
from multiprocessing import cpu_count
from lib import utils
from lib.gb2fasta import gb2fasta
from lib.extract_translation_from_genbank import extract_translation_from_genbank


def __init__(self, infile, project):
        '''
        infile: path to the original GenBank file for this genome
        project: path to the project directory, which has all the derivative files
        '''
        self.infile = infile
        self.project = project
        self.species = ''
        self.accession = '' # accession of the longest record
        self.id = '' # id of the longest record
        self.accessions = [] # other accessions
        self.ids = [] # other ids
        self.records = None # points at the indexed GenBank record(s)
        self.contigs = 0 # count of records/contigs/plasmids
        self.genome_size = 0 # sum of lengths of all records
        self.fnafile = '' # path to .fna file
        self.as2file = '' # path to antismash2 file
        self.is_annotated = False # antismash2 annotation
        self.num_clusters = 0 # number of antismash2 clusters
        self.faafile = ''
        self.number2products = {} # numeric clusters to their products
        self.coords2numbers = {} # coords of cluster to its number
        self.cluster2genes = {} # number to a list
        self.gene2cluster = {} # gene ID to cluster number
        self.clustersizes = {} # cluster number to cluster size, bps
        self.is_genecluster_parsed = False
        # These values are for the parent to check.
        self.antismash2_reused = False
        self.antismash_warning_shown = False

        organism = {} # maps accessions to 'organism' field
        primary_length = 0 # current primary accession sequence length
        for r in SeqIO.parse(infile, 'genbank', generic_dna):
            # r.name is an accession (e.g. AF080235),
            # r.id is a versioned accession (e.g. AF080235.1)
            # We use r.id as more specific.
            logging.debug('record name and ID are %s and %s', r.name, r.id)
            if not r.id or r.id == '':
                logging.error('Genome %s has no usable ID!' % input)
                raise Exception('GenomeHasNoIDError')
            self.contigs += 1
            self.accessions.append(r.name)
            self.ids.append(r.id)
            organism[r.name] = r.annotations['organism']
            if self.id == '':
                self.accession = r.name
                self.id = r.id
                primary_length = len(r.seq)
            elif len(r.seq) > primary_length:
                self.accession = r.name
                self.id = r.id
                primary_length = len(r.seq)
            self.genome_size += len(r.seq)
        if '-' in self.id:
            logging.error('quickparanoid allows no dashes in filenames; genome "%s" has a dash!',
                          self.id)
            raise Exception('DashesInIDsAreForbiddenError')
        del primary_length, r
        # Remove primary accession from the list of all accessions.
        self.accessions.remove(self.accession)
        self.ids.remove(self.id)




        faafile = join(args.project, primary_id + '.faa')
        if antismash2_reused and args.force and exists(faafile):
            logging.debug('Reusing existing translations file.')
        else:
            extract_translation_from_genbank(as2file, faafile, False)
        del as2file, fnafile, faafile, infile
        del contigs, genome_size, organism, primary_accession, accessions
        del ids, primary_id


def gb2fna(self):
    # Write FASTA to ID.fna, do nothing if file exists.
    self.fnafile = join(self.project, self.id + '.fna')
    gb2fasta(self.infile, self.fnafile, self.id)


def run_antismash(self, force, antismash_warning_shown, no_extensions):
    '''
    force: will re-use existing antismash2 annotation
    antismash_warning_shown: flag, whether antismash2 re-use warning had already been shown
    no_extensions: run antismash with --no-extensions
    '''
    # Re-use antismash2 annotation if it exists.
    self.as2file = join(self.project, self.id + '.gbk')
    self.antismash2_reused = False
    if force and exists(self.as2file):
        if not antismash_warning_shown:
            self.antismash_warning_shown = True
            logging.warning('Reusing existing antismash2 annotation(s).')
            logging.warning('\t--no-extensions option will NOT be honored!')
            logging.warning('\tDo not use --force, or delete antismash2 *.gbk files to re-run antismash2 annotation.')
        # Re-using any further files is only possible if we do re-use antismash files.
        self.antismash2_reused = True
    else:
        output_folder = join(self.project, self.id)
        as2_options = ['run_antismash', '--outputfolder', output_folder]
        # Removed time-consuming/unnecessary options: smcogs, clusterblast,
        # subclusterblast.
        # Removed unnecessary options: --full-blast, --full-hmmer.
        # TODO: when this part is properly parallelized, set --cpus to 1?
        as2_options.extend(['--cpus', str(cpu_count()), '--verbose', '--all-orfs'])
        as2_options.extend(['--input-type', 'nucl', '--inclusive'])
        if no_extensions:
            # TODO: make this warning be shown only once, similar to antismash_warning_shown
            logging.warning("using --no-extensions option; cluster.py DOES NOT check if this option is supported!")
            as2_options.append('--no-extensions')
        as2_options.append(self.fnafile)
        logging.info('Running antismash2: %s', ' '.join(as2_options))
        # TODO: show output from child processes?...
        out, err, retcode = utils.execute(as2_options)
        if retcode != 0:
            logging.debug('antismash2 returned %d: %r while scanning %r',
                          retcode, err, self.fnafile)
        # antismash's algorithm for naming the output file:
        # basename = seq_records[0].id
        # output_name = path.join(options.outputfoldername, "%s.final.gbk" % basename)
        # However, seq_records[0].id can be dot-truncated by BioPython to length 12;
        # thus using wildcard match.
        antismash2_file = glob.glob(join(output_folder, '*.final.gbk'))[0]
        move(antismash2_file, self.as2file)
        rmtree(output_folder)
        del output_folder, out, err, retcode


def parse_gene_cluster_relations(self):
    pass
