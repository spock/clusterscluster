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
from multiprocessing import Process, Queue, cpu_count
from bx.intervals.intersection import Interval, IntervalTree
from lib import utils
from lib.gb2fasta import gb2fasta
from lib.extract_translation_from_genbank import extract_translation_from_genbank
# Extension trimming is implemented as --no-extensions in antismash2,
# making this import obsolete.
from lib import hmm_detection


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
        self.total_genes_in_clusters = 0 # number of genes in all clusters
        self.fnafile = '' # path to .fna file
        self.as2file = '' # path to antismash2 file
        self.is_annotated = False # antismash2 annotation
        self.num_clusters = 0 # number of antismash2 clusters
        self.clusters = [] # list of (numeric) cluster IDs
        self.faafile = ''
        self.number2products = {} # numeric clusters to their products
        self.coords2numbers = {} # coords of cluster to its number
        self.cluster2genes = {} # number to a list
        self.gene2clusters = {} # gene ID to a list of cluster numbers
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


def as2faa(self, force):
    self.faafile = join(self.project, self.id + '.faa')
    if self.antismash2_reused and force and exists(self.faafile):
        logging.debug('Reusing existing translations file.')
    else:
        extract_translation_from_genbank(self.as2file, self.faafile, False)


def gb2fna(self):
    # Write FASTA to ID.fna, do nothing if file exists.
    self.fnafile = join(self.project, self.id + '.fna')
    gb2fasta(self.infile, self.fnafile, self.id)


def run_antismash(self, force, antismash_warning_shown, no_extensions,
                  cores = cpu_count()):
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
        as2_options.extend(['--cpus', str(cores), '--verbose', '--all-orfs'])
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


def parse_cluster_number(note):
    '''
    Given a list of items from the "note" field of the GenBank feature,
    created by Antismash2, return cluster number.
    '''
    for i in note:
        if i.startswith('Cluster number: '):
            return int(i[16:])


def parse_gene_cluster_relations(self, args):
    '''
    Assign genes to clusters, clusters to genes, and fill other useful maps and
    lists.
    args: original arguments namespace from parent/main (args.trim etc)
    '''
    rulesdict = hmm_detection.create_rules_dict()
    logging.info('Parsing antismash2 genbank file %s.', self.as2file)

    # Dict of per-SeqRecord Interval trees of clusters, i.e.
    # clustertrees[record_number] = IntervalTree()
    clustertrees = {}
    # List of all records. TODO: make this an index, and put into self?
    records = list(SeqIO.parse(self.infile, "genbank", generic_dna))
    print('Parsing clusters and assigning genes to them in %s (%s):' %
          (self.id, self.species))
    if args.trim:
        logging.info('\tgetting extension sizes for diff. cluster types from antismash2 config')
        logging.info('\tNote: cluster coordinates are shown after trimming, except for skipped putative clusters.')
    # Populate clusters tree and dict with (start, end) as keys.
    for r in records:
        record_number = records.index(r)
        clustertrees[record_number] = IntervalTree()
        for f in r.features:
            if f.type == 'cluster':
                cluster_number = parse_cluster_number(f.qualifiers['note'])
                start = int(f.location.start.position)
                end = int(f.location.end.position)
                if args.skipp and f.qualifiers['product'][0] == 'putative':
                    logging.debug('\tskipping putative cluster #%s at (%s, %s)',
                                  cluster_number, start, end)
                    continue
                # Putative clusters have neither extensions nor rules for them.
                if args.trim and f.qualifiers['product'][0] != 'putative':
                    # Use cluster type to get extension size, including composite types.
                    extension = hmm_detection.get_extension_size_by_cluster_type(f.qualifiers['product'][0],
                                                                                 rulesdict)
                    # If cluster is shorter than double extension, do not trim.
                    if 2 * extension >= (end - start):
                        logging.info('\tnot trimming #%s (%s) - shorter than double extension 2 * %s',
                                     cluster_number, f.qualifiers['product'][0],
                                     extension)
                    else:
                        # If cluster is at the genome edge - skip extension trimming.
                        if start == 0:
                            logging.info('\tnot trimming left-extension for #%s (%s) - at the genome start',
                                         cluster_number, f.qualifiers['product'][0])
                        else:
                            start += extension
                        if end == len(r):
                            logging.info('\tnot trimming right-extension for #%s (%s) - at the genome end',
                                         cluster_number, f.qualifiers['product'][0])
                        else:
                            end -= extension
                        try:
                            assert start < end
                        except:
                            logging.exception('trimming extension failed for: %s (%s), extension %s, ori start %s, new start %s, ori end %s, new end %s',
                                              f.qualifiers['product'][0],
                                              cluster_number, extension,
                                              f.location.start.position,
                                              start,
                                              f.location.end.position,
                                              end)
                            raise
                clustertrees[record_number].add_interval(Interval(start,
                                                                  end))
                self.cluster2genes[cluster_number] = []
                self.clusters.append(cluster_number)
                self.number2products[cluster_number] = f.qualifiers['product'][0]
                # Start/end coords theoretically can collide between
                # multiple loci of the same genome: check for this.
                try:
                    assert (start, end) not in self.coords2numbers
                except:
                    logging.exception("Cluster %s coordinates (%s, %s) collision in %s.",
                                      cluster_number, start, end, self.id)
                    raise
                self.coords2numbers[(start, end)] = cluster_number
                self.clustersizes[cluster_number] = end - start
                logging.info('''\t('%s', %s): %s [%s, %s], %s bp long''',
                             self.id, cluster_number, f.qualifiers['product'][0],
                             start, end, end-start)
    logging.debug('\tAssign genes to biosynthetic clusters')
    for r in records:
        record_number = records.index(r)
        for f in r.features:
            if f.type == 'CDS':
                # 'cl' is a list of Intervals, each has 'start' and 'end' attributes.
                cl = clustertrees[record_number].find(int(f.location.start.position),
                                                      int(f.location.end.position))
                if len(cl) > 0:
                    # Determine which qualifier type to use for gene name.
                    if 'locus_tag' in f.qualifiers:
                        qualifier = 'locus_tag'
                    elif 'gene' in f.qualifiers:
                        qualifier = 'gene'
                    gene_name = f.qualifiers[qualifier][0] + ':' + self.id
                    # Check for gene_name overlaps between loci.
                    try:
                        assert gene_name not in self.gene2clusters
                    except:
                        logging.exception("%s is not unique across the records of %s.",
                                          gene_name, self.id)
                        raise
                    # One gene may belong to more than one cluster.
                    logging.debug('\t\tgene at (%s, %s) overlaps with %s cluster(s), 1st is at (%s, %s)',
                                  f.location.start.position,
                                  f.location.end.position, len(cl),
                                  cl[0].start, cl[0].end)
                    self.total_genes_in_clusters += 1
                    self.gene2clusters[gene_name] = []
                    for cluster in cl:
                        cluster_serial = self.coords2numbers[(cluster.start, cluster.end)]
                        self.cluster2genes[cluster_serial].append(gene_name)
                        self.gene2clusters[gene_name].append(cluster_serial)
                    del gene_name, cluster_serial
    print('\t%s: %s clusters populated with %s genes' % (self.id,
                                                         len(self.cluster2genes),
                                                         self.total_genes_in_clusters))
    del clustertrees, records
