#!/usr/bin/python
'''
Given a bunch of GenBank files annotated with antismash,
find all paired occurrences of a regulator and a transporter
pointing in different directions (i.e. <-regulator| |transporter->
or <-transporter| |regulator->)
Lots of code borrowed from cluster.py, could probably be reorganized
into a library.
'''


from __future__ import print_function
import sys
import logging
import csv
from os.path import exists, basename#, join, dirname, realpath
from argparse import ArgumentParser

from lib.Genome import Genome


__version__ = 0.1
__date__ = '2014-10-28'
__updated__ = '2014-10-28'


def is_regulator(f):
    '''
    given feature 'f',
    return True if it is a regulator/repressor,
    False otherwise
    '''
    for substring in ('regulator', 'repressor', 'transcription'):
        if 'product' in f.qualifiers and substring in f.qualifiers['product'][0].lower():
            return True
        if 'label' in f.qualifiers and substring in f.qualifiers['label'][0].lower():
            return True
        if 'note' in f.qualifiers:
            for note in f.qualifiers['note']:
                if substring in note:
                    return True
    return False


def is_transporter(f):
    '''
    given feature 'f',
    return True if it is a transporter,
    False otherwise
    '''
    for substring in ('transport', 'abc', 'pump', 'transmembrane', 'trans-membrane'):
        if 'product' in f.qualifiers and substring in f.qualifiers['product'][0].lower():
            return True
        if 'label' in f.qualifiers and substring in f.qualifiers['label'][0].lower():
            return True
        if 'note' in f.qualifiers:
            for note in f.qualifiers['note']:
                if substring in note:
                    return True
    return False


def has_transporter_nearby(g, regulator, cluster_number):
    '''
    g: Genome object;
    f: feature (regulator/repressor)
    search for a transporter up- or down-stream from 'f', AND
    on the other strand; return either None, or the found feature
    '''
    # find the index of 'regulator' in orderstrands
    index_found = False # safeguard variable
    for regulator_index, current in enumerate(g.orderstrands[cluster_number]):
        # 'current' is a GeneOrder named tuple,
        # with .start, .strand, and .geneid = locus_tag:genome_id
        if current.start == regulator.location.start.position and current.strand == regulator.strand:
            # yes, this is our CDS!
            index_found = True
            break
    assert index_found
    # Examine previous index.
    if regulator_index > 0: # so that after -1 it will still be a valid index
        geneid = g.orderstrands[cluster_number][regulator_index - 1].geneid
        locus_tag = geneid.split(':')[0]
        prev = g.get_feature_by_locustag(locus_tag).feature
        if prev.strand != regulator.strand and is_transporter(prev):
            return prev
    # Examine next index.
    # TODO: almost identical to the code above, need to refactor.
    if regulator_index + 1 < len(g.orderstrands[cluster_number]): # so that after +1 it will still be in range
        geneid = g.orderstrands[cluster_number][regulator_index + 1].geneid
        locus_tag = geneid.split(':')[0]
        next_feature = g.get_feature_by_locustag(locus_tag).feature
        if next_feature.strand != regulator.strand and is_transporter(next_feature):
            return next_feature


def main():
    '''Command line options.'''
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)

    parser = ArgumentParser()
    parser.add_argument("-d", "--debug", dest="debug", action="store_true", default=False, help="set verbosity level to debug [default: %(default)s]")
    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true", default=False, help="report only warnings and errors [default: %(default)s]")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument("--prefix", default='out', help="output CSV files prefix [default: %(default)s]")
    parser.add_argument('--from-file', action = 'store', help='read paths to GenBank files (one per line) from the provided file')
    parser.add_argument(dest="paths", help="paths to the GenBank files with genomes to analyze", metavar="path", nargs='*')
    args = parser.parse_args()

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

    # Map genome IDs to Genome objects.
    inputs = dict()
    # make loaded genbanks stay in RAM until explicitly unloaded
    args.highmem = True

    # counters
    total_genes_in_clusters = 0
    gene_pairs = 0

    # prepare CSV for output
    csvout = open(args.prefix + '_gene_pairs.csv', 'w')
    writer = csv.writer(csvout, delimiter = '\t', quoting = csv.QUOTE_NONE)
    # Output header.
    header = ['filename', 'genomeID', 'species',
              'cluster_number', 'cluster_type', 'cluster_start', 'cluster_end',
              'regulator', 'regulator_start', 'regulator_end', 'regulator_strand',
              'transporter', 'transporter_start', 'transporter_end', 'transporter_strand']
    writer.writerow(header)

    # Catalog all genomes.
    for infile in args.paths:
        g = Genome(infile, '.', assume_infile_after_as2 = True)
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
        # perform actual processing
        g.load()
        g.parse_gene_cluster_relations(args)
        # iterate all clusters
        for cluster in g.clusters:
            # iterate all the CDS of a cluster
            for CDS in g.cluster2genes[cluster]:
                total_genes_in_clusters += 1
                # select proper record/contig and proper feature by index <- by locus_tag
                (record_index, feature_index) = g.CDS[CDS.split(':')[0]]
                regulator = g.records[ record_index ].features[ feature_index ]
                if is_regulator(regulator):
                    print('"%s" is a regulator...' % regulator.qualifiers['product'][0])
                    # check if there is a transporter nearby
                    transporter = has_transporter_nearby(g, regulator, cluster)
                    if transporter is not None:
                        print('--> Found transporter! "%s" and "%s"!' % (regulator.qualifiers['product'][0],
                                                                         transporter.qualifiers['product'][0]))
                        gene_pairs += 1
                        (c_record_index, c_feature_index) = g.clusterindex[cluster]
                        cluster_feature = g.records[c_record_index].features[c_feature_index]
                        row = [basename(infile),
                               g.id,
                               g.species,
                               cluster,
                               g.number2products[cluster],
                               cluster_feature.location.start.position,
                               cluster_feature.location.end.position,
                               regulator.qualifiers['product'][0],
                               regulator.location.start.position,
                               regulator.location.end.position,
                               regulator.strand,
                               transporter.qualifiers['product'][0],
                               transporter.location.start.position,
                               transporter.location.end.position,
                               transporter.strand]
                        writer.writerow(row)
                    else:
                        print('\t... which has no transporter nearby.')
        # clear RAM
        g.unload()
        # this might be unnecessary?..
        inputs[g.id] = g

    csvout.close()
    print('Total genes in clusters:', total_genes_in_clusters)
    print('Gene pairs found:', gene_pairs)

if __name__ == "__main__":
    sys.exit(main())
