#!/usr/bin/python
# encoding: utf-8
'''
cluster -- shortdesc

cluster is a description

It defines classes_and_methods

@author:     user_name

@copyright:  2013 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''


import sys
import os
import csv


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from Bio import SeqIO
from bx.intervals.intersection import Interval, IntervalTree


__all__ = []
__version__ = 0.1
__date__ = '2013-07-10'
__updated__ = '2013-07-10'


DEBUG = 0
TESTRUN = 0
PROFILE = 0


def find_ortho_cluster(g):
    'returns the number of the cluster, to which gene "g" belongs'
    return


def parse_cluster_number(note):
    'given a list of items from "note" field, return cluster number'
    for i in note:
        if i.startswith('Cluster number: '):
            return int(i[16:])


def process(paths):
    'main method'
    # TODO: check if output file exists to avoid overwriting it.
    outfile = paths[0]
    speciesfile = paths[1]
    orthofile = paths[2]
    genbanks = paths[3:]


    # Declare important variables.
    # Mapping of each gene to the clusterID it belongs to.
    gene2ortho = {}
    # List of all the genes in the cluster with given clusterID.
    ortho2genes = {}
    # List of recognized species.
    species = []
    # List of all clusters from all species.
    all_clusters = []
    # Map of per-species cluster numbers to their products
    numbers2products = {}
    # Dict of per-species GenBank records.
    genbank = {}
    # Dict of per-species interval trees of clusters.
    clustertrees = {}
    # Dict of per-species dicts of cluster-to-gene relations.
    # Each cluster dict uses a key = (species, number).
    cluster2genes = {}
    # Per-species mapping of cluster coordinates tuple to their numbers.
    coords2numbers = {}


    print 'Reading multiparanoid gene clusters:'
    with open(orthofile) as tsv:
        # Sample line:
        # 1    avermitilis_MA4680.faa    SAV_1680.BA000030    1    1.000    Kitasatospora_setae_DSM43861.faa-SirexAA_E.faa-albus_J1074.faa-avermitilis_MA4680.faa-cattleya_DSM46488.faa-coelicolor_A3_2.faa-flavogriseus_IAF45CD.faa-griseus_NBRC13350.faa-scabiei_87.22.faa-venezuelae_Shinobu_719.faa-violaceusniger_Tu4113.faa    diff. numbers
        # fieldnames = ['#clusterID', 'species', 'gene', 'is_seed_ortholog', 'confidence', 'species_in_cluster', 'tree_conflict']
        reader = csv.DictReader(tsv, delimiter = '\t')
        try:
            for row in reader:
                # Build gene-to-cluster-ID dict
                gene2ortho[row['gene']] = row['#clusterID']
                # Build cluster-ID-to-all-genes dict of lists.
                if row['#clusterID'] in ortho2genes:
                    ortho2genes[row['#clusterID']].append(row['gene'])
                else:
                    ortho2genes[row['#clusterID']] = [row['gene']]
        except csv.Error as e:
            sys.exit('file %s, line %d: %s' % (orthofile, reader.line_num, e))
    print '\ttotal lines read:', reader.line_num
    print '\ttotal entries in gene2ortho:', len(gene2ortho)
    print '\ttotal entries in ortho2genes:', len(ortho2genes)


    print 'Reading species list file:'
    for s in open(speciesfile):
        species.append(s.strip())
    print '\ttotal species:', len(species)
    species.sort(key=lambda s: s.lower())


    print 'Reading all genbank files:'
    # Detect which species it is by the first 5 characters.
    # FIXME: needs a better solution.
    for gb in genbanks:
        for s in species:
            if s[0:5] == gb[0:5]:
                print '\t%s corresponds to species %s' % (gb, s)
                break
        genbank[s] = SeqIO.read(gb, "genbank")
    print '\ttotal records parsed:', len(genbank)
#    print genbank.keys()


    print 'Parsing clusters and assigning genes to them:'
    for s in species:
        clustertrees[s] = IntervalTree()
        cluster2genes[s] = {}
        numbers2products[s] = {}
        coords2numbers[s] = {}
        # Populate clusters tree and dict with (start, end) as keys.
        for f in genbank[s].features:
            if f.type == 'cluster':
                start = int(f.location.start.position)
                end = int(f.location.end.position)
                clustertrees[s].add_interval(Interval(start, end))
                cluster_number = parse_cluster_number(f.qualifiers['note'])
                cluster2genes[s][cluster_number] = []
                all_clusters.append((s, cluster_number))
                numbers2products[s][cluster_number] = f.qualifiers['product'][0]
                coords2numbers[s][(start, end)] = cluster_number
#                print '\tadding cluster %s (%s) at (%s, %s)' % (cluster_number, f.qualifiers['product'][0], start, end)
        # Assign genes to each cluster.
        num_genes = 0
        for f in genbank[s].features:
            if f.type == 'CDS':
                # Determine which qualifier type to use for gene name.
                if 'locus_tag' in f.qualifiers:
                    qualifier = 'locus_tag'
                elif 'gene' in f.qualifiers:
                    qualifier = 'gene'
                # 'cl' is a list of Intervals, each has 'start' and 'end' attributes.
                cl = clustertrees[s].find(int(f.location.start.position), int(f.location.end.position))
                if len(cl) > 0:
                    # One gene may belong to more than one cluster.
#                    print 'gene at (%s, %s) overlaps with %s cluster(s), 1st is at (%s, %s)' % (f.location.start.position, f.location.end.position, len(cl), cl[0].start, cl[0].end)
                    num_genes += 1
                    for cluster in cl:
                        cluster2genes[s][coords2numbers[s][(cluster.start, cluster.end)]].append(f.qualifiers[qualifier][0] + '.' + genbank[s].name)
        print '\t%s: %s clusters populated with %s genes' % (s, len(cluster2genes[s]), num_genes)
    print '\tadded %s clusters from %s species' % (len(all_clusters), len(species))
    # Sort to ensure stable order of cluster pairs.
    all_clusters.sort(key=lambda s: s.lower())

    print 'Freeing memory.'
    del genbank
    del clustertrees


    print 'Constructing gene-based cluster links. This may take a while.'
    # Iterate all clusters.
    for c in all_clusters:
        # Iterate each gene in the current cluster.
        print 'cluster', c
        for g in cluster2genes[c[0]][c[1]]:
            print 'gene', g
            # List of clusters 'g' is linked to.
            g_map = []
            # Iterate list of all genes from the g's orthology cluster, except for 'g' itself.
            for xeno_gene in ortho2genes[gene2ortho[g]].remove(g):
                g_map.append()

        # weights calculation:
        # 1. get all clusters we have links to from our current cluster (just a union of all gene->cluster links within current cluster)
        # 2. iterate each of those linked clusters
        # 3. for each one, check each gene of the current cluster to see if it has a link out; increment counter if yes
        # 4. calculate weight for this link.
        #
        # link weight formula:
        # 1. if num_links/min(num_genes_cl1, num_genes_cl2) == 1.0, then weight is 1.0
        # 2. otherwise, weight = (num_links/2) * ( 1/min(num_genes_cl1, num_genes_cl2) + 1/max(num_genes_cl1, num_genes_cl2) )
        #
        # map of 'c' to other clusters, with weights
        #
        # resolve single-species multi-mapping clusters by weight (sp1.a->sp2.b 0.5, sp1.a->sp2.c 0.9, then sp1.a-> sp2.c)
        #
        # prune zero-weight links and links below the threshold (default threshold is 0)
        #
        # sort DESC from clusters with 11 links down to unique clusters with 0 links
        #
        # using cluster products, for each group of same-links clusters output tables for them:
        # Clusters with 11 links (a total of XXX):
        # sp1 sp2 sp3 ...
        # cl1 cl2 cl3 ...
        # ...
        #
        # After each such table, output a diagonal matrix of link weights between all possible cluster pairs.
        # This is done for the entire set of clusters from the summary table (so with 1 row of 11-linked clusters,
        # we'll show link weights between 55 pairs).
        # Cluster link weights:
        #    cl1 cl2 cl3 cl4
        # cl1 -  0.5 0.6 0.9
        # cl2     -  0.7 0.5
        # cl3         -  0.8
        # cl4             -
        # Once per species, show a table of intra-species cluster links, with weights;
        # looks just like the above weights table, but now only clusters from single
        # species are used.



#class CLIError(Exception):
#    '''Generic exception to raise and log different fatal errors.'''
#    def __init__(self, msg):
#        super(CLIError).__init__(type(self))
#        self.msg = "E: %s" % msg
#    def __str__(self):
#        return self.msg
#    def __unicode__(self):
#        return self.msg


def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2013 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

#    try:
        # Setup argument parser
    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
#        parser.add_argument("-r", "--recursive", dest="recurse", action="store_true", help="recurse into subfolders [default: %(default)s]")
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %(default)s]")
#        parser.add_argument("-i", "--include", dest="include", help="only include paths matching this regex pattern. Note: exclude is given preference over include. [default: %(default)s]", metavar="RE" )
#        parser.add_argument("-e", "--exclude", dest="exclude", help="exclude paths matching this regex pattern. [default: %(default)s]", metavar="RE" )
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument(dest="paths", help="paths to input files: out.txt config multiparanoid.txt 1.gb 2.gb ..", metavar="path", nargs='+')

    # Process arguments
    args = parser.parse_args()

    paths = args.paths
#        verbose = args.verbose
#        recurse = args.recurse
#        inpat = args.include
#        expat = args.exclude

#        if verbose > 0:
#            print("Verbose mode on")
#            if recurse:
#                print("Recursive mode on")
#            else:
#                print("Recursive mode off")

#        if inpat and expat and inpat == expat:
#            raise CLIError("include and exclude pattern are equal! Nothing will be processed.")

#        for inpath in paths:
#            ### do something with inpath ###
#            print(inpath)
    process(paths)
    return 0
#    except KeyboardInterrupt:
#        ### handle keyboard interrupt ###
#        return 0
#    except Exception, e:
#        if DEBUG or TESTRUN:
#            raise(e)
#        indent = len(program_name) * " "
#        sys.stderr.write(program_name + ": " + repr(e) + "\n")
#        sys.stderr.write(indent + "  for help use --help")
#        return 2

if __name__ == "__main__":
#    if DEBUG:
#        sys.argv.append("-h")
#        sys.argv.append("-v")
#        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'cluster_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())
