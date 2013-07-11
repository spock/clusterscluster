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


def parse_cluster_number(note):
    'given a list of items from "note" field, return cluster number'
    for i in note:
        if i.startswith('Cluster number: '):
            return int(i[16:])


def process(paths, threshold = 0.0):
    'main method'
    # TODO: check if output file exists to avoid overwriting it.
    outfile = paths[0]
    speciesfile = paths[1]
    orthofile = paths[2]
    genbanks = paths[3:]


    # Declare important variables.
    # Mapping of each gene to the clusterID(s) it belongs to.
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
    # Inverse of the above: per-species mapping of each gene to the cluster(s) it belongs to.
    gene2clusters = {}
    # Per-species mapping of cluster coordinates tuple to their numbers.
    coords2numbers = {}
    # 2-level nested dict of cluster pairs link weights, e.g. cluster_weights['A'] = {'B': 0.95, ...}
    cluster_weights = {}
    # Same as above, but without duplicate links to other species.
    weights_clean = {}
    # Same as cluster_weights, but only for intra-species links.
    weights_intra = {}
    # Mapping of record.names from GenBank files (LOCUS) to species.
    locus2species = {}


    print 'Reading multiparanoid gene clusters:'
    with open(orthofile) as tsv:
        # Sample line:
        # 1    avermitilis_MA4680.faa    SAV_1680.BA000030    1    1.000    Kitasatospora_setae_DSM43861.faa-SirexAA_E.faa-albus_J1074.faa-avermitilis_MA4680.faa-cattleya_DSM46488.faa-coelicolor_A3_2.faa-flavogriseus_IAF45CD.faa-griseus_NBRC13350.faa-scabiei_87.22.faa-venezuelae_Shinobu_719.faa-violaceusniger_Tu4113.faa    diff. numbers
        # fieldnames = ['#clusterID', 'species', 'gene', 'is_seed_ortholog', 'confidence', 'species_in_cluster', 'tree_conflict']
        reader = csv.DictReader(tsv, delimiter = '\t')
        try:
            for row in reader:
                # Understand both Multi- and Quick-paranoid files.
                if '#clusterID' in row:
                    # QuickParanoid
                    idname = '#clusterID'
                else:
                    # MultiParanoid
                    idname = 'clusterID'
                # Build gene-to-cluster-ID(s) dict
                if row['gene'] not in gene2ortho:
                    gene2ortho[row['gene']] = []
                gene2ortho[row['gene']].append(row[idname])
                # Build cluster-ID-to-all-genes dict of lists.
                if row[idname] in ortho2genes:
                    ortho2genes[row[idname]].append(row['gene'])
                else:
                    ortho2genes[row[idname]] = [row['gene']]
        except csv.Error as e:
            sys.exit('file %s, line %d: %s' % (orthofile, reader.line_num, e))
    print '\ttotal lines read:', reader.line_num
    print '\ttotal entries in gene2ortho:', len(gene2ortho)
    print '\ttotal entries in ortho2genes:', len(ortho2genes)
#    print ortho2genes['4258']
#    sys.exit()


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
        locus2species[genbank[s].name] = s
    print '\ttotal records parsed:', len(genbank)
#    print genbank.keys()


    print 'Parsing clusters and assigning genes to them:'
    for s in species:
        clustertrees[s] = IntervalTree()
        cluster2genes[s] = {}
        numbers2products[s] = {}
        coords2numbers[s] = {}
        gene2clusters[s] = {}
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
                # 'cl' is a list of Intervals, each has 'start' and 'end' attributes.
                cl = clustertrees[s].find(int(f.location.start.position), int(f.location.end.position))
                if len(cl) > 0:
                    # Determine which qualifier type to use for gene name.
                    if 'locus_tag' in f.qualifiers:
                        qualifier = 'locus_tag'
                    elif 'gene' in f.qualifiers:
                        qualifier = 'gene'
                    gene_name = f.qualifiers[qualifier][0] + '.' + genbank[s].name
                    # One gene may belong to more than one cluster.
#                    print 'gene at (%s, %s) overlaps with %s cluster(s), 1st is at (%s, %s)' % (f.location.start.position, f.location.end.position, len(cl), cl[0].start, cl[0].end)
                    num_genes += 1
                    gene2clusters[s][gene_name] = []
                    for cluster in cl:
                        cluster2genes[s][coords2numbers[s][(cluster.start, cluster.end)]].append(gene_name)
                        gene2clusters[s][gene_name].append((s, coords2numbers[s][(cluster.start, cluster.end)]))
        print '\t%s: %s clusters populated with %s genes' % (s, len(cluster2genes[s]), num_genes)
    print '\tadded %s clusters from %s species' % (len(all_clusters), len(species))
    # Sort by species.lower() to ensure stable order of cluster pairs.
    all_clusters.sort(key=lambda s: s[0].lower())


    print 'Freeing memory.'
    del genbank
    del clustertrees


    print 'Constructing gene-based cluster links.'
    # Simple counter of the number of cluster pairs we are processing.
    num_pairs = 0
    # Iterate all clusters.
    for c in all_clusters:
        # Iterate each gene in the current cluster.
#        print c,
        if c not in cluster_weights:
            cluster_weights[c] = {}
        s = c[0]
        cluster_number = c[1]
        # Per-gene lists of linked clusters.
        gene_links = {}
        # List of all clusters we have links to. Cannot use set, as clusters come as lists.
        all_links = []
        for g in cluster2genes[s][cluster_number]:
            # Iterate list of all genes from the g's orthology cluster, except for 'g' itself.
            if g not in gene2ortho:
                continue
#            print '\tgene', g
            gene_links[g] = []
            for orthoclust in gene2ortho[g]:
                # Here, *must* create a copy by slicing, as otherwise genes are removed from ortho2genes.
                ortho_genes_wo_g = ortho2genes[orthoclust][:]
                ortho_genes_wo_g.remove(g)
#                print 'ortho_genes_wo_g', ortho_genes_wo_g
                for xeno_gene in ortho_genes_wo_g:
                    # Extract LOCUS from gene name, find species from it.
                    xeno_species = locus2species[xeno_gene.rsplit('.')[-1]]
#                    print 'xeno gene and species:', xeno_gene, xeno_species
                    if xeno_gene in gene2clusters[xeno_species]:
                        # 'cluster' is a list of clusters
                        cluster = gene2clusters[xeno_species][xeno_gene]
                        gene_links[g].extend(cluster)
                        all_links.extend(cluster)
#                        print cluster
#                        print gene_links[g]
        # Uniquefy all_links.
        all_links = set(all_links)
#        print all_links

#        print 'Calculating link weights for', c
        c_genes = len(gene_links)
        # Iterate each linked cluster.
        for link in all_links:
#            print 'link', link
            # Total number of gene-level links between 'c' and 'link'.
            links_between = 0
            # Number of genes in the 'remote' cluster.
            link_genes = len(cluster2genes[link[0]][link[1]])
            for g, v in gene_links.iteritems():
                if link in v:
                    links_between += 1
            # link weight formula:
            # 1. if links_between/min(c_genes, link_genes) == 1.0, then weight is 1.0
            # 2. otherwise, weight = (links_between/2.0) * ( 1/min(c_genes, link_genes) + 1/max(c_genes, link_genes) )
            weight = float(links_between) / max(min(c_genes, link_genes), links_between)
            if weight == 0.0:
                print 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
                continue
            num_pairs += 1
            # If c_genes == link_genes, then recalculation will not change anything.
            if weight < 0.9999 and c_genes != link_genes:
#                print '\tweight', weight
                weight_old = weight
                weight = 0.5 * links_between * ( 1.0 / min(c_genes, link_genes) + 1.0 / max(c_genes, link_genes) )
                try:
                    assert weight < weight_old
                except:
                    print 'old', weight_old, 'new', weight
                    print 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
                    raise
#                print '\tweight', weight
#            if weight >= 1.0:
#                print 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
            if link in cluster_weights[c]:
                cluster_weights[c][link] = max(weight, cluster_weights[c][link])
                print 'this should never happen'
                print 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
                sys.exit()
            else:
                if weight >= threshold:
                    cluster_weights[c][link] = weight
        if len(cluster_weights[c]) == 0:
            del cluster_weights[c]
    print
    print 'Total cluster pair weights calculated: %s' % num_pairs


    print 'Resolving multi-maps to single species by weight, removing self-links, applying threshold.'
    print '%s initial link seeds' % len(cluster_weights)
    for c1 in cluster_weights.iterkeys():
        # If sp1.a->sp2.b=0.5, sp1.a->sp2.c=0.9, then sp1.a->sp2.c.
        # Group all linked clusters by species.
        by_species = {}
        # Populate by-species group.
        for c2 in cluster_weights[c1].iterkeys():
            if c2[0] not in by_species:
                by_species[c2[0]] = []
            by_species[c2[0]].append((cluster_weights[c1][c2], c2))
        # Iterate all groups, filling weights_clean and weights_intra.
#        if c1 == ('scabiei_87.22.faa', 15):
#            print 'by_species', by_species
        for s, v in by_species.iteritems():
            # Intra.
            if s == c1[0]:
                for onec in v:
                    if c1 not in weights_intra:
                        weights_intra[c1] = {}
                    weights_intra[c1][onec[1]] = onec[0]
            # Normal case.
            elif len(v) == 1:
                if c1 not in weights_clean:
                    weights_clean[c1] = {}
                weights_clean[c1][v[0][1]] = v[0][0]
            # Multi-map case.
            elif len(v) > 1:
                best = sorted(v)[-1]
                if c1 not in weights_clean:
                    weights_clean[c1] = {}
                weights_clean[c1][best[1]] = best[0]
#        if c1 == ('scabiei_87.22.faa', 15):
#            print 'intra', weights_intra[c1]
#            print 'clean', weights_clean[c1]
    print '\t%s self-link seeds separated' % len(weights_intra)
    print '\t%s clean non-redundant weight seeds obtained' % len(weights_clean)


    print 'Grouping into clusters with 11, 10, ... links.'
    # Dict grouping weights_clean by the number of links inside.
    by_count = {}
    # Init.
    for i in range(1, len(species) + 1):
        by_count[i] = {}
    for c1, nested_dict in weights_clean.iteritems():
        # +1 for the c1, which is not counted.
        group_len = len(nested_dict) + 1
        try:
            assert group_len <= len(species)
        except:
            print 'group_len', group_len
            print 'c1', c1
            print 'nested', nested_dict
            raise
        by_count[group_len][c1] = nested_dict
    print 'Summary:'
    for i in range(len(species), 0, -1):
        print '\t%s: %s groups' % (i, len(by_count[i]))

    # using cluster products, for each group of same-links clusters output tables for them:
    # Clusters with 11 links (a total of XXX):
    # sp1 sp2 sp3 ...
    # cl1 cl2 cl3 ...
    # ...

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
    parser.add_argument('--threshold', action = 'store', type=float, default = 0.0, help='cluster links with weight below this one will be discarded [default: %(default)s]')

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
    process(paths, args.threshold)
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
