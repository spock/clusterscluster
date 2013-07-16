#!/usr/bin/python
# encoding: utf-8
'''
cluster -- finds biosynthetic clusters linked by orthologs in multiple genomes

cluster is a tool which processes MultiParanoid/QuickParanoid output for
multiple genomes, together with the GenBank files of those genomes, and
calculates "links" between biosynthetic clusters in those genomes. A "link"
between any two clusters exists, if at least one pair of genes in those
clusters belong to the same cluster of orthologs. Each link also has a weight
assigned to it, which tells the fraction of the orthologs between the two clusters.

It defines classes_and_methods

@author:     bogdan

@copyright:  2013 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''


# TODO:
# - calculate "cluster weight" as an average of all link weights
# - add cluster-level threshold (by cluster weight)
# - add an option to ignore putative clusters
# - trim clusters on both ends for the length which antismash2 adds "for insurance"


import sys
import os
import csv


from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from Bio import SeqIO
from bx.intervals.intersection import Interval, IntervalTree
from multiprocessing import Process, Queue, cpu_count


import hmm_detection


__all__ = []
__version__ = 0.3
__date__ = '2013-07-10'
__updated__ = '2013-07-15'


TESTRUN = 0
PROFILE = 0


class Tee(object):
    '''Write to both sys.stdout and provided file/mode, just like shell's "tee" command.'''
    def __init__(self, name, mode):
        self.file = open(name, mode)
        self.stdout = sys.stdout
        sys.stdout = self
    def __del__(self):
        self.close()
    def write(self, data):
        self.file.write(data)
        self.stdout.write(data)
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    def close(self):
        if self.stdout is not None:
            sys.stdout = self.stdout
            self.stdout = None
        if self.file is not None:
            self.file.close()
            self.file = None


def print_cluster_numbers_row(s2c, species, tee):
    first = True
    for s in species:
        if first:
            first = False
        else:
            tee.write('\t')
        if s in s2c:
            tee.write(str(s2c[s]))
    tee.write('\n')


def parse_cluster_number(note):
    'given a list of items from "note" field, return cluster number'
    for i in note:
        if i.startswith('Cluster number: '):
            return int(i[16:])


def print_species_header(species, tee):
    'Print table header.'
    first = True
    for s in species:
        if first:
            first = False
        else:
            tee.write('\t')
        tee.write(s)
    tee.write('\n')


def dedup_clusters(source, target):
    'removes identical clusters from the supplied "source" dict of dicts, writes unique clusters to "target"'
    # List of clusters to skip when creating target.
    skiplist = []
    # FIXME: counter of different weights. Must fix different weights, not count them.
    weights_differ = 0
    total_weight_diffs = 0.0
    min_diff = 1.0
    max_diff = 0.0
    for c1, n in source.iteritems():
        if c1 in skiplist:
            continue
        target[c1] = n
        curclust = [c1]
        curclust.extend(n.keys())
        curclust.sort()
        for sub in n.iterkeys():
            if sub not in source:
                continue
            subflat = [sub]
            subflat.extend(source[sub].keys())
            subflat.sort()
            if subflat == curclust:
                # Check a single weight (e.g. between c1 and sub), just in case
                try:
                    assert source[c1][sub] == source[sub][c1]
                except:
                    weights_differ += 1
                    delta = abs(source[sub][c1] - source[c1][sub])
                    total_weight_diffs += delta
                    min_diff = min(min_diff, delta)
                    max_diff = max(max_diff, delta)
                    # FIXME
#                    print 'c1', c1
#                    print 'source', source[c1]
#                    print 'sub', sub
#                    print 'duplicate', source[sub]
#                    raise
                skiplist.append(sub)
    print '\tfound and removed', len(set(skiplist)), 'duplicate clusters;'
    print '\tnow only', len(target), 'seed clusters remain.'
    print '\tfound %s differing weights; min/avg/max differences are %s, %s, %s.' % (weights_differ, min_diff, total_weight_diffs / weights_differ, max_diff)


def process(config, paranoid, paths, threshold = 0.0, prefix = 'out_', trim = True, skipp = False):
    '''Main method which does all the work.
    "paranoid" is the path to multi/quick-paranoid output file.
    "paths" is a list of paths to genbank files we want to compare.
    "threshold" is a biocluster-biocluster link weight threshold.
    "prefix" is prepended to all output files.
    "trim", if True, causes antismash2 clusters to lose non-core extensions at
    both ends of the cluster (these are hard-coded in antismash2). Note that for composite-type
    clusters only the shorter of the extensions will be trimmed (e.g. bacteriocin extension
    for the backteriocin-t1pks cluster).
    "skipp" means "skip putative clusters", if set.'''

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
    # As above, but after removing identical clusters of clusters.
    weights_dedup = {}
    # Same as cluster_weights, but only for intra-species links.
    weights_intra = {}
    weights_intra_dedup = {}
    # Mapping of record.names from GenBank files (LOCUS) to species.
    locus2species = {}


    if verbose > 0:
        print 'Reading quick/multi-paranoid gene clusters:'
    with open(paranoid) as tsv:
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
                # Build gene-to-ortho-cluster-ID(s) dict. row['gene'] is the gene ID.
                if row['gene'] not in gene2ortho:
                    gene2ortho[row['gene']] = []
                gene2ortho[row['gene']].append(row[idname])
                # Build cluster-ID-to-all-genes dict of lists. row[idname] is the cluster number.
                if row[idname] not in ortho2genes:
                    ortho2genes[row[idname]] = []
                ortho2genes[row[idname]].append(row['gene'])
        except csv.Error as e:
            sys.exit('file %s, line %d: %s' % (paranoid, reader.line_num, e))
    if verbose > 0:
        print '\ttotal lines read:', reader.line_num
        print '\ttotal entries in gene2ortho:', len(gene2ortho)
        print '\ttotal entries in ortho2genes:', len(ortho2genes)


    if verbose > 0:
        print 'Reading species list file:'
    for s in open(config):
        species.append(s.strip())
    species.sort(key = lambda s: s.lower())
    if verbose > 0:
        print '\ttotal species:', len(species)


    if verbose > 0:
        print 'Reading all genbank files in parallel:'
    task_queue = Queue()
    done_queue = Queue()

    # Simple helper
    def worker(tasks, done):
        while not tasks.empty():
            (s, gb) = tasks.get()
            done.put((s, SeqIO.read(gb, "genbank")))

    # Detect which species it is by the first 5 characters.
    # FIXME: needs a better solution.
    for gb in paths:
        for s in species:
            if s[0:5] == gb[0:5]:
                if verbose > 1:
                    print '\t%s corresponds to species %s' % (gb, s)
                break
        task_queue.put((s, gb))
    for i in range(cpu_count()):
        Process(target=worker, args=(task_queue, done_queue)).start()
    for i in range(len(paths)):
        s, rec = done_queue.get()
        genbank[s] = rec
        locus2species[genbank[s].name] = s
    if verbose > 0:
        print '\ttotal records parsed:', len(genbank)
    task_queue.close()
    done_queue.close()


    if verbose > 0:
        print 'Parsing clusters and assigning genes to them:'
    if verbose > 1:
        print '\tgetting extension sizes for diff. cluster types from antismash2 config'
        rulesdict = hmm_detection.create_rules_dict()
    for s in species:
        if verbose > 2:
            print '\tprocessing', s
        clustertrees[s] = IntervalTree()
        cluster2genes[s] = {}
        numbers2products[s] = {}
        coords2numbers[s] = {}
        gene2clusters[s] = {}
        # Populate clusters tree and dict with (start, end) as keys.
        if trim and verbose > 1:
            print '\tNote: cluster coordinates are shown after trimming, except for skipped putative clusters.'
        for f in genbank[s].features:
            if f.type == 'cluster':
                cluster_number = parse_cluster_number(f.qualifiers['note'])
                start = int(f.location.start.position)
                end = int(f.location.end.position)
                if skipp and f.qualifiers['product'][0] == 'putative':
                    if verbose > 1:
                        print '\tskipping putative cluster #%s at (%s, %s)' % \
                                (cluster_number, start, end)
                    continue
                if trim:
                    # Use cluster type to get extension size, including composite types.
                    extension = hmm_detection.get_extension_size_by_cluster_type(f.qualifiers['product'][0], rulesdict)
                    # If cluster is at the genome edge - skip extension trimming.
                    if start == 0:
                        if verbose > 1:
                            print '\tnot trimming left-extension for #%s (%s) - at the genome start' % \
                                   (parse_cluster_number(f.qualifiers['note']), f.qualifiers['product'][0])
                    else:
                        start += extension
                    if end == len(genbank[s]):
                        if verbose > 1:
                            print '\tnot trimming right-extension for #%s (%s) - at the genome end' % \
                                   (parse_cluster_number(f.qualifiers['note']), f.qualifiers['product'][0])
                    else:
                        end -= extension
                    try:
                        assert start < end
                    except:
                        print 'trimming extension failed for: '
                        print f.qualifiers['product'][0], '(%s)' % parse_cluster_number(f.qualifiers['note'])
                        print 'extension', extension
                        print 'ori start %s, new start %s' % (f.location.start.position, start)
                        print 'ori  end %s, new  end %s' % (f.location.end.position, end)
                        raise
                clustertrees[s].add_interval(Interval(start, end))
                cluster2genes[s][cluster_number] = []
                all_clusters.append((s, cluster_number))
                numbers2products[s][cluster_number] = f.qualifiers['product'][0]
                coords2numbers[s][(start, end)] = cluster_number
                if verbose > 1:
                    print '\tadding cluster %s (%s) at (%s, %s)' % (cluster_number, f.qualifiers['product'][0], start, end)
        if verbose > 1:
            print '\tNow assigning genes to biosynthetic clusters'
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
                    if verbose > 2:
                        print '\t\tgene at (%s, %s) overlaps with %s cluster(s), 1st is at (%s, %s)' % (f.location.start.position, f.location.end.position, len(cl), cl[0].start, cl[0].end)
                    num_genes += 1
                    gene2clusters[s][gene_name] = []
                    for cluster in cl:
                        cluster2genes[s][coords2numbers[s][(cluster.start, cluster.end)]].append(gene_name)
                        gene2clusters[s][gene_name].append((s, coords2numbers[s][(cluster.start, cluster.end)]))
        if verbose > 1:
            print '\t%s: %s clusters populated with %s genes' % (s, len(cluster2genes[s]), num_genes)
    # Extra verification.
    for s in species:
        for num, cl in cluster2genes[s].iteritems():
            try:
                assert len(cl) > 0
            except:
                print 'cluster %s of %s has 0 genes' % (num, s)
                raise
    if verbose > 0:
        print '\tadded %s clusters from %s species' % (len(all_clusters), len(species))
    # Sort by species.lower() to ensure stable order of cluster pairs.
    all_clusters.sort(key = lambda s: s[0].lower())
    sys.exit()

    if verbose > 2:
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
            # Iterate list of all genes from the g's orthology cluster, *including* 'g' itself - w/o it we lose 1 link.
            if g not in gene2ortho:
                continue
#            print '\tgene', g
            gene_links[g] = []
            for orthoclust in gene2ortho[g]:
                # Here, *must* create a copy by slicing, as otherwise genes are removed from ortho2genes.
                ortho_genes = ortho2genes[orthoclust][:]
#                ortho_genes.remove(g)
#                print 'ortho_genes_wo_g', ortho_genes_wo_g
                for xeno_gene in ortho_genes:
                    # Extract LOCUS from gene name, find species from it.
                    xeno_species = locus2species[xeno_gene.rsplit('.')[-1]]
#                    print 'xeno gene and species:', xeno_gene, xeno_species
                    if xeno_gene in gene2clusters[xeno_species]:
                        # 'cluster' is a list of clusters
                        cluster = gene2clusters[xeno_species][xeno_gene]
                        for clust in cluster:
                            if clust != c:
                                gene_links[g].append(clust)
                                all_links.append(clust)
        # Uniquefy all_links.
        all_links = set(all_links)


#        print 'Calculating link weights for', c
        c_genes = len(cluster2genes[c[0]][c[1]])
        # Iterate each linked cluster.
        for link in all_links:
#            print 'link', link
            # Total number of gene-level links between 'c' and 'link'.
            links_between = 0
            # Number of genes in the 'remote' cluster.
            link_genes = len(cluster2genes[link[0]][link[1]])
            for v in gene_links.itervalues():
                if link in v:
                    links_between += 1
            # link weight formula:
            # 1. if links_between/min(c_genes, link_genes) == 1.0, then weight is 1.0
            # 2. otherwise, weight = (links_between/2.0) * ( 1/min(c_genes, link_genes) + 1/max(c_genes, link_genes) )
            weight = float(links_between) / max(min(c_genes, link_genes), links_between)
            #custom = [('avermitilis_MA4680.faa', 21), ('cattleya_DSM46488.faa', 42)]
#            custom = [('coelicolor_A3_2.faa', 34), ('SirexAA_E.faa', 34)]
#            if (link in custom and c in custom):
#                print 'c', c, 'link', link
#                print 'c_genes (c)', c_genes, cluster2genes[c[0]][c[1]]
#                print 'link_genes (link)', link_genes, cluster2genes[link[0]][link[1]]
#                print 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
            if weight == 0.0:
                print 'c', c, 'link', link, 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
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
            try:
                assert weight <= 1.0
            except:
                print 'weight', weight
                print 'links_between', links_between, 'c_genes', c_genes, 'link_genes', link_genes
                raise
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


    print 'Resolving multi-maps to single species by weight, and removing links to self.'
    print '\t%s initial link seeds' % len(cluster_weights)
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
    print '\t%s self-link seeds separated' % len(weights_intra)
    print '\t%s clean non-redundant weight seeds obtained' % len(weights_clean)


    print 'De-duplicating inter-species clusters:'
    dedup_clusters(source = weights_clean, target = weights_dedup)
    print 'De-duplicating intra-species clusters:'
    dedup_clusters(source = weights_intra, target = weights_intra_dedup)


    print 'Grouping into clusters with 11, 10, ... links.'
    # Dict grouping weights_dedup by the number of links inside.
    by_count = {}
    # Init.
    for i in range(1, len(species) + 1):
        by_count[i] = {}
    for c1, nested_dict in weights_dedup.iteritems():
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


    #del weights_clean, weights_dedup, weights_intra


    print 'Summary:'
    for i in range(len(species), 0, -1):
        print '\t%s: %s groups' % (i, len(by_count[i]))

    # All output tables have headers, row names, and tab-separated columns.

    # Using cluster products, for each group of same-links-count clusters output tables for them:
    # Clusters with 11 links (a total of XXX):
    # sp1 sp2 sp3 ...
    # cl1 cl2 cl3 ...
    # ...
    for i in range(len(species), 0, -1):
        if len(by_count[i]) == 0:
            continue
        fname = prefix + '_' + str(i) + '_links.csv'
        print 'Groups of size', i, '(%s)' % len(by_count[i]), 'will be written to file', fname
        # From here on output goes to both stdout and the file.
        #tee = Tee(fname, 'a')
        tee = open(fname, 'w')
        print_species_header(species, tee)
        # Print data rows.
        group_s2c = {}
        for c1, nested in by_count[i].iteritems():
            # Map species of the current group to cluster numbers.
            s2c = {}
            s2c[c1[0]] = c1[1]
            for (s, n) in nested.iterkeys():
                s2c[s] = n
            first = True
            for s in species:
                if first:
                    first = False
                else:
                    tee.write('\t')
                if s in s2c:
                    tee.write(numbers2products[s][s2c[s]] + ' (' + str(s2c[s]) + ')')
            tee.write('\n')
            group_s2c[c1] = s2c
        # After each such table, output a diagonal matrix of link weights between all
        # possible cluster pairs in each of the groups above.
        # Cluster link weights:
        #    cl1 cl2 cl3 cl4
        # cl1 -  0.5 0.6 0.9
        # cl2     -  0.7 0.5
        # cl3         -  0.8
        # cl4             -
        tee.write('\nCluster pairs link weights\n')
        for c1, nested in by_count[i].iteritems():
            print_species_header(species, tee)
            print_cluster_numbers_row(group_s2c[c1], species, tee)
            # Keep track of numbers already shown in the 1st column.
            printed = []
            for s_r in species[:]:
                printed.append(s_r)
                if s_r not in group_s2c[c1]:
                    # Do not print entirely empty rows.
                    continue
                # Print cluster number in the 1st column.
                tee.write(str(group_s2c[c1][s_r]))
                first = True
                for s_c in species[:]:
                    if first:
                        # Skip firstmost record entirely.
                        first = False
                        continue
                    tee.write('\t')
                    if s_c in printed:
                        # Don't print the lower diagonal of values.
                        continue
                    # Find an existing link.
                    row_key = (s_r, group_s2c[c1][s_r])
                    col_key = (s_c, group_s2c[c1][s_c])
                    if row_key in by_count[i]:
                        tee.write(round(by_count[i][row_key][col_key], 2))
                    elif col_key in by_count[i]:
                        tee.write(round(by_count[i][col_key][row_key], 2))
                    else:
                        print 'Neither', row_key, 'nor', col_key, 'found in cluster_weights.'
                tee.write('\n')
        tee.close()
#        del tee

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

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-v", "--verbose", dest="verbose", action="count", default=0, help="set verbosity level, up -vvv [default: %(default)s]")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true", default=False, help="set verbosity level to maximal [default: %(default)s]")
    parser.add_argument('-V', '--version', action='version', version=program_version_message)
    parser.add_argument("--no-trim", dest="no_trim", action="store_true", default=False, help="do not trim away antismash2 cluster extensions [default: %(default)s]")
    parser.add_argument("--skip-putative", dest="skipp", action="store_true", default=False, help="exclude putative clusters from the analysis [default: %(default)s]")
    parser.add_argument("--prefix", default='out', help="output CSV files prefix [default: %(default)s]")
    parser.add_argument('--threshold', action = 'store', type=float, default = 0.0, help='cluster links with weight below this one will be discarded [default: %(default)s]')
    parser.add_argument(dest="config", help="path to plain-text species list file", metavar="config")
    parser.add_argument(dest="paranoid", help="path multiparanoid/quickparanoid sqltable file", metavar="sqltable")
    parser.add_argument(dest="paths", help="paths to GenBank files annotated with antismash2", metavar="path", nargs='+')

    # Process arguments
    args = parser.parse_args()
    trim = not args.no_trim

    global verbose
    if args.debug:
        verbose = 3
    else:
        verbose = args.verbose
#    recurse = args.recurse
#    inpat = args.include
#    expat = args.exclude
#    if inpat and expat and inpat == expat:
#        raise CLIError("include and exclude pattern are equal! Nothing will be processed.")
#    for inpath in paths:
#        ### do something with inpath ###
#        print(inpath)
    process(prefix=args.prefix, config=args.config, paranoid=args.paranoid,
            paths=args.paths, threshold=args.threshold, trim=trim, skipp=args.skipp)
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
