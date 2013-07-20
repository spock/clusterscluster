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


import sys
import os
import csv
import itertools
from pprint import pprint

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


def bin_key(weight):
    'for a given weight, return appropriate weight_bin key'
    if weight <= 0.05: return 0.05
    if weight <= 0.10: return 0.1
    if weight <= 0.15: return 0.15
    if weight <= 0.20: return 0.2
    if weight <= 0.25: return 0.25
    if weight <= 0.30: return 0.3
    if weight <= 0.35: return 0.35
    if weight <= 0.40: return 0.4
    if weight <= 0.45: return 0.45
    if weight <= 0.50: return 0.5
    if weight <= 0.55: return 0.55
    if weight <= 0.60: return 0.6
    if weight <= 0.65: return 0.65
    if weight <= 0.70: return 0.7
    if weight <= 0.75: return 0.75
    if weight <= 0.80: return 0.8
    if weight <= 0.85: return 0.85
    if weight <= 0.90: return 0.9
    if weight <= 0.95: return 0.95
    if weight <= 1.00: return 1.0


def process(config, paranoid, paths, threshold = 0.0, prefix = 'out_', trim = True,
            skipp = False, strict = False, ortho = False, sizes = False, names = False):
    '''Main method which does all the work.
    "paranoid" is the path to multi/quick-paranoid output file.
    "paths" is a list of paths to genbank files we want to compare.
    "threshold" is a biocluster-biocluster link weight threshold.
    "prefix" is prepended to all output files.
    "trim", if True, causes antismash2 clusters to lose non-core extensions at
    both ends of the cluster (these are hard-coded in antismash2). Note that for composite-type
    clusters only the shorter of the extensions will be trimmed (e.g. bacteriocin extension
    for the backteriocin-t1pks cluster).
    "skipp" means "skip putative clusters", if set.
    "ortho", if True, will only use orthology clusters which do not have'''

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
    # Two lists of 2-tuples with weight > threshold links between clusters.
    intra_one = []
    inter_one = []
    # Sizes of clusters, in basepairs. clustersizes[(species, number)] = 65450
    clustersizes = {}


    def get_gene_links_to_bioclusters(gene):
        if gene not in gene2ortho:
            return []
        bioclusters = []
        for orthoclust in gene2ortho[gene]:
            if verbose > 3:
                print 'gene', gene, 'belongs to ortho-cluster', orthoclust
            for xeno_gene in ortho2genes[orthoclust]:
                # Extract LOCUS from gene name, find species from it.
                xeno_species = locus2species[xeno_gene.rsplit('.')[-1]]
                # Check if xeno_gene belongs to any biosynthetic clusters.
                if xeno_gene in gene2clusters[xeno_species]:
                    bioclusters.extend(gene2clusters[xeno_species][xeno_gene])
                    if verbose > 3:
                        print '\txeno_gene', xeno_gene, 'belongs to', gene2clusters[xeno_species][xeno_gene]
        return bioclusters


    def calculate_links(c1, c2):
        links = 0
        for gene in cluster2genes[c1[0]][c1[1]]:
            if c2 in get_gene_links_to_bioclusters(gene):
                links += 1
        return links


    def calculate_weight(c1, c2):
        c1_genes = len(cluster2genes[c1[0]][c1[1]])
        c2_genes = len(cluster2genes[c2[0]][c2[1]])
        links1 = calculate_links(c1, c2)
        links2 = calculate_links(c2, c1)
        if verbose > 3:
            print '\tlinks1 =', links1, 'and links2 =', links2, 'for', c1, c2, '(%s and %s genes)' % (c1_genes, c2_genes)
        links = max(links1, links2)
        if links > max(c1_genes, c2_genes):
            links = max(c1_genes, c2_genes)
        if strict:
            if links > c1_genes or links > c2_genes:
                links = min(c1_genes, c2_genes)
                if verbose > 3:
                    print 'strict mode, new links is', links
        try:
            assert links <= c1_genes or links <= c2_genes
        except:
            if verbose > 2:
                print 'c1', c1, 'c2', c2
                print 'c1_genes (c1)', c1_genes, cluster2genes[c1[0]][c1[1]]
                print 'c2_genes (c2)', c2_genes, cluster2genes[c2[0]][c2[1]]
                print 'links', links, 'c1_genes', c1_genes, 'c2_genes', c2_genes
            raise
        # Link weight formula.
        # Multiple versions were tried/examined.
        # 1. weight = links/min(c1_genes, c2_genes)
        # Unfortunately, this one causes 1-gene clusters to have too many weight-1.0 links.
        # 2. weight = (links/2.0) * ( 1/min(c1_genes, c2_genes) + 1/max(c1_genes, c2_genes) )
        # (or the equivalent weight = 0.5 * links * ( 1.0 / c1_genes + 1.0 / c2_genes ) )
        # was giving values > 1.0, as 'links' can be > than c1 or c2.
        # 3. "safe" derivative of #2,
        # weight = 0.5* ( min(c1_genes, links) / c1_genes + min(c2_genes, links) / c2_genes )
        # was too optimistic: 0.53 for links = 1, c1 = 1, c2 = 20.
        # The least biased formula is below. It properly penalizes if the c1_genes and
        # c2_genes are too different, and is also safe against links > c1_genes or links > c2_genes.
        both = float(c1_genes + c2_genes)
        if sizes:
            # Use relative cluster sizes as link contribution weight.
            size1 = clustersizes[c1]
            size2 = clustersizes[c2]
            total_size = float(size1 + size2)
            weight = (size1/total_size) * min(c1_genes, links)/both + (size2/total_size) * min(c2_genes, links)/both
            if verbose > 1:
               print 'weight, ratio1, raw weight 1, ratio2, raw weight 2:', weight, size1/total_size, min(c1_genes, links)/both, size2/total_size, min(c2_genes, links)/both
        else:
            weight = min(c1_genes, links)/both + min(c2_genes, links)/both
        if weight == 0.0:
            if verbose > 3:
                print 'c1', c1, 'c2', c2, 'links', links, 'c1_genes', c1_genes, 'c2_genes', c2_genes
            return weight
        try:
            assert weight <= 1.0
        except:
            print 'weight', weight
            print 'links', links, 'c1_genes', c1_genes, 'c2_genes', c2_genes
            raise
        if weight >= threshold:
            if c1[0] == c2[0]:
                intra_one.append((round(weight, 2), c1, c2))
            else:
                inter_one.append((round(weight, 2), c1, c2))
        return weight


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
                if (ortho and row['tree_conflict'] != 'No') or (names and row['tree_conflict'] == 'diff. names'):
                    continue
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
            try:
                # _nowait: otherwise race conditions emerge, when tasks.empty() evaluates
                # to False, but then some other worker grabs the work and leaves tasks empty.
                (s, gb) = tasks.get_nowait()
            except: # Queue.Empty
                return
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
    if trim and verbose > 1:
        print '\tNote: cluster coordinates are shown after trimming, except for skipped putative clusters.'
    for s in species:
        if verbose > 2:
            print '\tprocessing', s
        clustertrees[s] = IntervalTree()
        cluster2genes[s] = {}
        numbers2products[s] = {}
        coords2numbers[s] = {}
        gene2clusters[s] = {}
        # Populate clusters tree and dict with (start, end) as keys.
        for f in genbank[s].features:
            if f.type == 'cluster':
                cluster_number = parse_cluster_number(f.qualifiers['note'])
                start = int(f.location.start.position)
                end = int(f.location.end.position)
                if skipp and f.qualifiers['product'][0] == 'putative':
                    if verbose > 2:
                        print '\tskipping putative cluster #%s at (%s, %s)' % \
                                (cluster_number, start, end)
                    continue
                # Putative clusters have neither extensions nor rules for them.
                if trim and f.qualifiers['product'][0] != 'putative':
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
                clustersizes[(s, cluster_number)] = end - start
                if verbose > 1:
                    print '\t(%s, %s): %s [%s, %s], %s bp long' % (s, cluster_number, f.qualifiers['product'][0], start, end, end-start)
        if verbose > 2:
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


    if verbose > 2:
        print 'Freeing memory.'
    del genbank
    del clustertrees


    if verbose > 0:
        print 'Constructing gene-based cluster links.'
    pairs = itertools.combinations(all_clusters, 2)
    # Simple counters of the number of cluster pairs/weights we are processing.
    num_pairs = 0
    # Storage for bins to show weights distribution.
    weight_bins = {}
    for x in range(5, 105, 5):
        # 20 bins, step 0.05; key means "less than", e.g. bin '100' has values greater
        # than 95, and smaller than or equal to 100.
        weight_bins[x/100.0] = 0
    for (c1, c2) in pairs:
        num_pairs += 1
        weight = calculate_weight(c1, c2)
        if verbose > 3:
            print '\tassigned weight', weight, 'to', c1, 'and', c2
        elif verbose > 2 and weight > 0.0:
            print '\tassigned weight', weight, 'to', c1, 'and', c2
        weight_bins[bin_key(weight)] += 1
        if weight >= threshold:
            if c1 not in cluster_weights:
                cluster_weights[c1] = {}
            if c2 not in cluster_weights:
                cluster_weights[c2] = {}
            cluster_weights[c1][c2] = weight
            cluster_weights[c2][c1] = weight
    if verbose > 1:
        print '\tGenerated %s pairs between %s clusters' % (num_pairs, len(all_clusters))
    if verbose > 0:
        print 'Distribution of %s weights (bins of size 0.05)' % num_pairs
        print 'bin\tweights\tgraph'
        # Maximal 'width' of 1 bar.
        height = 120
        for x in range(5, 105, 5):
            x = x/100.0
            val = weight_bins[x]
            bar = '#' * int(round(height*val/num_pairs))
            print '%s\t%s\t%s' %(x, val, bar)


    # Manual data examination had shown that many multi-links between species
    # have high scores of up 1.0. Removing them may adversely affect the
    # estimation of new clusters potential. This is why multi-link resolution
    # is commented out below.
    if verbose > 0:
#        print 'Resolving multi-maps to single species by weight, and removing intra-species links.'
        print 'Removing intra-species links.'
    # List of cluster pairs we had already iterated, to avoid double-processing.
    skip_list = []
    # Counters for actual remaining cluster link pairs.
    num_pairs = 0
    num_pairs_intra = 0
    for c1 in cluster_weights.iterkeys():
        # Group all linked clusters by species.
        by_species = {}
        # Populate by-species group.
        for c2 in cluster_weights[c1].iterkeys():
            if (c1, c2) in skip_list:
                continue
            skip_list.append((c1, c2))
            skip_list.append((c2, c1))
            if c2[0] not in by_species:
                by_species[c2[0]] = []
            by_species[c2[0]].append((cluster_weights[c1][c2], c2))
        # Iterate all groups, filling weights_clean and weights_intra.
        for s, v in by_species.iteritems():
            if s == c1[0]: # Intra case: link(s) to own cluster(s).
                if verbose > 3:
                    print 'intra-species (c1, c2, list):', c1, c2, v
                if c1 not in weights_intra:
                    weights_intra[c1] = {}
                for onec in v:
                    if onec[1] not in weights_intra:
                        weights_intra[onec[1]] = {}
                    weights_intra[c1][onec[1]] = onec[0]
                    weights_intra[onec[1]][c1] = onec[0] # mirror
                num_pairs_intra += len(v)
            elif len(v) == 1: # Normal case: 1 linked cluster in other species.
                if c1 not in weights_clean:
                    weights_clean[c1] = {}
                if v[0][1] not in weights_clean:
                    weights_clean[v[0][1]] = {}
                if verbose > 3:
                    print 'one-to-one:', c1, v[0]
                weights_clean[c1][v[0][1]] = v[0][0]
                weights_clean[v[0][1]][c1] = v[0][0] # mirror
                num_pairs += 1
            elif len(v) > 1: # Multi-map case.
                # If sp1.a->sp2.b=0.5, sp1.a->sp2.c=0.9, then sp1.a->sp2.c.
                # See detailed comment above for the reason of commenting out.
#                best = sorted(v)[-1]
#                weights_clean[c1][best[1]] = best[0]
                if c1 not in weights_clean:
                    weights_clean[c1] = {}
                if verbose > 3:
                    print 'one-to-many:', c1, v
#                    print 'best:', best
                num_pairs += len(v)
                for onec in v:
                    if onec[1] not in weights_clean:
                        weights_clean[onec[1]] = {}
                    weights_clean[c1][onec[1]] = onec[0]
                    weights_clean[onec[1]][c1] = onec[0] # mirror
    if verbose > 0:
        print '\t%s unique intra-species links' % num_pairs_intra
    if verbose > 0:
        print '\t%s unique inter-species links' % num_pairs
    del skip_list


    def clusters_of_clusters():
        'not used at the moment, possibly not fully functional - not reviewed'
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


    if verbose > 1:
        print 'All pairs of clusters with link weight over %s (inter-species).' % threshold
        pprint(inter_one)
        print
        print 'All pairs of clusters with link weight over %s (intra-species).' % threshold
        pprint(intra_one)
        print


    def get_unique_clusters(s, allowed_species = False):
        '''returns the quantity of unique clusters in the provided species "s".
        Linked clusters must be from "allowed_species", if specified;
        otherwise, all species are searched for linked clusters.'''
        # Counter of non-unique clusters.
        non_unique = 0
        for u in cluster2genes[s]:
            if (s, u) in weights_clean:
                links_to = weights_clean[(s, u)].keys()
                for link in links_to:
                    if link[0] == s:
                        continue
                    if allowed_species and link[0] not in allowed_species:
                        continue
                    non_unique += 1
                    break
        return len(cluster2genes[s]) - non_unique


    def graph_unique_change_when_adding(reverse = True):
        print 'Graph of the change of unique clusters fraction with each new added genome'
        if not reverse:
            print '(reversed: from genomes with less clusters to genomes with more)'
        # List of tuples (number_of_clusters, species), for sorting.
        numclust_species = []
        for s in species:
            numclust_species.append((len(cluster2genes[s]), s))
        # 'reverse' means DESC
        numclust_species.sort(reverse=reverse)
        # List of the species we are allowed to look for linked clusters in.
        allowed_species = []
        # Horizontal bar "height" (length).
        height = 90
        first = True
        for (numclust, s) in numclust_species:
            if first: # no need to calculate anything, ratio is 1.0
                first = False
                ratio = 1.0
            else:
                unique = get_unique_clusters(s, allowed_species)
                ratio = round(float(unique) / numclust, 2)
            allowed_species.append(s)
            bar = '#' * int(round(height*ratio))
            bar = bar.ljust(height)
            print '%s\t%s\t%s' % (bar, ratio, s)


    graph_unique_change_when_adding()
#    graph_unique_change_when_adding(False)


    def cumulative_growth(reverse = True):
        '''Shows expected and observed growth of the number of unique clusters
        with each new genome'''
        print 'Graph of the ratio of observed/expected total unique clusters.'
        if not reverse:
            print '(reversed: from genomes with less clusters to genomes with more)'
        # List of tuples (number_of_clusters, species), for sorting.
        numclust_species = []
        for s in species:
            numclust_species.append((len(cluster2genes[s]), s))
        # 'reverse' means DESC
        numclust_species.sort(reverse=reverse)
        # All-important counters.
        total_expected = 0
        total_observed = 0
        # List containing all the data for the table, as tuples
        # (observed, expected, ratio, species)
        table =[]
        # Horizontal bar "height" (length).
        height = 90
        # Draw a reference 100% line.
        print '%s\t%s\t%s' % ('#' * height, 1.0, 'Reference')
        for (numclust, s) in numclust_species:
            total_expected += len(cluster2genes[s])
            total_observed += get_unique_clusters(s)
            ratio = round(float(total_observed) / total_expected, 2)
            table.append((total_observed, total_expected, round(ratio, 2), s))
            bar = '#' * int(round(height*ratio))
            bar = bar.ljust(height)
            print '%s\t%s\t%s' % (bar, ratio, s)
        print 'Data table'
        print 'Obs.\tExp.\tRatio\tGenome'
        for item in table:
            print '%s\t%s\t%s\t%s' % item

    cumulative_growth()
#    cumulative_growth(False)


    # Finally, show a list of all species, stating the number of unique
    # clusters they have.
    print 'Graph of the percentage of unique clusters in each of the genomes.'
    height = 60
    # Draw a reference 100% line.
    print '%s\t%s\t%s\t%s\t%s' % ('Bar'.ljust(height), 'Unique', 'Total', 'Ratio', 'Genome')
    print '%s\t%s\t%s\t%s\t%s' % ('#' * height, 1.0, 1.0, 1.0, 'Reference')
    for s in species:
        unique = get_unique_clusters(s)
        total = len(cluster2genes[s])
        ratio = round(unique / float(total), 2)
        bar = '#' * int(round(height*ratio))
        bar = bar.ljust(height)
        print '%s\t%s\t%s\t%s\t%s' % (bar, unique, total, ratio, s)


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
    parser.add_argument("--strict", dest="strict", action="store_true", default=False, help="weight between clusters with 5 and 10 genes will never exceed 0.5 [default: %(default)s]")
    parser.add_argument("--use-sizes", dest="use_sizes", action="store_true", default=False, help="each cluster's contribution to link weight is scaled by relative cluster sizes; can be combined with --strict [default: %(default)s]")
    parser.add_argument("--no-names", dest="no_names", action="store_true", default=False, help="only use ortho-clusters which do not have diff.names problems [default: %(default)s]")
    parser.add_argument("--no-problems", dest="no_problems", action="store_true", default=False, help="only use ortho-clusters which do not have [diff.names/diff.numbers] problems [default: %(default)s]")
    parser.add_argument("--prefix", default='out', help="output CSV files prefix [default: %(default)s]")
    parser.add_argument('--threshold', action = 'store', type=float, default = 0.0, help='cluster links with weight below this one will be discarded [default: %(default)s]')
    parser.add_argument('--height', action = 'store', type=int, default = 50, help='bar heights for text graphs [default: %(default)s]')
    parser.add_argument(dest="config", help="path to plain-text species list file", metavar="config")
    parser.add_argument(dest="paranoid", help="path multiparanoid/quickparanoid sqltable file", metavar="sqltable")
    parser.add_argument(dest="paths", help="paths to GenBank files annotated with antismash2", metavar="path", nargs='+')

    # Process arguments
    args = parser.parse_args()
    trim = not args.no_trim

    global verbose
    global height
    height = args.height
    if args.debug:
        verbose = 10
    else:
        verbose = args.verbose

#    if inpat and expat and inpat == expat:
#        raise CLIError("include and exclude pattern are equal! Nothing will be processed.")

    print 'Used arguments and options:'
    pprint(args)
    process(prefix=args.prefix, config=args.config, paranoid=args.paranoid,
            paths=args.paths, threshold=args.threshold, trim=trim, skipp=args.skipp,
            strict=args.strict, ortho=args.no_problems, sizes=args.use_sizes,
            names=args.no_names)
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
