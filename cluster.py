#!/usr/bin/python
# encoding: utf-8
'''
cluster -- find biosynthetic clusters linked by orthologs in multiple genomes.

First, cluster processes all input genbank files to associate genome sequence
with ID, species and strain information. If, during the analysis, a duplicate ID
is found, warning is issued and a numeric dot-suffix (e.g. ".1") is appended to
the ID. If all components are identical (ID, species, strain, and genome
length/contig count), then an exception is raised.
For each input file, nucleotide sequence is extracted as a [multi-]fasta file,
and annotated using antismash (--cpus 1, but up to as many parallel instances as
there are cores, likely using queues). This step ensures that there are no
differences in the methods of gene finding. Antismash annotation, depending on
the options, will or will not extend identified clusters. Translations of
antismash-annotated CDS are then extracted to a protein multi-fasta file.
At the end of the 1st stage, program maintains a relation between genome
attributes (ID, etc) and 2 files: antismash-generated genbank, and .faa.

Second, when all the input files have been processed as described, an all-vs-all
pblast is performed between all possible genome pairs, and also each genome is
pblasted against itself. InParanoid is applied to all pairs of genomes to find
orthologs.

Third, multiparanoid or quickparanoid is applied to generate a single table of
orthologous groups.

Finally, cluster processes MultiParanoid/QuickParanoid output for
multiple genomes, together with the GenBank files of those genomes, and
calculates "links" between biosynthetic clusters in those genomes. A "link"
between any two clusters exists, if at least one pair of genes in those
clusters belongs to the same cluster of orthologs. Each link also has a weight
assigned to it, which tells the fraction of the orthologs between the two clusters.

cluster outputs CSV files with collected information for further analysis.
'''


from __future__ import print_function
import sys
import logging
import os
import glob
import csv

from pprint import pprint
from os import mkdir, symlink, getcwd, chdir, remove
from os.path import exists, join, dirname, realpath, basename#, splitext
from shutil import move
from argparse import ArgumentParser
from multiprocessing import Process, Queue, cpu_count
from itertools import permutations, combinations, product, combinations_with_replacement

from lib import utils
from lib.ClusterPair import ClusterPair, Cluster
from lib.Genome import Genome
from lib.MultiParanoid import MultiParanoid


__all__ = []
__version__ = 0.6
__date__ = '2013-07-10'
__updated__ = '2014-04-11'


def print_cluster_numbers_row(s2c, species, tee):
    '''
    Print a row of tab-separated cluster numbers.
    s2c: species to clusters dict
    species: list of identifiers
    tee: file-like handle
    '''
    first = True
    for s in species:
        if first:
            first = False
        else:
            tee.write('\t')
        if s in s2c:
            tee.write(str(s2c[s]))
    tee.write('\n')


def print_species_header(species, tee):
    '''
    Print table header: a row of tab-separated genome identifiers
    species: list of identifiers
    tee: file-like handle
    '''
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
    return 1.0

#
# Set of functions used by process()
#


def calculate_weight(c1, c2, cluster2genes, clustersizes, intra_one, inter_one,
                     numbers2products, mp, gene2clusters, args, inputs):
    '''
    Given 2 biosynthetic clusters - c1 and c2 - calculate the weight of the
    link between them.
    '''
    c1_genes = len(cluster2genes[c1[0]][c1[1]])
    c2_genes = len(cluster2genes[c2[0]][c2[1]])
    links1 = float(min(calculate_links(c1, c2, cluster2genes, mp, gene2clusters, inputs), c1_genes))
    links2 = float(min(calculate_links(c2, c1, cluster2genes, mp, gene2clusters, inputs), c2_genes))
    logging.debug('\tlinks1 = %s and links2 = %s for %s and %s (%s and %s genes)',
                  links1, links2, c1, c2, c1_genes, c2_genes)
    if args.strict:
        # No link can be larger than the number of genes in the 2nd cluster.
        if links1 > c2_genes:
            links1 = float(c2_genes)
            logging.debug('strict mode, new links1 is %s', links1)
        if links2 > c1_genes:
            links2 = float(c1_genes)
            logging.debug('strict mode, new links2 is %s', links2)
        try:
            assert links1 <= min(c1_genes, c2_genes) and links2 <= min(c1_genes, c2_genes)
        except:
            logging.exception('c1: %s ; c2: %s', c1, c2)
            logging.exception('c1_genes: %s (c1: %s)', c1_genes, cluster2genes[c1[0]][c1[1]])
            logging.exception('c2_genes: %s (c2: %s)', c2_genes, cluster2genes[c2[0]][c2[1]])
            logging.exception('links1: %s ; links2: %s ; c1_genes: %s ; c2_genes: %s',
                              links1, links2, c1_genes, c2_genes)
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
    size1 = float(clustersizes[c1])
    size2 = float(clustersizes[c2])
    if args.use_sizes:
        # Use relative physical cluster sizes as link contribution weight.
        total_size = float(size1 + size2)
        weight = 1 / total_size * ( size1 * links1 / c1_genes + size2 * links2 / c2_genes )
        if weight > 0:
            logging.debug('\tweight %s\t%s of %s\t+\t%s of %s', round(weight, 2),
                          round(size1/total_size, 2), round(links1/c1_genes, 2),
                          round(size2/total_size, 2), round(links2/c2_genes, 2))
    else:
        both = float(c1_genes + c2_genes)
        weight = links1/both + links2/both
    if args.scale:
        weight2 = weight * min(size1, size2) / max(size1, size2)
        logging.debug('\tscaled weight %s to %s', round(weight, 2), round(weight2, 2))
        weight = weight2
        del weight2
    try:
        assert weight <= 1.0001 # precision allowance
    except:
        logging.exception('weight: %s', weight)
        logging.exception('links1: %s ; links2: %s ; c1_genes: %s ; c2_genes: %s',
                          links1, links2, c1_genes, c2_genes)
        raise
    if weight >= args.threshold:
        if c1[0] == c2[0]: #intra-species link
            intra_one.append((c1, numbers2products[c1[0]][c1[1]],
                              len(cluster2genes[c1[0]][c1[1]]), clustersizes[c1],
                              c2, numbers2products[c2[0]][c2[1]],
                              len(cluster2genes[c2[0]][c2[1]]),
                              clustersizes[c2], round(weight, 2)))
        else:
            inter_one.append((c1, numbers2products[c1[0]][c1[1]],
                              len(cluster2genes[c1[0]][c1[1]]), clustersizes[c1],
                              c2, numbers2products[c2[0]][c2[1]],
                              len(cluster2genes[c2[0]][c2[1]]),
                              clustersizes[c2], round(weight, 2)))
    return weight


def get_unique_clusters(s, cluster2genes, weights_clean, allowed_species = False):
    '''
    Return the quantity of unique clusters in the provided species "s".
    Linked clusters must be from "allowed_species", if specified;
    otherwise, all species are searched for linked clusters.
    '''
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


def graph_unique_change_when_adding(species, cluster2genes, weights_clean, reverse = True):
    print("Graph of the change of unique clusters' fraction with each new added genome.")
    if not reverse:
        print('(reversed: from genomes with less clusters to genomes with more)')
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
        # FIXME: some species have 0 clusters (numclust = 0), foresee this!
        if first: # no need to calculate anything, ratio is 1.0
            first = False
            ratio = 1.0
        else:
            unique = get_unique_clusters(s, cluster2genes, weights_clean, allowed_species)
            ratio = round(float(unique) / numclust, 2)
        allowed_species.append(s)
        bar = '#' * int(round(height*ratio))
        bar = bar.ljust(height)
        print('%s\t%s\t%s' % (bar, ratio, s))

def cumulative_growth(species, cluster2genes, weights_clean, reverse = True):
    '''
    Show expected and observed growth of the number of unique clusters
    with each new added genome.
    '''
    print('Graph of the ratio of observed/expected total unique clusters.')
    if not reverse:
        print('(reversed: from genomes with less clusters to genomes with more)')
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
    print('%s\t%s\t%s' % ('#' * height, 1.0, 'Reference'))
    for (_, s) in numclust_species:
        total_expected += len(cluster2genes[s])
        total_observed += get_unique_clusters(s, cluster2genes, weights_clean)
        ratio = round(float(total_observed) / total_expected, 2)
        table.append((total_observed, total_expected, round(ratio, 2), s))
        bar = '#' * int(round(height*ratio))
        bar = bar.ljust(height)
        print('%s\t%s\t%s' % (bar, ratio, s))
    print('Data table')
    print('Obs.\tExp.\tRatio\tGenome')
    for item in table:
        print('%s\t%s\t%s\t%s' % item)

def clusters_of_clusters(species, weights_clean, numbers2products, args):
    'not used at the moment, possibly not fully functional - not reviewed'
    print('Grouping into clusters with 11, 10, ... links.')
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
            logging.exception('group_len %s', group_len)
            logging.exception('c1 %s', c1)
            logging.exception('nested %s', nested_dict)
            raise
        by_count[group_len][c1] = nested_dict
    print('Summary:')
    for i in range(len(species), 0, -1):
        print('\t%s: %s groups' % (i, len(by_count[i])))

    # All output tables have headers, row names, and tab-separated columns.

    # Using cluster products, for each group of same-links-count clusters output tables for them:
    # Clusters with 11 links (a total of NUMBER):
    # sp1 sp2 sp3 ...
    # cl1 cl2 cl3 ...
    # ...
    for i in range(len(species), 0, -1):
        if len(by_count[i]) == 0:
            continue
        fname = args.prefix + '_' + str(i) + '_links.csv'
        print('Groups of size', i, '(%s)' % len(by_count[i]),
              'will be written to file', fname)
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
                        print('Neither', row_key, 'nor', col_key,
                              'found in cluster_weights.')
                tee.write('\n')
        tee.close()
#        del tee
        # Once per species, show a table of intra-species cluster links, with weights;
        # looks just like the above weights table, but now only clusters from single
        # species are used.

#
# /Set of functions used by process()
#


def process(all_clusters, inputs, paranoid, args):
    # TODO: simplify, split up this function
    '''
    Analyze quickparanoid results. 'args' contains:
    "paranoid" is the path to quickparanoid output file.
    "args.paths" is a list of paths to genbank files we want to compare.
    "args.threshold" is a biocluster-biocluster link weight threshold.
    "args.prefix" is prepended to all output files.
    "args.trim", if True, causes antismash2 clusters to lose non-core extensions at
    both ends of the cluster (these are hard-coded in antismash2). Note that for composite-type
    clusters only the shorter of the extensions will be trimmed (e.g. bacteriocin extension
    for the backteriocin-t1pks cluster).
    "args.skipp" means "skip putative clusters", if set.
    "args.no_tree_problems", if True, will only use orthology clusters which do not have any problems in the tree_conflict column.
    "args.no_name_problems": as above, but only for the 'diff. names' problem in the tree_conflict column.
    "args.use_sizes" will weigh each clusters contribution to final weight according to clusters length proportion of total length, in bp.
    "args.scale" will scale down weight by the ratio of physical cluster lengths: min(size1, size2)/max(size1, size2).
    '''
    # TODO: deprecate skipp? (can be filtered out at the analysis step)

    # Declare important variables.
    # List of recognized species IDs.
    species = inputs.keys()
    # 2-level nested dict of cluster pairs link weights,
    # e.g. cluster_weights['A'] = {'B': 0.95, ...}
    # TODO: replace with numpy matrix/array
    cluster_weights = {}
    # Same as above, but without duplicate links to other species.
    weights_clean = {} # ????
    # Same as cluster_weights, but only for intra-species links.
    weights_intra = {}
    # Two lists of 2-tuples with weight > threshold links between clusters.
    intra_one = []
    inter_one = []

    # Storage for bins to show weights distribution.
    weight_bins = {}
    for x in range(5, 105, 5):
        # 20 bins, step 0.05; key means "less than", e.g. bin '100' has values greater
        # than 95, and smaller than or equal to 100.
        weight_bins[x/100.0] = 0
    for (c1, c2) in pairs:
        num_pairs += 1
        weight = calculate_weight(c1, c2, cluster2genes, clustersizes, intra_one,
                                  inter_one, numbers2products, mp, gene2clusters, args, inputs)
        logging.debug('\tassigned weight %s to %s and %s', round(weight, 2), c1, c2)
        weight_bins[bin_key(weight)] += 1
        if weight >= args.threshold:
            if c1 not in cluster_weights:
                cluster_weights[c1] = {}
            if c2 not in cluster_weights:
                cluster_weights[c2] = {}
            cluster_weights[c1][c2] = weight
            cluster_weights[c2][c1] = weight
    print('\tGenerated %s pairs between %s clusters' % (num_pairs, len(all_clusters)))
    print('Distribution of %s weights (bins of size 0.05)' % num_pairs)
    print('bin\tweights\tgraph')
    # Maximal 'width' of 1 bar.
    height = 120
    for x in range(5, 105, 5):
        x = x/100.0
        val = weight_bins[x]
        bar = '#' * int(round(height*val/num_pairs))
        print('%s\t%s\t%s' % (x, val, bar))


    # Manual data examination had shown that many multi-links between species
    # have high scores of up to 1.0. Removing them may adversely affect the
    # estimation of new clusters potential. This is why multi-link resolution
    # is commented out below.
#    print 'Resolving multi-maps to single species by weight, and removing intra-species links.'
    print('Removing intra-species links.')
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
                logging.debug('intra-species (c1, c2, list): %s, %s, %s', c1, c2, v)
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
                logging.debug('one-to-one: %s, %s', c1, v[0])
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
                logging.debug('one-to-many: %s, %s', c1, v)
#                    print 'best:', best
                num_pairs += len(v)
                for onec in v:
                    if onec[1] not in weights_clean:
                        weights_clean[onec[1]] = {}
                    weights_clean[c1][onec[1]] = onec[0]
                    weights_clean[onec[1]][c1] = onec[0] # mirror
    print('\t%s unique intra-species links' % num_pairs_intra)
    print('\t%s unique inter-species links' % num_pairs)
    del skip_list

    print('All pairs of clusters with link weight over %s (inter-species).' % args.threshold)
#        pprint(inter_one)
    print('Species,\tNumber,\tType,\tGenes,\tLength,\tSpecies,\tNumber,\tType,\tGenes,\tLength,\tWeight')
    for (cl1, t1, g1, l1, cl2, t2, g2, l2, w) in inter_one:
        print('%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s\n' %
              (cl1[0], cl1[1], t1, g1, l1, cl2[0], cl2[1], t2, g2, l2, w), end='')
    print('All pairs of clusters with link weight over %s (intra-species).' % args.threshold)
#        pprint(intra_one)
    print('Species,\tNumber,\tType,\tGenes,\tLength,\tSpecies,\tNumber,\tType,\tGenes,\tLength,\tWeight')
    for (cl1, t1, g1, l1, cl2, t2, g2, l2, w) in intra_one:
        print('%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s,\t%s\n' %
              (cl1[0], cl1[1], t1, g1, l1, cl2[0], cl2[1], t2, g2, l2, w), end='')
    print()

    graph_unique_change_when_adding(species, cluster2genes, weights_clean)
#    graph_unique_change_when_adding(species, cluster2genes, weights_clean, False)

    cumulative_growth(species, cluster2genes, weights_clean)
#    cumulative_growth(species, cluster2genes, weights_clean, False)

    # Finally, show a list of all species, stating the number of unique
    # clusters they have.
    print('Graph of the percentage of unique clusters in each of the genomes.')
    height = 60
    # Draw a reference 100% line.
    print('%s\t%s\t%s\t%s\t%s' % ('Bar'.ljust(height), 'Unique', 'Total', 'Ratio', 'Genome'))
    print('%s\t%s\t%s\t%s\t%s' % ('#' * height, 1.0, 1.0, 1.0, 'Reference'))
    for s in inputs.keys():
        unique = get_unique_clusters(s, cluster2genes, weights_clean)
        total = len(cluster2genes[s])
        ratio = round(unique / float(total), 2)
        bar = '#' * int(round(height*ratio))
        bar = bar.ljust(height)
        print('%s\t%s\t%s\t%s\t%s' % (bar, unique, total, ratio, s))


def preprocess_input_files(inputs, args):
    '''
    inputs: dict[genome.id] = Genome, is populated by this function.
    args.paths: list of genbank files to process.
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
                logging.debug('Append _1, _2 etc to the ID until it becomes unique.')
                for i in range(1, 10):
                    new_id = g.id + '_' + str(i)
                    if new_id in inputs:
                        continue
                    else:
                        g.id = new_id
                        del new_id
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
    del _, p
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
    for k in list(inputs.iterkeys()):
        del inputs[k]
    for _ in range(total_genomes):
        g = done_queue.get()
        if g != None:
            inputs[g.id] = g
        # Populate all_clusters.
        #for c in g.clusters:
        #    all_clusters.append(cluster(g.id, c))
    # 6. Wait for processes to finish, close queues.
    logging.debug("Joining processes.")
    for p in workers_list:
        p.join()
    del p, workers_list, total_genomes
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
    del _, p
    # 3. Populate the tasks queue.
    logging.debug("Populating the queue.")
    counter = 0
    for _ in faafiles:
        counter += 1
        tasks.put((counter, _)) # will block until all-qsize items are consumed
    del _, counter
    # 4. Add STOP messages.
    logging.debug("Adding %s STOP messages to task_queue.", num_workers)
    for _ in range(num_workers):
        tasks.put((0, 'STOP'))
    del _
    # 5. Wait for all processes to finish, close the queue.
    for p in workers:
        p.join()
    tasks.close()
    del p, total_genomes, workers

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
    del _, p
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
    del pair, counter
    # 4. Add STOP messages.
    logging.debug("Adding %s STOP messages to task_queue.", num_workers)
    for _ in range(num_workers):
        tasks.put((0, 'STOP', ''))
    del _
    # 5. Wait for all processes to finish, close the queue.
    for p in workers:
        p.join()
    tasks.close()
    del p, total_permutations, workers

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
        del _, p
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
        del _
        # 5. Wait for all processes to finish, close the queue.
        for p in workers:
            p.join()
        tasks.close()
        del p, total_combinations, workers, num_workers
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


# multiparanoid, reportedly, cannot handle more than ~20 species.
# multiparanoid requires program file editing for input and output paths.
#def run_multiparanoid(inputs, args):
#    '''
#    multiparanoid.pl -species Actinosynnema_mirum_DSM43827.faa+Amycolatopsis_mediterranei_S699.faa+Kitasatospora_setae_KM6054.faa
#    '''
#    pass


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
    del _

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
    parser.add_argument("--cutoff", dest="cutoff", type = float, default = 0.4, help="protein identity cut-off when aligning with usearch [default: %(default)s]")
    parser.add_argument("--skip-putative", dest="skipp", action="store_true", default=False, help="exclude putative clusters from the analysis [default: %(default)s]")
    parser.add_argument("--skip-orthology", action="store_true", default=False, help="do not run any orthology analysis [default: %(default)s]")
    parser.add_argument("--strict", dest="strict", action="store_true", default=False, help="weight between clusters with 5 and 10 genes will never exceed 0.5 [default: %(default)s]")
    parser.add_argument("--scale", dest="scale", action="store_true", default=False, help="scale link weight down by a factor of min(size1, size2)/max(size1, size2) [default: %(default)s]")
    parser.add_argument("--use-sizes", dest="use_sizes", action="store_true", default=False, help="each cluster's contribution to link weight is scaled by relative cluster sizes; can be combined with --strict [default: %(default)s]")
    parser.add_argument("--no-name-problems", dest="no_name_problems", action="store_true", default=False, help="only use ortho-clusters which do not have diff.names tree_conflict problems [default: %(default)s]")
    parser.add_argument("--no-tree-problems", dest="no_tree_problems", action="store_true", default=False, help="only use ortho-clusters which do not have [diff.names/diff.numbers] tree_conflict problems [default: %(default)s]")
    parser.add_argument('--emulate-inparanoid', action = 'store_true', default = False, help='only print generated inparanoid commands, do not run; exit after inparanoid blasting; suppress some normal output [default: %(default)s]')
    parser.add_argument("--prefix", default='out', help="output CSV files prefix [default: %(default)s]")
    parser.add_argument("--project", default='cluster_project', help="put all the project files into this directory [default: %(default)s]")
    parser.add_argument('--force', action = 'store_true', default = False, help='insist on re-using existing project directory (this will re-use existing intermediate files) [default: %(default)s]')
    parser.add_argument('--no-extensions', action = 'store_true', default = False, help='pass --no-extensions option to the modified antismash2 (see README for details) [default: %(default)s]')
    parser.add_argument('--threshold', action = 'store', type=float, default = 0.0, help='cluster links with weight below this one will be discarded [default: %(default)s]')
    parser.add_argument('--height', action = 'store', type=int, default = 50, help='bar heights for text graphs [default: %(default)s]')
    parser.add_argument('--from-file', action = 'store', help='read paths to GenBank files (one per line) from the provided file')
    parser.add_argument(dest="paths", help="paths to the GenBank files with genomes to analyze", metavar="path", nargs='*')
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
                del _
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
            logging.warning('Re-using existing project directory "%s" as requested.' % args.project)
        else:
            logging.error('Specified project directory "%s" already exists! Use --force to continue anyway.' % args.project)
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

    # Counter of submitted cluster pairs.
    submitted_tasks = 0
    # Counter of processed cluster pairs.
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

    # Preparing for parallel ClusterPair processing.
    # 0. Define queues.
    tasks = Queue() # maxsize = cpu_count() * 2
    done  = Queue()
    # 1. Define worker.
    def cp_processor(tasks, done, mp, genomes, args):
        '''Fully process single ClusterPair from 'tasks', put CSV row into 'done'.'''
        while True:
            try:
                task = tasks.get()
            except:
                logging.exception('tasks queue empty in cp_processor?')
                raise
            if task == 'STOP':
#                logging.debug('STOP found, exiting.')
                break
            cl1, cl2, g1, g2 = task
            cp = ClusterPair(cl1, cl2)
            # FIXME: add back support for args.skip_orthology
            # Calculate the number of orthologous links. mp, genomes, args are all taken from context.
            cp.assign_orthologous_link(mp, genomes, args)
            if cp.link1 > 0 or cp.link2 > 0:
                # Calculate gene-level and clusterpair-average protein identities in clusters.
                cp.CDS_identities(g1, g2, args.cutoff, args.fulldp)
                if cp.avg_identity[0] > 0:
#                    logging.debug('Average identity of %s and %s is %s ',
#                                  cp.gc1, cp.gc2, cp.avg_identity)
                    #print(cp.protein_identities)
                    # TODO: Optional, depends on args: end-trim non-similar genes?
                    if cp.avg_identity[0] > 2:
                        # Calculate gene order and orientation (strandedness) preservation (for similar genes).
                        cp.gene_order(g1, g2)
                    # Calculate predicted domains order preservation within similar genes.
                    #cp.domains(genomes)
                    # Calculate cluster-level nucleotide identity.
                    #cp.nucleotide_similarity(genomes)
                    #cluster_pairs.append(cp) # may not need to collect all these cluster pairs, simply dump them and forget
                    row = [int(cp.intra), g1.id, g2.id,
                           g1.species, g2.species, cp.c1, cp.c2,
                           g1.number2products[cp.c1],
                           g2.number2products[cp.c2],
                           cp.num_c1_genes(g1), cp.num_c2_genes(g2),
                           g1.clustersizes[cp.c1],
                           g2.clustersizes[cp.c2], cp.link1, cp.link2,
                           cp.avg_identity[0], round(cp.avg_identity[1], 1),
                           round(cp.pearson, 2), round(cp.kendall, 2),
                           round(cp.spearman, 2)]
                    done.put(row)
#                    logging.debug('row')
                    continue
            done.put(None)
#            logging.debug('None')
    # 2. Start workers.
    workers_list = []
    logging.debug('Starting %s cp_processors.', cpu_count())
    for _ in range(cpu_count()):
        p = Process(target=cp_processor, args=(tasks, done, mp, genomes, args))
        p.start()
        workers_list.append(p)
    del p

    for g1, g2 in combinations_with_replacement(genomes.keys(), r = 2):
        combinations_counter += 1
        print('%s / %s\t' % (combinations_counter, total_combinations), genomes[g1].id, genomes[g2].id)
        if prev_g1 and g1 != prev_g1:
            genomes[prev_g1].unload()
            prev_g1 = g1
        genomes[g1].load()
        genomes[g2].load() # if g1 == g2, this will do nothing (g1 already loaded)
        genome1 = genomes[g1]
        genome2 = genomes[g2]
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
        # 3. Populate tasks.
#        logging.debug('Populating tasks queue.')
        # Genome-pair specific counter of submitted tasks.
        local_submitted = 0
        for cl1, cl2 in pairslist:
            local_submitted += 1
            # Have to pass entire genomes dict, as cp_processor is defined before we .load() any genomes.
            tasks.put((cl1, cl2, genome1, genome2))
#            logging.debug('Local-submitted %s.', local_submitted)
        submitted_tasks += local_submitted
        # 4. Collect results and write them to file.
#        logging.debug('Collecting results.')
        for _ in range(local_submitted):
            row = done.get()
#            logging.debug('Collected task %s of %s local.', _ + 1, local_submitted)
            if row == None: # no result for this cluster pair
                continue
            cluster_pairs_counter += 1
            writer.writerow(row)
#        logging.debug('Done collecting results.')

        if g1 != g2:
            genomes[g2].unload()
        #del genome1, genome2

    # 5. Add STOP messages.
    logging.debug('Adding %s STOP messages to tasks.', cpu_count())
    for _ in range(cpu_count()):
        tasks.put('STOP')
    # 6. Wait for processes to finish, close queues.
    logging.debug('Joining processes.')
    for p in workers_list:
        p.join(timeout = 10)
    del p, workers_list
    tasks.close()
    done.close()

    print('Processed %s of %s cluster pairs.' % (cluster_pairs_counter,
                                                 submitted_tasks))
    del submitted_tasks, cluster_pairs_counter
    csvout.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
