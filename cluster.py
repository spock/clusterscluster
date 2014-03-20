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

from pprint import pprint
from os import mkdir, rename, symlink, getcwd, chdir, remove
from os.path import exists, join, dirname, realpath, basename#, splitext
from shutil import rmtree
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Alphabet import generic_dna#, generic_protein
from bx.intervals.intersection import Interval, IntervalTree
from multiprocessing import Process, Queue, cpu_count
from itertools import permutations, combinations

from lib import MultiParanoid as MP
from lib.gb2fasta import gb2fasta
from lib import utils
from lib.extract_translation_from_genbank import extract_translation_from_genbank
# TODO: extension trimming will be implemented by modifying antismash2
# configuration, obsoleting this import.
from lib import hmm_detection


__all__ = []
__version__ = 0.4
__date__ = '2013-07-10'
__updated__ = '2014-03-20'


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
    '''Given a list of items from the "note" field of the GenBank feature,
    created by Antismash2, return cluster number.'''
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
    return 1.0

#
# Set of functions used by process()
#

def get_gene_links_to_bioclusters(gene, mp, gene2clusters, inputs):
    '''For the given gene,
    - find to which orthology group/cluster it belongs,
    - check if other genes from that group are a part of some biosynthetic clusters,
    - return the list of those biosynthetic clusters'''
    if gene not in mp.gene2ortho:
        return []
    bioclusters = []
    for orthoclust in mp.gene2ortho[gene]:
        logging.debug('gene %s belongs to ortho-cluster %s', gene, orthoclust)
        for xeno_gene in mp.ortho2genes[orthoclust]:
            # Extract LOCUS from gene name, find species from it.
            #xeno_species = inputs[xeno_gene.rsplit('.')[-1]]
            xeno_species = xeno_gene.rsplit(':')[-1]
            # Check if xeno_gene belongs to any biosynthetic clusters.
            if xeno_gene in gene2clusters[xeno_species]:
                bioclusters.extend(gene2clusters[xeno_species][xeno_gene])
                logging.debug('\txeno_gene %s belongs to %s', xeno_gene,
                              gene2clusters[xeno_species][xeno_gene])
    return bioclusters


def calculate_links(c1, c2, cluster2genes, mp, gene2clusters, inputs):
    '''Given biosynthetic clusters c1 and c2, count the number of unique
    orthologoues gene pairs ("links") between them.'''
    links = 0
    for gene in cluster2genes[c1[0]][c1[1]]:
        if c2 in get_gene_links_to_bioclusters(gene, mp, gene2clusters, inputs):
            links += 1
    return links


def calculate_weight(c1, c2, cluster2genes, clustersizes, intra_one, inter_one,
                     numbers2products, mp, gene2clusters, args, inputs):
    '''Given 2 biosynthetic clusters - c1 and c2 - calculate the weight of the
    link between them.'''
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
    try:
        assert weight <= 1.0001 # precision allowance
    except:
        logging.exception('weight: %s', weight)
        logging.exception('links1: %s ; links2: %s ; c1_genes: %s ; c2_genes: %s',
                          links1, links2, c1_genes, c2_genes)
        raise
    if weight >= args.threshold:
        if c1[0] == c2[0]:
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


def graph_unique_change_when_adding(species, cluster2genes, reverse = True):
    print('Graph of the change of unique clusters fraction with each new added genome.')
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
        if first: # no need to calculate anything, ratio is 1.0
            first = False
            ratio = 1.0
        else:
            unique = get_unique_clusters(s, allowed_species)
            ratio = round(float(unique) / numclust, 2)
        allowed_species.append(s)
        bar = '#' * int(round(height*ratio))
        bar = bar.ljust(height)
        print('%s\t%s\t%s' % (bar, ratio, s))

def cumulative_growth(species, cluster2genes, reverse = True):
    '''Shows expected and observed growth of the number of unique clusters
    with each new genome'''
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
    for (numclust, s) in numclust_species:
        total_expected += len(cluster2genes[s])
        total_observed += get_unique_clusters(s)
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


def parse_gene_cluster_relations(inputs, args, cluster2genes, gene2clusters, numbers2products, coords2numbers, clustersizes, all_clusters):
    '''
    All of: cluster2genes, gene2clusters, numbers2products, coords2numbers,
            clustersizes, all_clusters
    will be updated upon execution of this function.
    '''
    # Dict of per-species GenBank records, e.g. genbank[id] = [rec1, rec2, ...]
    genbank = {}
    # Dict of per-species interval trees of clusters, one per record.
    clustertrees = {}

    species = inputs.keys()
    species_count = len(species)
    logging.info('Parsing genbank files in parallel.')
    task_queue = Queue() # record.id + path to antismashed genbank, or 'STOP'.
    done_queue = Queue()

    # Parallel worker.
    def worker(tasks, done):
        try:
            (k, gb) = tasks.get() # By default, there is no timeout.
        except: # Queue.Empty
            logging.warning("worker(tasks, done) encountered an exception trying tasks.get()")
            raise
        # Exit if 'STOP' element is found.
        if gb == 'STOP':
            logging.debug("STOP found, exiting.")
            return
        # FIXME: add all of parsing in here.
        done.put((k, list(SeqIO.parse(gb, "genbank", generic_dna))))

    logging.info("Starting %s GenBank parse workers.", cpu_count())
    for _ in range(min(cpu_count(), species_count)):
        # TODO: do I need to join() these processes later?
        Process(target=worker, args=(task_queue, done_queue)).start()
    del _

    logging.debug("Putting ID-path tuples into task_queue")
    for k, v in inputs.iteritems():
        # submit record.id and /path/to/genbank
        task_queue.put((k, v['as2file']))
    del k, v

    logging.debug("Adding %s STOP messages to task_queue.", cpu_count())
    for _ in range(cpu_count()):
        task_queue.put(0, 'STOP')
    del _

    for _ in range(len(species)):
        k, rec_list = done_queue.get()
        genbank[k] = rec_list
    del _, rec_list

    print('\ttotal records parsed: %s' % len(genbank))
    task_queue.close()
    done_queue.close()


    # FIXME: merge into the genbank parsing worker,
    # to remove genbanks from memory when done.
    print('Parsing clusters and assigning genes to them:')
    logging.info('\tgetting extension sizes for diff. cluster types from antismash2 config')
    rulesdict = hmm_detection.create_rules_dict()
    if args.trim:
        logging.info('\tNote: cluster coordinates are shown after trimming, except for skipped putative clusters.')
    for s in species:
        logging.debug('\tprocessing %s (%s)', s, inputs[s]['species'])
        clustertrees[s] = {}
        cluster2genes[s] = {}
        numbers2products[s] = {}
        coords2numbers[s] = {}
        gene2clusters[s] = {}
        # Populate clusters tree and dict with (start, end) as keys.
        for r in genbank[s]: # genbank[s] is a list of records
            record_number = genbank[s].index(r)
            clustertrees[s][record_number] = IntervalTree()
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
                        extension = hmm_detection.get_extension_size_by_cluster_type(f.qualifiers['product'][0], rulesdict)
                        # If cluster is shorter than double extension, do not trim.
                        if 2 * extension >= (end - start):
                            logging.info('\tnot trimming #%s (%s) - shorter than double extension 2 * %s',
                                         parse_cluster_number(f.qualifiers['note']),
                                         f.qualifiers['product'][0], extension)
                        else:
                            # If cluster is at the genome edge - skip extension trimming.
                            if start == 0:
                                logging.info('\tnot trimming left-extension for #%s (%s) - at the genome start',
                                             parse_cluster_number(f.qualifiers['note']), f.qualifiers['product'][0])
                            else:
                                start += extension
                            if end == len(r):
                                logging.info('\tnot trimming right-extension for #%s (%s) - at the genome end',
                                             parse_cluster_number(f.qualifiers['note']), f.qualifiers['product'][0])
                            else:
                                end -= extension
                            try:
                                assert start < end
                            except:
                                logging.exception('trimming extension failed for: %s (%s), extension %s, ori start %s, new start %s, ori end %s, new end %s',
                                                  f.qualifiers['product'][0],
                                                  parse_cluster_number(f.qualifiers['note']),
                                                  extension, f.location.start.position,
                                                  start, f.location.end.position,
                                                  end)
                                raise
                    clustertrees[s][record_number].add_interval(Interval(start, end))
                    cluster2genes[s][cluster_number] = []
                    all_clusters.append((s, cluster_number))
                    numbers2products[s][cluster_number] = f.qualifiers['product'][0]
                    # Start/end coords theoretically can collide between
                    # multiple loci of the same genome: check for this.
                    assert (start, end) not in coords2numbers[s]
                    coords2numbers[s][(start, end)] = cluster_number
                    clustersizes[(s, cluster_number)] = end - start
                    logging.info('''\t('%s', %s): %s [%s, %s], %s bp long''', s, cluster_number,
                                 f.qualifiers['product'][0], start, end, end-start)
        logging.info('\tAssign genes to biosynthetic clusters')
        num_genes = 0
        for r in genbank[s]: # genbank[s] is a list of records
            record_number = genbank[s].index(r)
            for f in r.features:
                if f.type == 'CDS':
                    # 'cl' is a list of Intervals, each has 'start' and 'end' attributes.
                    cl = clustertrees[s][record_number].find(int(f.location.start.position), int(f.location.end.position))
                    if len(cl) > 0:
                        # Determine which qualifier type to use for gene name.
                        if 'locus_tag' in f.qualifiers:
                            qualifier = 'locus_tag'
                        elif 'gene' in f.qualifiers:
                            qualifier = 'gene'
                        gene_name = f.qualifiers[qualifier][0] + '.' + s
                        # Check for gene_name overlaps between loci.
                        assert gene_name not in gene2clusters[s]
                        # One gene may belong to more than one cluster.
                        logging.debug('\t\tgene at (%s, %s) overlaps with %s cluster(s), 1st is at (%s, %s)',
                                      f.location.start.position, f.location.end.position,
                                      len(cl), cl[0].start, cl[0].end)
                        num_genes += 1
                        gene2clusters[s][gene_name] = []
                        for cluster in cl:
                            cluster2genes[s][coords2numbers[s][(cluster.start, cluster.end)]].append(gene_name)
                            gene2clusters[s][gene_name].append((s, coords2numbers[s][(cluster.start, cluster.end)]))
                        del gene_name
        print('\t%s: %s clusters populated with %s genes' % (s, len(cluster2genes[s]), num_genes))
        del num_genes
    # Extra verification.
    for s in species:
        for num, cl in cluster2genes[s].iteritems():
            try:
                assert len(cl) > 0
            except:
                print('cluster %s of %s has 0 genes' % (num, s))
                raise
    print('\tadded %s clusters from %s species' % (len(all_clusters), species_count))
    # Sort by species.lower() to ensure stable order of cluster pairs.
    # FIXME: to implement proper sorting, may need to use faafiles instead of r.id as s/species.
    all_clusters.sort(key = lambda s: s[0].lower())

    logging.debug('Freeing memory.')
    del genbank, clustertrees

#
# /Set of functions used by process()
#

def process(inputs, paranoid, args):
    '''
    Analyze quickparanoid results. 'args' contains:
    "paranoid" is the path to quickparanoid output file.
    "paths" is a list of paths to genbank files we want to compare.
    "threshold" is a biocluster-biocluster link weight threshold.
    "prefix" is prepended to all output files.
    "trim", if True, causes antismash2 clusters to lose non-core extensions at
    both ends of the cluster (these are hard-coded in antismash2). Note that for composite-type
    clusters only the shorter of the extensions will be trimmed (e.g. bacteriocin extension
    for the backteriocin-t1pks cluster).
    "skipp" means "skip putative clusters", if set.
    "no_tree_problems", if True, will only use orthology clusters which do not have any problems in the tree_conflict column.
    "no_name_problems": as above, but only for the 'diff. names' problem in the tree_conflict column.
    "use_sizes" will weigh each clusters contribution to final weight according to clusters length proportion of total length, in bp.
    "scale" will scale down weight by the ratio of physical cluster lengths: min(size1, size2)/max(size1, size2).
    '''
    # FIXME: re-think this entire function.
    # TODO: modularize?...
    # TODO: replace species with either inputs.keys() or faafiles.
    # TODO: replace paths with args.paths, or inputs?
    # TODO: replace 'trim' support with no-extensions support (?).
    # TODO: deprecate skipp? (can be filtered out at the analysis step)
    # TODO: add memory cleanup everywhere :)
    # TODO: switch from record.name to record.id.

    mp = MP.MultiParanoid(paranoid, args.no_tree_problems, args.no_name_problems)

    # Declare important variables.
    # List of recognized species IDs.
    species = inputs.keys()
    # List of all clusters from all species.
    all_clusters = []
    # Map of per-species cluster numbers to their products
    numbers2products = {}
    # Per-species mapping of cluster coordinates tuple to their numbers.
    coords2numbers = {}
    # 2-level nested dict of cluster pairs link weights,
    # e.g. cluster_weights['A'] = {'B': 0.95, ...}
    cluster_weights = {}
    # Same as above, but without duplicate links to other species.
    weights_clean = {}
    # Same as cluster_weights, but only for intra-species links.
    weights_intra = {}
    # Mapping of record.id from GenBank files (LOCUS) to species.
    # DEPRECATED: use inputs[primary_id]['species']
    #locus2species = {}
    # Dict of per-species dicts of cluster-to-gene relations.
    # Each cluster dict uses a key = (species, number).
    cluster2genes = {}
    # Inverse of the above: per-species mapping of each gene to the cluster(s) it belongs to.
    gene2clusters = {}
    # Two lists of 2-tuples with weight > threshold links between clusters.
    intra_one = []
    inter_one = []
    # Sizes of clusters, in basepairs. clustersizes[(species, number)] = 65450
    clustersizes = {}

    # Read genbank files, assign genes to clusters, clusters to genes, and
    # fill other useful maps and lists.
    parse_gene_cluster_relations(inputs, args, cluster2genes, gene2clusters,
                                 numbers2products, coords2numbers, clustersizes,
                                 all_clusters)

    print('Constructing gene-based cluster links.')
    pairs = combinations(all_clusters, 2)
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
    # have high scores of up 1.0. Removing them may adversely affect the
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

    print('All pairs of clusters with link weight over %s (inter-species).', args.threshold)
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

    graph_unique_change_when_adding()
#    graph_unique_change_when_adding(False)

    cumulative_growth()
#    cumulative_growth(False)

    # Finally, show a list of all species, stating the number of unique
    # clusters they have.
    print('Graph of the percentage of unique clusters in each of the genomes.')
    height = 60
    # Draw a reference 100% line.
    print('%s\t%s\t%s\t%s\t%s' % ('Bar'.ljust(height), 'Unique', 'Total', 'Ratio', 'Genome'))
    print('%s\t%s\t%s\t%s\t%s' % ('#' * height, 1.0, 1.0, 1.0, 'Reference'))
    for s in species:
        unique = get_unique_clusters(s)
        total = len(cluster2genes[s])
        ratio = round(unique / float(total), 2)
        bar = '#' * int(round(height*ratio))
        bar = bar.ljust(height)
        print('%s\t%s\t%s\t%s\t%s' % (bar, unique, total, ratio, s))


def preprocess_input_files(inputs, args):
    '''
    inputs: dictionary to populate
    args.paths: list of genbank files to process
    '''
    # FIXME: split into 2 parts:
    # - linear processing of genbank IDs
    # - parallel translation and antismashing
    for infile in args.paths:
        # TODO: convert this to a worker to run in parallel.
        contigs = 0 # number of fragments of the genome (can be contigs, plasmids, chromosomes, etc)
        genome_size = 0 # total size of all contigs
        organism = {} # maps accessions to 'organism' field
        primary_accession = '' # for the accession of the longest record
        primary_id = '' # as above
        primary_length = 0 # current primary accession sequence length
        accessions = [] # other accessions, e.g. plasmids/contigs/etc
        ids = [] # as above
        # Extract key information.
        for r in SeqIO.parse(infile, 'genbank', generic_dna):
            # r.name is an accession (e.g. AF080235),
            # r.id is a versioned accession (e.g. AF080235.1)
            # We use r.id as more specific.
            logging.debug('record name and ID are %s and %s', r.name, r.id)
            if not r.id or r.id == '':
                logging.error('Genome %s has no usable ID!' % input)
                raise Exception('GenomeHasNoIDError')
            contigs += 1
            accessions.append(r.name)
            ids.append(r.id)
            organism[r.name] = r.annotations['organism']
            if primary_id == '':
                primary_accession = r.name
                primary_id = r.id
                primary_length = len(r.seq)
            elif len(r.seq) > primary_length:
                primary_accession = r.name
                primary_id = r.id
                primary_length = len(r.seq)
            genome_size += len(r.seq)
        if '-' in primary_id:
            logging.error('quickparanoid allows no dashes in filenames; genome "%s" has a dash!',
                          primary_id)
            raise Exception('DashesInIDsAreForbiddenError')
        del primary_length, r
        # Remove primary accession from the list of all accessions.
        accessions.remove(primary_accession)
        ids.remove(primary_id)

        # Check for duplicate ID.
        if primary_id in inputs:
            # Check for an exact duplicate.
            if (inputs[primary_id]['contigs'] == contigs and
                inputs[primary_id]['genome_size'] == genome_size and
                inputs[primary_id]['species'] == organism[primary_accession]):
                logging.error('Found identical genomes of "%s" with ID %s, %s bp long (in %s fragments).',
                              organism[primary_accession], primary_id, genome_size, contigs)
                raise Exception('IdenticalGenomesError')
            else:
                logging.debug('Append _1, _2 etc to the ID until it becomes unique.')
                for i in range(1, 10):
                    new_id = primary_id + '_' + str(i)
                    if new_id in inputs:
                        continue
                    else:
                        primary_id = new_id
                        del new_id

        # TODO: comment out or remove unneeded fields when the program is finished.
        inputs[primary_id] = {}
        inputs[primary_id]['contigs'] = contigs
        inputs[primary_id]['genome_size'] = genome_size
        inputs[primary_id]['infile'] = infile
        inputs[primary_id]['species'] = organism[primary_accession]
        inputs[primary_id]['accessions'] = accessions
        inputs[primary_id]['ids'] = ids

        # Write FASTA to ID.fna, do nothing if file exists.
        fnafile = join(args.project, primary_id + '.fna')
        inputs[primary_id]['fnafile'] = fnafile
        gb2fasta(infile, fnafile, primary_id)

        # Re-use antismash2 annotation if it exists.
        as2file = join(args.project, primary_id + '.gbk')
        inputs[primary_id]['as2file'] = as2file
        antismash2_reused = False
        if args.force and exists(as2file):
            logging.warning('Reusing existing antismash2 annotation; --no-extensions option will NOT be honored!')
            logging.warning('Do not use --force, or delete antismash2 *.gbk files to re-run antismash2 annotation.')
            # Re-using any further files is only possible if we do re-use antismash files.
            antismash2_reused = True
        else:
            output_folder = join(args.project, primary_id)
            as2_options = ['run_antismash', '--outputfolder', output_folder]
            # Removed time-consuming/unnecessary options: smcogs, clusterblast,
            # subclusterblast.
            # Removed unnecessary options: --full-blast, --full-hmmer.
            as2_options.extend(['--cpus', '1', '--verbose', '--all-orfs'])
            as2_options.extend(['--input-type', 'nucl'])
            if args.no_extensions:
                logging.warning("using --no-extensions option; cluster.py DOES NOT check if this option is supported!")
                as2_options.append('--no-extensions')
            as2_options.append(fnafile)
            logging.info('Running antismash2: %s', ' '.join(as2_options))
            # TODO: show output from child processes?...
            out, err, retcode = utils.execute(as2_options)
            if retcode != 0:
                logging.debug('antismash2 returned %d: %r while scanning %r',
                              retcode, err, fnafile)
            # antismash's algorithm for naming the output file:
            # basename = seq_records[0].id
            # output_name = path.join(options.outputfoldername, "%s.final.gbk" % basename)
            logging.debug('as2file (target): %s', as2file)
            antismash2_file = join(output_folder, primary_id + '.final.gbk')
            logging.debug('antismash2 file (source): %s', antismash2_file)
            rename(antismash2_file, as2file)
            rmtree(output_folder)
            inputs[primary_id]['as2file'] = as2file
            del output_folder, out

        faafile = join(args.project, primary_id + '.faa')
        inputs[primary_id]['faafile'] = faafile
        if antismash2_reused and args.force and exists(faafile):
            logging.warning('Reusing existing translations file!')
        else:
            extract_translation_from_genbank(as2file, faafile, False)
        del as2file, fnafile, faafile, infile
        del contigs, genome_size, organism, primary_accession, accessions
        del ids, primary_id
    # When all the workers have finished: make sure they exit, using a keyword.

def prepare_inparanoid(inputs, args):
    '''
    Set up directory and .faa-files symlink for running InParanoid and
    quickparanoid. Prepare list of .faa-files w/o paths, and return it.
    '''
    # Collect all faafile names into a list without paths.
    faafiles = []
    for _ in inputs.itervalues():
        faafiles.append(basename(_['faafile']))
    del _
    faafiles.sort(key = lambda s: s.lower())

    # Put everything inparanoid-related into a subdir.
    inparanoidir = realpath(join(args.project, 'inparanoid'))
    if not (args.force and exists(inparanoidir)):
        mkdir(inparanoidir)

    # Symlink all faa files there.
    for _ in faafiles:
        # TODO: can parallelize.
        linkname = join(inparanoidir, _)
        if not (args.force and exists(linkname)):
            symlink(join('..', _), linkname)
    del linkname, _

    return inparanoidir, faafiles

def run_inparanoid(inparanoidir, faafiles):
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

    for _ in faafiles:
        blast_single = ['inparanoid.pl', '--blast-only', _]
        logging.info('Running blast on a single genome: %s',
                     ' '.join(blast_single))
        out, err, retcode = utils.execute(blast_single)
        if retcode != 0:
            logging.debug('inparanoid returned %d: %r while blasting %r',
                          retcode, err, _)
    del _, blast_single, out, err, retcode

    for pair in permutations(faafiles, 2):
        blast_pair = ['inparanoid.pl', '--blast-only', pair[0], pair[1]]
        logging.info('Running blast on genome pair: %s',
                     ' '.join(blast_pair))
        out, err, retcode = utils.execute(blast_pair)
        if retcode != 0:
            logging.debug('inparanoid returned %d: %r while blasting %r and %r',
                          retcode, err, pair[0], pair[1])
    del blast_pair, pair, out, err, retcode

    # FIXME: make run in parallel.
    for pair in combinations(faafiles, 2):
        inparanoid_pair = ['inparanoid.pl', pair[0], pair[1]]
        logging.info('Running inparanoid analysis: %s', ' '.join(inparanoid_pair))
        out, err, retcode = utils.execute(inparanoid_pair)
        if retcode != 0:
            logging.debug('inparanoid returned %d: %r while analyzing %r and %r',
                          retcode, err, pair[0], pair[1])
    del inparanoid_pair, pair, out, err, retcode

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

    # Alternative path finding method: dirname(sys.argv[0])
    quickparanoid = join(dirname(realpath(__file__)), 'quickparanoid')
    curr_path = getcwd()
    logging.debug("Remembering current directory %s, changing to %s.",
                  curr_path, quickparanoid)
    chdir(quickparanoid)

    # quickparanoid requires a trailing slash.
    qp = ['./qp', inparanoidir + os.sep, configfile, project]
    logging.info('Running quickparanoid analysis: %s', ' '.join(qp))
    out, err, retcode = utils.execute(qp)
    if retcode != 0:
        logging.debug('quickparanoid returned %d: %r while analyzing %r in %r',
                      retcode, err, configfile, project)
    del out, err, retcode

    # Move generated 'project' and 'projects' executables to the {project} directory.
    rename(join(quickparanoid, project), join(curr_path, project, project))
    rename(join(quickparanoid, project + 's'),
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
    # FIXME: --trim is broken (on by default, impossible to turn off).
    parser.add_argument("--trim", dest="trim", action="store_false", default=True, help="trim away antismash2 cluster extensions [default: %(default)s]")
    parser.add_argument("--skip-putative", dest="skipp", action="store_true", default=False, help="exclude putative clusters from the analysis [default: %(default)s]")
    parser.add_argument("--strict", dest="strict", action="store_true", default=False, help="weight between clusters with 5 and 10 genes will never exceed 0.5 [default: %(default)s]")
    parser.add_argument("--scale", dest="scale", action="store_true", default=False, help="scale link weight down by a factor of min(size1, size2)/max(size1, size2) [default: %(default)s]")
    parser.add_argument("--use-sizes", dest="use_sizes", action="store_true", default=False, help="each cluster's contribution to link weight is scaled by relative cluster sizes; can be combined with --strict [default: %(default)s]")
    parser.add_argument("--no-name-problems", dest="no_name_problems", action="store_true", default=False, help="only use ortho-clusters which do not have diff.names tree_conflict problems [default: %(default)s]")
    parser.add_argument("--no-tree-problems", dest="no_tree_problems", action="store_true", default=False, help="only use ortho-clusters which do not have [diff.names/diff.numbers] tree_conflict problems [default: %(default)s]")
    parser.add_argument("--prefix", default='out', help="output CSV files prefix [default: %(default)s]")
    parser.add_argument("--project", default='cluster_project', help="put all the project files into this directory [default: %(default)s]")
    parser.add_argument('--force', action = 'store_true', default = False, help='insist on re-using existing project directory (this will re-use existing intermediate files) [default: %(default)s]')
    parser.add_argument('--no-extensions', action = 'store_true', default = False, help='pass --no-extensions option to the modified antismash2 (see README for detais) [default: %(default)s]')
    parser.add_argument('--threshold', action = 'store', type=float, default = 0.0, help='cluster links with weight below this one will be discarded [default: %(default)s]')
    parser.add_argument('--height', action = 'store', type=int, default = 50, help='bar heights for text graphs [default: %(default)s]')
    # FIXME: species will be read from input genbank files
#    parser.add_argument(dest="species", help="path to the plain-text species list file", metavar="species")
    # FIXME: multi/quick-paranoid results will be generated by the program
#    parser.add_argument(dest="paranoid", help="path to the multiparanoid/quickparanoid sqltable file", metavar="sqltable")
    parser.add_argument(dest="paths", help="paths to the GenBank files with genomes to analyze", metavar="path", nargs='+')
    args = parser.parse_args()

    # Remove any possible duplicate input filenames.
    args.paths = list(set(args.paths))

    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    else:
        level = logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s::%(levelname)s::%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

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

    # "Catalog" all input files.
    inputs = {} # Map genome ID to other properties.

    if len(args.paths) == 1:
        logging.warning("Single input file specified, program will exit after preprocessing.")
    preprocess_input_files(inputs, args)
    if len(args.paths) == 1: # single input - exit
        sys.exit(3)
    inparanoidir, faafiles = prepare_inparanoid(inputs, args)
    run_inparanoid(inparanoidir, faafiles)

    result_path = run_quickparanoid(inparanoidir, faafiles, args.project)
    del inparanoidir, faafiles

    # Process result_path.
    process(inputs, result_path, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
