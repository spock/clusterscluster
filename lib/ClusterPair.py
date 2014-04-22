from __future__ import print_function
import logging
import itertools
from scipy import stats
from tempfile import NamedTemporaryFile
from collections import namedtuple
from utils import usearch, SymKeyDict
from lib.Genome import GeneOrder, FeatureTuple


cluster = namedtuple('Cluster', ['genome', 'number'])


class ClusterPair(object):
    '''
    A pair of potentially similar biosynthetic clusters.
    '''


    def __init__(self, gc1, gc2):
        '''
        gc1, gc2: genome-clusternumber namedtuple pairs for clusters 1 and 2.
        '''
        self.gc1 = gc1
        self.gc2 = gc2
        self.g1 = gc1.genome
        self.c1 = gc1.number
        self.g2 = gc2.genome
        self.c2 = gc2.number
        # orthology links
        self.link1 = 0.0
        self.link2 = 0.0
        # orthology weight
        self.weight = 0.0
        # is this cluster pair intra-species?
        if self.g1 == self.g2:
            self.intra = True
        else:
            self.intra = False
        # Counters of genes in c1 and c2.
        self.c1_genes = 0
        self.c2_genes = 0
        # Gene-level protein identities, identities[(g1, g2)] = float,
        # symmetric (setting (g1, g2) makes (g2, g1) also accessible);
        # g1 and g2 are locus_tag:genome_id strings.
        self.protein_identities = SymKeyDict()
        # Single tuple for average gene-level protein identity,
        # (number_of_gene_pairs, average_identity)
        self.avg_identity = None
        # Dict of genome1 genes mapped to genome2 genes (by locus_tag:genomeid),
        # only for genes which have above-cutoff identity.
        self.gene1_to_gene2 = {}


    def num_c1_genes(self, genomes):
        '''Return the number of genes in cluster c1.'''
        if self.c1_genes == 0:
            self.c1_genes = len(genomes[self.g1].cluster2genes[self.c1])
        return self.c1_genes


    def num_c2_genes(self, genomes):
        '''Return the number of genes in cluster c2.'''
        if self.c2_genes == 0:
            self.c2_genes = len(genomes[self.g2].cluster2genes[self.c2])
        return self.c2_genes


    def get_gene_links_to_bioclusters(self, gene, mp, genomes):
        '''
        For the given `gene`,
        - find to which orthology group/cluster from `mp` it belongs,
        - check if other genes from that group are a part of some biosynthetic clusters,
        - return the list of those biosynthetic clusters, where
        - each is a Cluster(genome, number) named tuple.
        '''
        if gene not in mp.gene2ortho:
#            logging.debug('Gene %s is NOT in mp.gene2ortho.', gene)
            return [] # `gene` does not belong to any group of orthologs
#        logging.debug('Gene %s is in mp.gene2ortho.', gene)
        bioclusters = []
        for orthoclust in mp.gene2ortho[gene]:
#            logging.debug('gene %s belongs to ortho-cluster %s', gene, orthoclust)
            for xeno_gene in mp.ortho2genes[orthoclust]:
                # Extract LOCUS from gene name, find species from it.
                xeno_species = xeno_gene.rsplit(':')[-1]
                # Check if xeno_gene belongs to any biosynthetic clusters.
                xeno_gene2clusters = genomes[xeno_species].gene2clusters
                if xeno_gene in xeno_gene2clusters:
                    for clnum in xeno_gene2clusters[xeno_gene]:
                        xeno_cluster = cluster(xeno_species, clnum)
                        bioclusters.append(xeno_cluster)
#                        logging.debug('\txeno_gene %s belongs to %s', xeno_gene, xeno_cluster)
        return bioclusters


    def calculate_links(self, which, mp, genomes):
        '''
        Count the number of unique orthologoues gene pairs ("links") between
        clusters in this pair.
        mp: QuickParanoid object.
        genomes: dict of all known genomes.
        'which': which genome/cluster is the basis for calculation;
                 result will change depending on this!
        '''
        links = 0
        if which == 1:
            logging.debug("Looking at %s genes in cluster %s (2nd is %s)...",
                          len(genomes[self.g1].cluster2genes[self.c1]),
                          self.c1, self.gc2)
            for gene in genomes[self.g1].cluster2genes[self.c1]:
                # 'gene' is a 'locus_tag:genome_id' string.
                if self.gc2 in self.get_gene_links_to_bioclusters(gene, mp, genomes):
                    links += 1
        if which == 2:
            logging.debug("Looking at %s genes in cluster %s (1st is %s)...",
                          len(genomes[self.g2].cluster2genes[self.c2]),
                          self.c2, self.gc1)
            for gene in genomes[self.g2].cluster2genes[self.c2]:
                if self.gc1 in self.get_gene_links_to_bioclusters(gene, mp, genomes):
                    links += 1
        return links


    def assign_orthologous_link(self, mp, genomes, args):
        '''
        Assign orthologous link to this cluster pair.
        '''
        link1 = float(min(self.calculate_links(1, mp, genomes),
                           self.num_c1_genes(genomes)))
        link2 = float(min(self.calculate_links(2, mp, genomes),
                           self.num_c2_genes(genomes)))
        logging.debug('\tlink1= %s and link2= %s for %s and %s, %s/%s genes',
                      link1, link2, self.gc1, self.gc2, self.num_c1_genes(genomes),
                      self.num_c2_genes(genomes))
        if args.strict:
            # No link can be larger than the number of genes in the 2nd cluster.
            if link1 > self.num_c2_genes(genomes):
                link1 = float(self.num_c2_genes(genomes))
                logging.debug('strict mode, new link1 is %s', link1)
            if link2 > self.num_c1_genes(genomes):
                link2 = float(self.num_c1_genes(genomes))
                logging.debug('strict mode, new link2 is %s', link2)
            try:
                assert link1 <= min(self.num_c1_genes(genomes),
                                    self.num_c2_genes(genomes)) and link2 <= min(self.num_c1_genes(genomes),
                                                                                 self.num_c2_genes(genomes))
            except:
                logging.exception('c1: %s ; c2: %s', self.gc1, self.gc2)
                logging.exception('c1_genes: %s (c1: %s)', self.num_c1_genes(genomes),
                                  genomes[self.g1].cluster2genes[self.c1])
                logging.exception('c2_genes: %s (c2: %s)', self.num_c2_genes(genomes),
                                  genomes[self.g2].cluster2genes[self.c2])
                logging.exception('link1: %s ; link2: %s ; c1_genes: %s ; c2_genes: %s',
                                  link1, link2, self.num_c1_genes(genomes),
                                  self.num_c2_genes(genomes))
                raise
        self.link1 = link1
        self.link2 = link2


    def CDS_identities(self, genomes):
        '''
        Calculate per-gene-pair protein identity, and also average
        protein identity for cluster pair. Calculate average protein identity.
        Populate lists self.genes_with_pairs[genome] with geneids which have pairs.
        Save all values to self.
        '''
        GP = namedtuple('GenePair', ['identity', 'g1', 'g2'])
        # Get genelists for both clusters.
        gl1 = genomes[self.g1].cluster2genes[self.c1]
        gl2 = genomes[self.g2].cluster2genes[self.c2]
        # Iterate all pairs of genes of this cluster pair.
        gene_pairs = []
        # Make /tmp file and open it for writing.
        with NamedTemporaryFile(mode='w', dir='/tmp') as seqfile:
            # FIXME: avoid comparing cl. 2 to 9, and then later 9 to 2, when genome is the same
            for gene1, gene2 in itertools.product(gl1, gl2):
                #logging.debug('gene1, gene2: %s, %s', gene1, gene2)
                seq1 = genomes[self.g1].get_protein(gene1.split(':')[0])
                seqfile.write(">%s\n%s\n" % (gene1, seq1))
                seq2 = genomes[self.g2].get_protein(gene2.split(':')[0])
                seqfile.write(">%s\n%s\n" % (gene2, seq2))
                assert len(seq1) > 0 and len(seq2) > 0
            # seqfile is open for writing when usearch runs, so better at least flush it
            seqfile.flush()
            # Run usearch, once. TODO: do FULL later for chosen genepairs?
            results = usearch(seqfile.name, full = False)
            if results == '':
                # empty result: nothing above the cut-off, no similar gene pairs
                self.avg_identity = (0, 0.0)
                seqfile.close()
                return
            if results == None:
                # error in usearch
                logging.exception('usearch failed, see output above')
                raise Exception('UsearchFailedError')
            # Parse results into a list of tuples, each tuple is 1 gene pair.
            results_list = results.strip('\n').split('\n')
            for row in results_list:
                try:
                    # query, target should be from g1, g2, respectively
                    query, target, identity = row.split('\t')
                    gene_pairs.append(GP(identity, query, target))
                    self.gene1_to_gene2[query] = target
                except ValueError:
                    print('row:')
                    print(row)
                    logging.exception('Failed to parse usearch output.')
                    raise Exception('UsearchOutputParseError')
            del results, results_list, identity, query, target
        # Sort gene pairs by identity, descending order.
        gene_pairs.sort(reverse = True)
#        print(gene_pairs)
        # Two lists to check that we have not yet seen genes from c1 and c2.
        seen_1 = []
        seen_2 = []
        # Counter of gene pairs.
        num_pairs = 0
        # Cumulative sum of identities, for average.
        sum_identities = 0.0
        for pair in gene_pairs:
            if pair.g1 in seen_1 or pair.g2 in seen_2:
                # at least one of the genes already has a pair
                logging.debug('either %s or %s already has a pair', pair.g1, pair.g2)
                continue
            # else: save a new gene identity pair!
            self.protein_identities[(pair.g1, pair.g2)] = pair.identity
            seen_1.append(pair.g1)
            seen_2.append(pair.g2)
            num_pairs += 1
            sum_identities += float(pair.identity)
        self.avg_identity = (num_pairs, sum_identities / num_pairs)
        del gene_pairs, seen_1, seen_2, sum_identities, num_pairs


    def gene_order(self, genomes):
        '''
        Are similar genes in the same order or not?
        '''
        # FIXME: some gene indices are present more than once!!!! Examples:
        # {Cluster(genome='GMON_MP13002', number=6): [15, 14, 13, 10, 13], Cluster(genome='KRI_MP130_17', number=9): [6, 7, 8, 9, 10]}
        # {Cluster(genome='SBA_MP131_18', number=35): [2, 21, 23, 24, 25], Cluster(genome='SBA_MP131_18', number=5): [18, 23, 20, 20, 20]}
        # {Cluster(genome='SBA_MP131_18', number=10): [1, 4, 6, 9], Cluster(genome='SBA_MP131_18', number=33): [22, 22, 22, 22]}
        # For clusters 1 and 2, for genes which have identity pairs,
        # build a list of per-genome gene indexes (i.e. 0, 1, 2...)
        # (from lower to higher cluster coordinate).
        inds = {self.gc1: [], self.gc2: []}
        # 'gene' is a GeneOrder namedtuple, (start, strand, locus_tag:genome_id).
        for gene in genomes[self.g1].orderstrands[self.c1]: # iterate all cluster 1 genes
            if gene.geneid in self.gene1_to_gene2: # check if this gene has an identity pair
                # Offset indices by 1 (they start from 0, might be problematic for stats).
#                print(genomes[self.g1].orderstrands[self.c1])
                inds[self.gc1].append(genomes[self.g1].orderstrands[self.c1].index(gene) + 1)
                # Get the corresponding gene position from cluster 2.
                gene2_name = self.gene1_to_gene2[gene.geneid]
                gene2_feature = genomes[self.g2].get_feature_by_locustag(gene2_name.split(':')[0]).feature
                gene2 = GeneOrder(int(gene2_feature.location.start.position), gene2_feature.strand, gene2_name)
                inds[self.gc2].append(genomes[self.g2].orderstrands[self.c2].index(gene2) + 1)
#                print(genomes[self.g2].orderstrands[self.c2])
        print(inds) # FIXME: debug only
        # spearman: monotonicity (same gene order, disregarding distances), less sensitive to outliers
        # pearson: linearity (same gene order AND similar distances); not really applicable: "a measure of the linear relationship between two continuous random variables"
        # kendall: unlike spearman and pearson, here each point contributes equally
        # magnitude: kendall < spearman
        print('S:', stats.spearmanr(inds[self.gc1], inds[self.gc2])[0]) # FIXME: debug only
        print('P:', stats.pearsonr(inds[self.gc1], inds[self.gc2])[0]) # FIXME: debug only
        print('K:', stats.kendalltau(inds[self.gc1], inds[self.gc2])[0]) # FIXME: debug only
        # BioPython also has these stat routines:
#        import numpy as np
#        from Bio import Cluster
#        array1 = np.array(inds[self.g1])
#        array2 = np.array(inds[self.g2])
#        self.spearman = 1 - Cluster.distancematrix((array1, array2), dist="s")[1][0]
#        self.pearson = 1 - Cluster.distancematrix((array1, array2), dist="c")[1][0]
#        self.kendall = 1 - Cluster.distancematrix((array1, array2), dist="k")[1][0]
        # for strands: split every list in 2 (by strand), run stats on list pairs, pick the highest score


    def domains(self, genomes):
        '''
        Calculate collinearity of predicted domains and their substrates
        in similar genes.
        '''
        pass

    def nucleotide_similarity(self, genomes):
        pass
