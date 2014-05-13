from __future__ import print_function
import logging
import itertools
from scipy import stats
from tempfile import NamedTemporaryFile
from collections import namedtuple
from lib.Genome import GeneOrder
from lib.MultiParanoid import gene2species


Cluster = namedtuple('Cluster', ['genome', 'number'])
GP = namedtuple('GenePair', ['identity', 'g1', 'g2'])
Avg_Identity = namedtuple('Avg_Identity', ['num_gene_pairs', 'avg_identity'])


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
        self.link1 = 0
        self.link2 = 0
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
        # g1 and g2 are locus_tag:genome_id strings. It is assigned later as
        # SymKeyDict, so here it is declared as a simple dict.
        self.protein_identities = {}
        # Single tuple for average gene-level protein identity,
        # (number_of_gene_pairs, average_identity)
        self.avg_identity = Avg_Identity(0, 0.0)
        # Dict of genome1 genes mapped to genome2 genes (by locus_tag:genomeid),
        # only for genes which have above-cutoff identity.
        self.gene1_to_gene2 = {}
        # Correlation coefficients for gene order preservation.
        self.kendall = 0.0
        self.pearson = 0.0
        self.spearman = 0.0


    # TODO: these custom hash/eq are not sufficient for caching: number of genes in the clusters can be different.
    def __hash__(self):
        return hash((self.g1, self.c1, self.g2, self.c2))


    def __eq__(self, other):
        return (self.g1, self.c1, self.g2, self.c2) == (other.g1, other.c1, other.g2, other.c2)


    def num_c1_genes(self, genome1):
        '''Return the number of genes in cluster c1. genome1 = genomes[g1].'''
        if self.c1_genes == 0:
            self.c1_genes = len(genome1.cluster2genes[self.c1])
        return self.c1_genes


    def num_c2_genes(self, genome2):
        '''Return the number of genes in cluster c2. genome2 = genomes[g2].'''
        if self.c2_genes == 0:
            self.c2_genes = len(genome2.cluster2genes[self.c2])
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
        # FIXME: collapse calls to gene2ortho and ortho2genes into gene2genes:
        # gene2genes[gene.genome][gene.id] = [other.gene1, other.gene2, ...]
        # TODO: maybe store gene names not as genome:gene, but as (genome, gene)?
        for orthoclust in mp.gene2ortho[gene]:
#            logging.debug('gene %s belongs to ortho-cluster %s', gene, orthoclust)
            for xeno_gene in mp.ortho2genes[orthoclust]:
                # Extract LOCUS from gene name, find species from it.
                xeno_species = xeno_gene.rsplit(':')[-1]
                # Check if xeno_gene belongs to any biosynthetic clusters.
                xeno_gene2clusters = genomes[xeno_species].gene2clusters
                if xeno_gene in xeno_gene2clusters:
                    for clnum in xeno_gene2clusters[xeno_gene]:
                        xeno_cluster = Cluster(xeno_species, clnum)
                        bioclusters.append(xeno_cluster)
#                        logging.debug('\txeno_gene %s belongs to %s', xeno_gene, xeno_cluster)
        return bioclusters


    def calculate_links(self, mp, genomes):
        '''
        Count the number of unique orthologoues gene pairs ("links") between
        clusters in this pair.
        mp: QuickParanoid object.
        genomes: dict of all known genomes.
        'which': which genome/cluster is the basis for calculation;
                 result will change depending on this!
        '''
        links1 = 0
        links2 = 0
        logging.debug("Looking at %s genes in cluster %s (2nd is %s)...",
                      len(genomes[self.g1].cluster2genes[self.c1]),
                      self.c1, self.gc2)
        for gene in genomes[self.g1].cluster2genes[self.c1]: # FIXME: too slow
            # 'gene' is a 'locus_tag:genome_id' string.
            # TODO: convert gene from 'locus_tag:genome_id' to (locus_tag, genome_id) in:
            #  - Genome.cluster2genes
            #  - ClusterPair calculate_links
            # FIXME: optimize get_gene_links_to_bioclusters so that it stops as soon as self.gc2 is found! same below
            # FIXME: get_gene_links_to_bioclusters should return either 0 or 1, depending on whether clusters share a link.
            if self.gc2 in self.get_gene_links_to_bioclusters(gene, mp, genomes):
                links1 += 1
        logging.debug("Looking at %s genes in cluster %s (1st is %s)...",
                      len(genomes[self.g2].cluster2genes[self.c2]),
                      self.c2, self.gc1)
        for gene in genomes[self.g2].cluster2genes[self.c2]: # FIXME: too slow
            if self.gc1 in self.get_gene_links_to_bioclusters(gene, mp, genomes):
                links2 += 1
        return links1, links2


    def assign_orthologous_link(self, mp, genomes, args):
        '''
        Assign orthologous link to this cluster pair.
        '''
        link1, link2 = self.calculate_links(mp, genomes)
        # Make sure the number of links never exceeds the number of genes.
        genes_1 = self.num_c1_genes(genomes[self.g1])
        genes_2 = self.num_c2_genes(genomes[self.g2])
        link1 = min(link1, genes_1)
        link2 = min(link2, genes_2)
        logging.debug('\tlink1= %s and link2= %s for %s and %s, %s/%s genes',
                      link1, link2, self.gc1, self.gc2, genes_1, genes_2)
        if args.strict:
            # No link can be larger than the number of genes in the 2nd cluster.
            if link1 > genes_2:
                link1 = genes_2
                logging.debug('strict mode, new link1 is %s', link1)
            if link2 > genes_1:
                link2 = genes_1
                logging.debug('strict mode, new link2 is %s', link2)
            try:
                assert link1 <= min(genes_1, genes_2) and link2 <= min(genes_1,
                                                                       genes_2)
            except:
                logging.exception('c1: %s ; c2: %s', self.gc1, self.gc2)
                logging.exception('c1_genes: %s (c1: %s)', genes_1,
                                  genomes[self.g1].cluster2genes[self.c1])
                logging.exception('c2_genes: %s (c2: %s)', genes_2,
                                  genomes[self.g2].cluster2genes[self.c2])
                logging.exception('link1: %s ; link2: %s ; c1_genes: %s ; c2_genes: %s',
                                  link1, link2, genes_1, genes_2)
                raise
        self.link1 = link1
        self.link2 = link2


    def pre_CDS_identities(self, g1, g2):
        '''
        Prepare for running usearch: generate/open interlaced temporary file,
        return it.
        g1 and g2 are genomes[g1] and genomes[g2], respectively.
        '''
        # Get genelists for both clusters.
        gl1 = g1.cluster2genes[self.c1]
        gl2 = g2.cluster2genes[self.c2]
        # Make /tmp file and open it for writing.
        seqfile = NamedTemporaryFile(mode='w', dir='/tmp', delete = False)
        # Iterate all pairs of genes of this cluster pair, write to file.
        for gene1, gene2 in itertools.product(gl1, gl2):
            #logging.debug('gene1, gene2: %s, %s', gene1, gene2)
            seq1 = g1.get_protein(gene1.split(':')[0])
            seqfile.write(">%s\n%s\n" % (gene1, seq1))
            seq2 = g2.get_protein(gene2.split(':')[0])
            seqfile.write(">%s\n%s\n" % (gene2, seq2))
            assert len(seq1) > 0 and len(seq2) > 0
        seqfile.close()
        logging.debug('Created, opened, written to and closed seqfile %s.',
                      seqfile.name)
        return seqfile.name


    def gene_order(self, genome1, genome2):
        '''
        Are similar genes in the same order or not?
        '''
        # For clusters 1 and 2, for genes which have identity pairs,
        # build a list of per-genome gene indexes (i.e. 0, 1, 2...)
        # (from lower to higher cluster coordinate).
        inds = {self.gc1: [], self.gc2: []}
        # 'gene' is a GeneOrder namedtuple, (start, strand, locus_tag:genome_id).
        for gene_index, gene in enumerate(genome1.orderstrands[self.c1]): # iterate all cluster 1 genes
            if gene.geneid in self.gene1_to_gene2: # check if this gene has an identity pair
                # Offset indices by 1 (they start from 0, might be problematic for stats).
#                print(genome1.orderstrands[self.c1])
                inds[self.gc1].append(gene_index + 1)
                # Get the corresponding gene position from cluster 2.
                gene2_name = self.gene1_to_gene2[gene.geneid]
                gene2_feature = genome2.get_feature_by_locustag(gene2_name.split(':')[0]).feature
                gene2 = GeneOrder(int(gene2_feature.location.start.position), gene2_feature.strand, gene2_name)
                inds[self.gc2].append(genome2.orderstrands[self.c2].index(gene2) + 1)
#                print(genomes[self.g2].orderstrands[self.c2])
#        print(inds) # FIXME: debug only
        # spearman: monotonicity (same gene order, disregarding distances), less sensitive to outliers
        # pearson: linearity (same gene order AND similar distances); not really applicable: "a measure of the linear relationship between two continuous random variables"
        # kendall: unlike spearman and pearson, here each point contributes equally
        # magnitude: kendall < spearman
        self.pearson  = stats.pearsonr(inds[self.gc1], inds[self.gc2])[0]
        self.kendall  = stats.kendalltau(inds[self.gc1], inds[self.gc2])[0]
        self.spearman = stats.spearmanr(inds[self.gc1], inds[self.gc2])[0]
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
