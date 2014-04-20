from __future__ import print_function
import logging
import itertools
import csv
from tempfile import NamedTemporaryFile
from collections import namedtuple
from utils import usearch


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
            logging.debug('gene %s belongs to ortho-cluster %s', gene, orthoclust)
            for xeno_gene in mp.ortho2genes[orthoclust]:
                # Extract LOCUS from gene name, find species from it.
                xeno_species = xeno_gene.rsplit(':')[-1]
                # Check if xeno_gene belongs to any biosynthetic clusters.
                xeno_gene2clusters = genomes[xeno_species].gene2clusters
                if xeno_gene in xeno_gene2clusters:
                    for clnum in xeno_gene2clusters[xeno_gene]:
                        xeno_cluster = cluster(xeno_species, clnum)
                        bioclusters.extend(xeno_cluster)
                        logging.debug('\txeno_gene %s belongs to %s', xeno_gene,
                                      xeno_cluster)
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
                # FIXME: never finds anything
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
        protein identity for cluster pair.
        '''
        GP = namedtuple('GenePair', ['identity', 'g1', 'g2'])
        # Get genelists for both clusters.
        gl1 = genomes[self.g1].cluster2genes[self.c1]
        gl2 = genomes[self.g2].cluster2genes[self.c2]
        # Iterate all pairs of genes of this cluster pair.
        gene_pairs = []
        # Make /tmp file and open it for writing.
        with NamedTemporaryFile(mode='w', dir='/tmp') as seqfile:
            for gene1, gene2 in itertools.product(gl1, gl2):
                seq1 = genomes[self.g1].get_protein(gene1)
                seqfile.write(">%s\n%s\n" % (gene1, seq1))
                seq2 = genomes[self.g2].get_protein(gene2)
                seqfile.write(">%s\n%s\n" % (gene2, seq2))
                assert len(seq1) > 0 and len(seq2) > 0
            # Run usearch.
            results = usearch(seqfile)
            print(results)
            return
            # Use csv to parse results into a list of tuples
            reader = csv.reader([results])
            for row in reader:
                gene_pairs.append(GP(row[2], row[0], row[1])) # identity, query, target
        # Sort gene pairs by identity.
#        gene_pairs.sort()
#        for pair in pairs:
#            only allow 1, best similarity pair for each gene pair
#        if g1 in similarities, or g2 in similarities, (g1 in s and g2 in s[g1]) or (g2 in s and g1 in s[g2])
#            continue # similarity[g1][g2] = 0.6, symmetric
#        else:
#            similarities[g1] = {g2: similarity}
#            similarities[g2] = {g1: similarity}
