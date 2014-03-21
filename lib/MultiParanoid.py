'''
See MultiParanoid class docs.
'''


from __future__ import print_function
import csv
import logging

# debug, info, warning, error, critical


class MultiParanoid(object):
    '''
    Incapsulates all access to the Multi/Quick-Paranoid results file.
    Parses the file once upon initialization.
    '''

    def __init__(self, path, no_tree_conflicts = False, no_name_conflicts = False):
        '''
        Initialize variables and then parse the file.
        path: to the multi/quick-paranoid sqltable file.
        no_tree_conficts: ignore lines where tree_conflict != 'No' (includes no_name_conflicts)
        no_name_conflicts: ignore lines where tree_conflict = 'diff. names' (less strict)
        '''
        self.path = path
        self.no_tree_conflicts = no_tree_conflicts
        self.no_name_conficts = no_name_conflicts
        # Full list of species found in the file.
        self.species = set()
        # Mapping of each gene to the orthology clusterID(s) it belongs to.
        self.gene2ortho = {}
        # List of all the genes in the orthology cluster with given clusterID.
        self.ortho2genes = {}
        self.parse()

    def get_species_list(self):
        '''
        Returns the list of species found during parsing.
        '''
        return list(self.species)

    def count_common_genes(self):
        '''
        Returns count of genes/orthologs common to (shared by) all the species in the
        multiparanoid file.
        '''
        return self.count_common_genes_subset(self.species)

    def count_common_genes_subset(self, species):
        """
        Counts genes/orthologs shared by the specified subset of species.
        'species' is a list, each item is the GenBank LOCUS of the respective genome
        (e.g. ['BA000030']).
        """
        # FIXME: numbers returned by this method differ from (are lower than) the original
        # InParanoid output files.
        # Number of genomes we are looking for.
        N = len(species)
        logging.info('We are looking for genes common to %s genomes.', N)
        # Counter of common genes.
        common = 0
        for cluster in self.ortho2genes.itervalues():
            if len(cluster) < N:
                logging.debug('cluster is too short, only %s genes', len(cluster))
                continue
            logging.debug('Cluster is long enough: %s genomes!', len(cluster))
            found_species = set()
            for gene in cluster:
                if self.gene2species(gene) in species:
                    found_species.add(self.gene2species(gene))
            if len(found_species) < N:
                continue
            else:
                common += 1
        return common

    def gene2species(self, gene):
        '''
        Given a complex gene ID like SAV_1680:BA000030, return the part after the
        colon (BA000030 in this example), which is the LOCUS of the gene's genome.
        '''
        return gene[gene.find(':')+1:]

    def parse(self):
        logging.info('Reading quick/multi-paranoid gene clusters:')
        with open(self.path) as tsv:
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
                    self.species.add(self.gene2species(row['gene']))
                    # Build gene-to-ortho-cluster-ID(s) dict. row['gene'] is the gene ID.
                    if self.no_tree_conflicts and row['tree_conflict'] != 'No':
                        continue
                    if self.no_name_conficts and row['tree_conflict'] == 'diff. names':
                        continue
                    if row['gene'] not in self.gene2ortho:
                        self.gene2ortho[row['gene']] = []
                    self.gene2ortho[row['gene']].append(row[idname])
                    # Build cluster-ID-to-all-genes dict of lists. row[idname] is the cluster number.
                    if row[idname] not in self.ortho2genes:
                        self.ortho2genes[row[idname]] = []
                    self.ortho2genes[row[idname]].append(row['gene'])
            except csv.Error as e:
                logging.exception('file %s, line %d: %s', self.path, reader.line_num, e)
                raise
        logging.info('\ttotal lines read: %s', reader.line_num)
        logging.info('\ttotal species encountered: %s', len(self.species))
        logging.info('\tlist of these species: %s', ', '.join(self.species))
        logging.info('\ttotal entries in gene2ortho: %s', len(self.gene2ortho))
        logging.info('\ttotal entries in ortho2genes: %s', len(self.ortho2genes))
