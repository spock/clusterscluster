from __future__ import print_function
import csv
import logging
from collections import defaultdict


def gene2species(gene):
    '''
    Given a complex gene ID like SAV_1680:BA000030, return tuple (gene, species) like
    (SAV_1680, BA000030), where SAV_1680 is a locus_tag, and BA000030 is a genome id.
    '''
    pos = gene.find(':')
    return (gene[:pos], gene[pos + 1:])


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
        self.gene2ortho = defaultdict(list)
        # List of all the genes in the orthology cluster with given clusterID.
        self.ortho2genes = defaultdict(list)
        # Dict mapping every gene to the set of all the other genes in the
        # same orthology cluster:
        # gene2genes[gene.genome][gene.id] = [other.gene1, other.gene2, ...]
        self.gene2genes = defaultdict(lambda: defaultdict(set))
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

    def parse(self):
        logging.info('Reading quick/multi-paranoid gene clusters:')
        with open(self.path) as tsv:
            # Sample line:
            # 1    avermitilis_MA4680.faa    SAV_1680.BA000030    1    1.000    Kitasatospora_setae_DSM43861.faa-SirexAA_E.faa-albus_J1074.faa-avermitilis_MA4680.faa-cattleya_DSM46488.faa-coelicolor_A3_2.faa-flavogriseus_IAF45CD.faa-griseus_NBRC13350.faa-scabiei_87.22.faa-venezuelae_Shinobu_719.faa-violaceusniger_Tu4113.faa    diff. numbers
            # fieldnames = ['#clusterID', 'species', 'gene', 'is_seed_ortholog', 'confidence', 'species_in_cluster', 'tree_conflict']
            reader = csv.DictReader(tsv, delimiter = '\t')
            try:
                for row in reader:
                    # FIXME: only do header check/assignment on the 1st row, not every row.
                    # Understand both Multi- and Quick-paranoid files.
                    if '#clusterID' in row:
                        # QuickParanoid
                        idname = '#clusterID'
                    else:
                        # MultiParanoid
                        idname = 'clusterID'
                    self.species.add(gene2species(row['gene'])[1])
                    # Skip some rows depending on arguments.
                    if self.no_tree_conflicts and row['tree_conflict'] != 'No':
                        continue
                    if self.no_name_conficts and row['tree_conflict'] == 'diff. names':
                        continue
                    # Build gene-to-ortho-cluster-ID(s) dict. row['gene'] is the gene ID.
                    self.gene2ortho[row['gene']].append(row[idname])
                    # Build cluster-ID-to-all-genes dict of lists. row[idname] is the cluster number.
                    self.ortho2genes[row[idname]].append(row['gene'])
            except csv.Error as e:
                logging.exception('file %s, line %d: %s', self.path, reader.line_num, e)
                raise
        logging.info('\tbuilding geneid-2-other-geneids dict...')
        for cluster in self.ortho2genes.itervalues():
            if len(cluster) < 2:
                # Do not create self-self clusters.
                logging.debug('skipping orthology cluster with 1 or 0 genes...')
                continue
            for gene in cluster:
                geneid, species = gene2species(gene)
                #logging.debug('gene2genes[%s] before adding: %s', species, self.gene2genes[species])
                setwogene = set(cluster)
                setwogene.remove(gene)
                # Use 'union' operator - allows adding entire lists, and is possibly faster than add.
                self.gene2genes[species][geneid] |= setwogene
                #logging.debug('Added %s to %s/%s.', setwogene, species, geneid)
                #logging.debug('gene2genes[%s] after adding: %s', species, self.gene2genes[species])
                #raise Exception('Oughooch')
        logging.info('\ttotal lines read: %s', reader.line_num)
        logging.info('\ttotal species encountered: %s', len(self.species))
        logging.info('\tlist of these species: %s', ', '.join(self.species))
        logging.info('\ttotal entries in gene2ortho: %s', len(self.gene2ortho))
        logging.info('\ttotal entries in ortho2genes: %s', len(self.ortho2genes))
        logging.info('\ttotal species in gene2genes: %s', len(self.gene2genes))
