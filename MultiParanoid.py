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
        # Mapping of each gene to the orthology clusterID(s) it belongs to.
        self.gene2ortho = {}
        # List of all the genes in the orthology cluster with given clusterID.
        self.ortho2genes = {}

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
        print('\ttotal lines read:', reader.line_num)
        print('\ttotal entries in gene2ortho:', len(self.gene2ortho))
        print('\ttotal entries in ortho2genes:', len(self.ortho2genes))
