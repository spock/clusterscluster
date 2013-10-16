'''
See MultiParanoid class docs.
'''


import csv
import logging

# debug, info, warning, error, critical


class MultiParanoid(object):
    '''
    Incapsulates all access to the Multi/Quick-Paranoid results file.
    Parses the file once upon initialization.
    '''

    def __init__(self, path):
        '''
        Initialize variables and then parse the file.
        '''
        self.path = path
        # Mapping of each gene to the orthology clusterID(s) it belongs to.
        self.gene2ortho = {}
        # List of all the genes in the orthology cluster with given clusterID.
        self.ortho2genes = {}

    def parse(self):
        if verbose > 0:
            print 'Reading quick/multi-paranoid gene clusters:'
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
                    if (ortho and row['tree_conflict'] != 'No') or (names and row['tree_conflict'] == 'diff. names'):
                        continue
                    if row['gene'] not in self.gene2ortho:
                        self.gene2ortho[row['gene']] = []
                    self.gene2ortho[row['gene']].append(row[idname])
                    # Build cluster-ID-to-all-genes dict of lists. row[idname] is the cluster number.
                    if row[idname] not in self.ortho2genes:
                        self.ortho2genes[row[idname]] = []
                    self.ortho2genes[row[idname]].append(row['gene'])
            except csv.Error as e:
                sys.exit('file %s, line %d: %s' % (self.path, reader.line_num, e))
        if verbose > 0:
            print '\ttotal lines read:', reader.line_num
            print '\ttotal entries in gene2ortho:', len(self.gene2ortho)
            print '\ttotal entries in ortho2genes:', len(self.ortho2genes)
