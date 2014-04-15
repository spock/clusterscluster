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
        self.links = 0
        # orthology weight
        self.weight = 0.0
        # is this cluster pair intra-species?
        if self.g1 == self.g2:
            self.intra = True
        else:
            self.intra = False
