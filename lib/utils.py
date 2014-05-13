import subprocess
import logging
from lib.ClusterPair import GP, Avg_Identity


# Method borrowed from antismash.utils
def execute(commands, inputs = None):
    "Execute commands in a system-independent manner"

    if inputs is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    proc = subprocess.Popen(commands, stdin=stdin_redir,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(input = inputs)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        logging.debug("%r %r returned %r" % (commands, inputs[:40], e))
        raise


def usearch(interleaved, cutoff = 0.4, full = False):
    '''
    Run usearch on interleaved sequences file, return output as text (3 columns:
    query ID, target ID, percent identity).
    interleaved: sequences to align, multifasta, 1st/2nd/1st/2nd/etc
    cutoff: report zero for values below the cutoff
    full: perform full Dynamic Programming search (-fulldp), usually improves
    percent identity.
    '''
    usearch = ['usearch']
    if full:
        usearch.append('-fulldp')
    usearch.extend(['-id', str(cutoff), '-pairs_global', interleaved, '-userfields',
                 'query+target+id', '-userout', '/dev/stdout', '-quiet'])
    out, err, retcode = execute(usearch)
    if retcode != 0:
        logging.error('usearch failed with %d: %r while searching %r, full output follows:\n%s',
                      retcode, err, interleaved, out)
        return None
    else:
        return out


def run_usearch(seqfilename, cutoff, fulldp):
    '''
    Wrapper for usearch(): run usearch() and parse/return results.
    seqfile: NamedTemporaryFile object
    cutoff: minimal identity for usearch
    fulldp: use full dynamic programming solution (~10x slower)
    Calculate per-gene-pair protein identity, and also average
    protein identity for cluster pair.
    Return value: avg_identity, gene1_to_gene2, protein_identities
    '''
    # Dict of genome1 genes mapped to genome2 genes (by locus_tag:genomeid),
    # only for genes which have above-cutoff identity.
    gene1_to_gene2 = {}
    # Gene-level protein identities, identities[(g1, g2)] = float,
    # symmetric (setting (g1, g2) makes (g2, g1) also accessible);
    # g1 and g2 are locus_tag:genome_id strings.
    protein_identities = SymKeyDict()
    gene_pairs = []
    results = usearch(seqfilename, cutoff, fulldp)
    if results == '':
        # Empty result: nothing above the cut-off, no similar gene pairs.
        return None
    if results == None:
        # Error in usearch().
        logging.exception('usearch failed, see output above')
        raise Exception('UsearchFailedError')
    # Parse results into a list of tuples, each tuple is 1 gene pair.
    results_list = results.strip('\n').split('\n')
    for row in results_list:
        try:
            # Query, target should be from g1, g2, respectively.
            query, target, identity = row.split('\t')
            gene_pairs.append(GP(identity, query, target))
        except ValueError:
            print('row:')
            print(row)
            logging.exception('Failed to parse usearch output.')
            raise Exception('UsearchOutputParseError')
    del results, results_list, identity, query, target
    # Sort gene pairs by identity, descending order.
    gene_pairs.sort(reverse = True)
    #print(gene_pairs)
    # Two sets to check that we have not yet seen genes from c1 and c2.
    seen_1 = set()
    seen_2 = set()
    # Counter of gene pairs.
    num_pairs = 0
    # Cumulative sum of identities, for average.
    sum_identities = 0.0
    for pair in gene_pairs:
        if pair.g1 in seen_1 or pair.g2 in seen_2:
            # at least one of the genes already has a pair
            #logging.debug('either %s or %s already has a pair', pair.g1, pair.g2)
            continue
        # else: save a new gene identity pair!
        gene1_to_gene2[pair.g1] = pair.g2
        protein_identities[(pair.g1, pair.g2)] = pair.identity
        seen_1.add(pair.g1)
        seen_2.add(pair.g2)
        num_pairs += 1
        sum_identities += float(pair.identity)
    avg_identity = Avg_Identity(num_pairs, sum_identities / num_pairs)
    del gene_pairs, seen_1, seen_2, sum_identities, num_pairs
    return avg_identity, gene1_to_gene2, protein_identities


class SymKeyDict(dict):
    '''
    "Symmetric key-tuple" dict. Should work for any n-tuple.
    Sorts the key-tuple before getting/setting values.
    Example:
    >>> test = SymKeyDict()
    >>> isinstance(test, dict)
    1: True
    >>> k1 = ('a', 'b')
    >>> k2 = ('b', 'a')
    >>> k3 = ('a', 'c')
    >>> k4 = ('c', 'a')
    >>> test[k1] = 'ab' # set value using key k1
    >>> test
    2: {('a', 'b'): 'ab'}
    >>> test[k2]        # get value using the equivalent k2 key
    3: 'ab'
    >>> test[k4] = 'ca' # set value using k4 key
    >>> test[k3]        # get value using k3 key
    4: 'ca'
    >>> test
    5: {('a', 'b'): 'ab', ('a', 'c'): 'ca'}
    >>> test[k2] = 'ba' # this overwrites the stored value
    >>> test[k1]        # proof
    'ba'
    >>> k2 in test      # though ('a', 'b') is stored as the key,
    True                # ('b', 'a') is also reported as a valid key
    '''

    def __init__(self, *arg, **kw):
        super(SymKeyDict, self).__init__(*arg, **kw)

    def _key_sorted(self, key):
        '''
        The only custom method here. Sorts and returns 'key'.
        '''
        return tuple(sorted(key))

    def __getitem__(self, key):
        return super(SymKeyDict, self).__getitem__(self._key_sorted(key))

    def __setitem__(self, key, value):
        super(SymKeyDict, self).__setitem__(self._key_sorted(key), value)

    def __delitem__(self, key):
        super(SymKeyDict, self).__delitem__(self._key_sorted(key))

    def __contains__(self, key):
        return super(SymKeyDict, self).__contains__(self._key_sorted(key))
