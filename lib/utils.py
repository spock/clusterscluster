import subprocess
import logging


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
    else:
        return out

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
    >>> test[k2] # get value using the equivalent k2 key
    3: 'ab'
    >>> test[k4] = 'ca' # set value using k4 key
    >>> test[k3] # get value using k3 key
    4: 'ca'
    >>> test
    5: {('a', 'b'): 'ab', ('a', 'c'): 'ca'}
    >>> test[k2] = 'ba' # this overwrites the stored value
    >>> test[k1] # proof
    'ba'
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
