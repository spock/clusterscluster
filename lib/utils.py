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


def usearch(s1, s2, cutoff = 0.4, full = False):
    '''
    Run usearch on 2 sequences, return identity.
    s1, s2: two sequences to align
    cutoff: report zero for values below the cutoff
    full: perform full Dynamic Programming search (-fulldp)
    '''
    