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
    args = ['usearch']
    if full:
        args.append('-fulldp')
    args.extend(['-id', cutoff, '-pairs_global', interleaved, '-userfields',
                 'query+target+id,', '-userout', '/dev/stdout'])
    out, err, retcode = execute(args)
    if retcode != 0:
        logging.error('usearch failed with %d: %r while searching %r, full output follows:\n%s',
                      retcode, err, interleaved, out)
    else:
        return out
