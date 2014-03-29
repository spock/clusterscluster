#!/usr/bin/python


from __future__ import print_function
import sys
import os
import subprocess
from os import getcwd, chdir, remove
from os.path import join, dirname, realpath, isdir, exists
from multiprocessing import cpu_count, Queue, Process


#print(sys.argv)
if len(sys.argv) < 2:
    print("At least one file with inparanoid commands must be provided!")
    print("Example: %s commands_0.txt [commands_1.txt [...]]" % sys.argv[0])
    sys.exit(2)


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
        print("%r %r returned %r" % (commands, inputs[:40], e))
        raise


# Add custom_inparanoid to PATH.
custom_inparanoid = join(dirname(realpath(__file__)), 'custom_inparanoid')
os.environ['PATH'] = custom_inparanoid + ':' + os.environ['PATH']
#print('PATH after prepending custom_inparanoid: %s' % os.environ['PATH'])
del custom_inparanoid


# Check that 285_analysis/inparanoid exists and is a directory.
curr_path = getcwd()
analysis = join(dirname(realpath(__file__)), '285_analysis', 'inparanoid')
if not (exists(analysis) and isdir(analysis)):
    print("Cannot find the 285_analysis/inparanoid directory! Aborting!")
    sys.exit(3)


# Setup workers.
num_workers = cpu_count()
qsize = 10 * num_workers
tasks = Queue(maxsize = qsize)
# 1. Define worker.
def worker(tasks):
    'parallel worker for single command'
    chdir(analysis)
    while True:
        try:
            # args: a list with 1 or 2 elements
            (counter, args) = tasks.get() # By default, there is no timeout.
        except: # Queue.Empty
            print("worker(tasks) encountered an exception trying tasks.get()")
            raise
        # Exit if 'STOP' element is found.
        if args[0] == 'STOP':
            print("STOP found, exiting.")
            break
        blast = ['inparanoid.pl', '--blast-only']
        blast.extend(args)
        print('Start blast %s: %s' % (counter, ' '.join(blast)) )
        out, err, retcode = execute(blast)
        print(' Done blast %s: %s' % (counter, ' '.join(blast)) )
        if retcode != 0:
            print('inparanoid returned %d: %r while blasting %r, full output follows:\n%s' %
                  (retcode, err, ' and '.join(args), out) )
# 2. Start workers.
workers = []
for _ in range(num_workers):
    p = Process(target=worker, args=(tasks,))
    workers.append(p)
    p.start()
del _, p


# 3. Populate the tasks queue. Main loop.
counter = 0 # serial number of the current task
for infile in sys.argv[1:]:
    if not exists(infile):
        print("%s does not exist, skipping" % infile)
        continue
    print("Processing %s." % infile)
    with open(infile) as inh:
        # Change directory.
        chdir(analysis)
        for line in inh:
            line = line.strip()
            # parse
            fragments = line.split(' ')
            args = fragments[2:]
            # formatdb input
            for f in args:
                # first check if the formatdb file exists
                if not exists(f + '.psq'):
                    formatdb = ['formatdb', '-i', f]
                    out, err, retcode = execute(formatdb)
                    print('formatted %s: %s' % (f, ' '.join(formatdb)) )
                    if retcode != 0:
                        print('formatdb returned %d: %r while formatting %r, full output follows:\n%s' %
                              (retcode, err, f, out))
            # populate the queue
            counter += 1
            tasks.put((counter, args)) # will block until all-qsize items are consumed
        # change directory back to where we were
        chdir(curr_path)
    # delete the fully processed file
#    print("Deleting %s." % infile)
#    remove(infile)


# 4. Add STOP messages.
print("Adding %s STOP messages to task_queue." % num_workers)
for _ in range(num_workers):
    tasks.put((0, ['STOP']))
del _
# 5. Wait for all processes to finish, close the queue.
for p in workers:
    p.join()
tasks.close()

# TODO: cleanup all *.pin, *.psq, *.phr leftover files.
