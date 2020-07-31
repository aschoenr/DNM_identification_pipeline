import sys
import os
import glob
import time
from multiprocessing import Process

import settings
import resource_usage


def check_args(args):
    if len(args) != 2:
        print '\n\tUSAGE: ' + args[0] + ' fam_id'
        print '\tNOTE: location of tools and root input/output dirs hardcoded in settings\n'
        sys.exit()

def check_indir(indir):
    if not os.path.isdir(indir):
        print '\n\tERROR: Invalid family ID (expecting ####_## from dir: ' + settings.indirroot + ')\n'
        sys.exit()

def check_outdir(outdir, file_descriptor='*'):

    if not os.path.isdir(outdir):
        os.system('mkdir ' + outdir)
    elif len(glob.glob(outdir + file_descriptor)) > 0:
        print '\n    WARNING: Output directory not empty (may overwrite existing data)'
        answer = ''

        while answer != 'y' and answer != 'n':
            answer = raw_input('    Do you want to continue (y/n): ')

        print ''

        if answer == 'n':
            sys.exit()

def all_done(l):
    for e in l:
        if not e:
            return False
    return True

def log_resource_use(outdir, fam_id, procs, status, jobs=None, func=None):
    logf = open(outdir + fam_id + '_res.log', 'a', 0)
    
    while(not all_done(status)):
        for i in range(len(procs)):
            p = procs[i]
            if p.is_alive():
                s = resource_usage.get_stats(p.pid)
                if s != '':
                    logf.write(s + '\n')
                    children = resource_usage.get_all_children(p.pid)
                    for c in children:
                        s = resource_usage.get_stats(p.pid)
                        if s != '':
                            logf.write(s + '\n')
            else:
                if status[i] == False:
                    status[i] = True
                    if jobs is not None and len(jobs) > 0:
                        job = jobs.pop(0)
                        #proc = Process(target=func, args=(job[0],job[1],job[2],job[3]))
                        proc = Process(target=func, args=tuple(job))
                        procs.append(proc)
                        status.append(False)
                        proc.start()
        time.sleep(60)
    logf.close()

    for proc in procs:
        proc.join()

def load_ped(ped_file):
    
    mother = ''
    father = ''
    daughters = []
    sons = []
    
    f = open(ped_file, 'r')
    for line in f:
        line = line.strip().split()
        if not(line[2] == '0' and line[3] == '0'):
            father = line[2]
            mother = line[3]
            if line[4] == '1':
                sons.append(line[1])
            else:
                daughters.append(line[1])
    f.close()

    return mother, father, daughters, sons
