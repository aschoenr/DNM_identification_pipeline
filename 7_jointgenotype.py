#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions


'''
Parallelism options

This tool can be run in multi-threaded mode using this option.

    TreeReducible (-nt)
'''

def jointgenotype(chr, indir, outdir):

    '''
    gatk-launch GenotypeGVCFs \
        -R data/ref/ref.fasta \
        -V gendb://my_database \
        -G StandardAnnotation -newQual \
        -O test_output.vcf 
    '''

    tmp_chr = chr
    if chr.startswith('chr0'):
        tmp_chr = 'chr' + chr[4:]

    logfile = outdir + chr + '.log'
    jgcmd  = settings.gatk + ' GenotypeGVCFs \\\n'
    jgcmd += '\t-R ' + settings.ref_genome + ' \\\n'
    jgcmd += '\t-V gendb://' + indir + tmp_chr + ' \\\n'
    jgcmd += '\t-G StandardAnnotation -new-qual \\\n'
    jgcmd += '\t-O ' + outdir + chr + '.vcf \\\n\t'
    #jgcmd += '\t--num-threads ' + str(settings.total_cores) + ' \\\n\t'

    #cmd  = 'time (' + jgcmd + '>>' + logfile + ' 2>&1) '
    cmd  = jgcmd + '>>' + logfile + ' 2>&1 '#) '


    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile)
    os.system(cmd)

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '6_GenomicsDBImport/' + fam_id + '/'
outdir = settings.indirroot + '7_jointgenotype/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)

procs = []
status = []
jobs = []

for chr in range(1,23):
    if chr < 10:
        jobs.append(['chr0' + str(chr), indir, outdir])
    else:
        jobs.append(['chr' + str(chr), indir, outdir])

jobs.append(['chrX', indir, outdir])

if os.path.isdir(indir + 'chrY'):
    jobs.append(['chrY', indir, outdir])

while len(jobs) > 0 and len(procs) < settings.num_hc_procs:
    job = jobs.pop(0)
    proc = Process(target=jointgenotype, args=tuple(job))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status, jobs, jointgenotype)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
