#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions

def genomicsDBImport(people, chr, indir,  outdir):

    '''
    gatk GenomicsDBImport \
        -V data/gvcfs/mother.g.vcf \
        -V data/gvcfs/father.g.vcf \
        -V data/gvcfs/son.g.vcf \
        --genomicsdb-workspace-path my_database \
        --intervals chr20
    '''

    logfile = outdir + chr

    gdbicmd  = settings.gatk + ' GenomicsDBImport \\\n'
    for p in people:
        gdbicmd += '\t-V ' + indir + p + '_' + chr + '.g.vcf.gz \\\n'

    gdbicmd += '\t--genomicsdb-workspace-path ' + outdir + chr + ' \\\n'
    gdbicmd += '\t--intervals ' + chr + ' \\\n\t'

    #cmd  = 'time (' + gdbicmd + ' 2>>' + logfile + '.log) '
    cmd  = gdbicmd + ' 2>>' + logfile + '.log '#) '

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile + '.log')
    os.system(cmd)

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '5_haplotypecaller/' + fam_id + '/'
outdir = settings.indirroot + '6_GenomicsDBImport/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)

procs = []
status = []
mother, father, daughters, sons = functions.load_ped(settings.ped_dir + fam_id + '.ped')
jobs = []

for chr in range(1,23):
    jobs.append([[father] + [mother] + sons + daughters, 'chr' + str(chr), indir, outdir])

jobs.append([[father] + [mother] + sons + daughters, 'chrX', indir, outdir])

#if len(sons) > 0:
jobs.append([[father] + sons,'chrY', indir, outdir])

while len(jobs) > 0 and len(procs) < settings.num_hc_procs:
    job = jobs.pop(0)
    proc = Process(target=genomicsDBImport, args=(job[0],job[1],job[2],job[3],))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status, jobs, genomicsDBImport)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()

