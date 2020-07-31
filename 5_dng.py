#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import functions 
import settings

def dng_dnm(fam_id, father, mother, kids, chr):

    # samtools mpileup -gDf hg19.fa s1.bam s2.bam s3.bam | dng dnm auto --ped sample.ped --bcf -
    # dngcmd  = samtools + ' mpileup -gDf ' + ref_genome_dir + 'Homo_sapiens_assembly38.fasta ' 

    ped_file = settings.ped_dir + fam_id + '.ped'

    dngcmd  = settings.samtools + ' mpileup -t DP -r ' + str(chr) + ' -gf \\\n'
    dngcmd += '\t' + settings.ref_genome + ' \\\n'
    for kid in kids:
        dngcmd += '\t' + indir + kid + '_recal_reads.bam \\\n'
    dngcmd += '\t' + indir + father + '_recal_reads.bam \\\n'
    dngcmd += '\t' + indir + mother + '_recal_reads.bam \\\n'
    dngcmd += '\t2>>' + outdir + fam_id + '_' + chr + '.log | \\\n' 

    dngcmd += settings.dng + ' dnm auto --ped ' + ped_file + ' --bcf - \\\n'
    dngcmd += '\t>'  + outdir + fam_id + '_' + chr + '.dnm \\\n'
    dngcmd += '\t2>>' + outdir + fam_id + '_' + chr + '.log' 

    #cmd = 'time (' + dngcmd + ')'
    cmd = dngcmd

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + outdir + fam_id + '_' + chr + '.log')
    os.system(cmd)


start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1] 

indir = settings.indirroot + '4_rbqs/' + fam_id + '/'
outdir = settings.indirroot + '5_dng/' + fam_id + '/'  

functions.check_indir(indir)
functions.check_outdir(outdir, '*.dnm')

mom, dad, daughters, sons = functions.load_ped(settings.ped_dir + fam_id + '.ped')
kids = daughters + sons

procs = []
status = []

'''
chromosomes = []
for i in range(1,23):
    chromosomes.append('chr' + str(i))
chromosomes.extend(['chrX', 'chrY'])

for chr in chromosomes: 
    proc = Process(target=dng_dnm, args=(fam_id, dad, mom, kids, chr,))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status)
'''

jobs = []
for i in range(1,23):
    jobs.append([fam_id, dad, mom, kids,'chr' + str(i)])
jobs.append([fam_id, dad, mom, kids,'chrX'])
jobs.append([fam_id, dad, mom, kids,'chrY'])

while len(jobs) > 0 and len(procs) < settings.total_cores:   #settings.num_hc_procs:
    job = jobs.pop(0)
    proc = Process(target=dng_dnm, args=tuple(job))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status, jobs, dng_dnm)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
