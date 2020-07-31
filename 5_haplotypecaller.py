#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions


def haplotypecaller(BAM, chr, indir, outdir):
  
    '''
     gatk --java-options "-Xmx4g" HaplotypeCaller  \
        -R Homo_sapiens_assembly38.fasta \
        -I input.bam \
        -O output.g.vcf.gz \
        -ERC GVCF
    '''

    sample_id = BAM.split('_recal_reads.bam')[0].split(indir)[1]

    logfile = outdir + sample_id + '_chr'  

    hccmd  = settings.gatk + ' --java-options "-Xmx4g" HaplotypeCaller \\\n'
    hccmd += '\t-R ' + settings.ref_genome + ' \\\n'
    if chr.endswith('*'):
        chr = chr[:-1]
        #hccmd += '\t-ploidy 1 \\\n'
    hccmd += '\t-L chr' + chr + ' \\\n'
    hccmd += '\t-I ' + BAM + ' \\\n'
    hccmd += '\t-O ' + outdir + sample_id + '_chr' + chr + '.g.vcf.gz -ERC GVCF \\\n'

    #cmd  = 'time (' + hccmd + '\t>>' + logfile + chr + '.log 2>&1) '
    cmd  = hccmd + '\t>>' + logfile + chr + '.log 2>&1 '#) '

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile + chr + '.log')
    os.system(cmd)


start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '4_rbqs/' + fam_id + '/'
outdir = settings.indirroot + '5_haplotypecaller/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)

BAMs = glob.glob(indir + '*.bam')
BAMs.sort()
num_threads = settings.total_cores / len(BAMs)

procs = []
status = []
males = set([])
mother, father, daughters, sons = functions.load_ped(settings.ped_dir + fam_id + '.ped')

males = set([father] + sons)

people = set([])
for BAM in BAMs:
    people.add(BAM)

jobs = []

for chr in range(1,23):
    for p in people:
        jobs.append([p, str(chr), indir, outdir])

for p in people:
    if p.split('_recal_reads.bam')[0].split(indir)[1] in males:
        jobs.append([p, 'X*', indir, outdir])
        jobs.append([p, 'Y*', indir, outdir])
    else:
        jobs.append([p, 'X', indir, outdir])


while len(jobs) > 0 and len(procs) < settings.num_hc_procs:
    job = jobs.pop(0)
    proc = Process(target=haplotypecaller, args=(job[0],job[1],job[2],job[3],))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status, jobs, haplotypecaller)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
