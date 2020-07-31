#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions

def mark_dups(BAM, indir, outdir):
    
    sample_id = BAM.split('.bam')[0].split(indir)[1]
    logfile = outdir + sample_id

    mdcmd  = 'java -jar ' + settings.picard + ' MarkDuplicates \\\n' 
    mdcmd += '\tINPUT=' + BAM + ' \\\n' 
    mdcmd += '\tOUTPUT=' + outdir + sample_id + '_dedup.bam \\\n'
    mdcmd += '\tMETRICS_FILE=' + outdir + sample_id + '_metrics.txt \\\n'
    mdcmd += '\tTMP_DIR=' + settings.tmp_dir  

    #cmd  = 'time (' + mdcmd + '\\\n\t>>' + logfile + '_mark_dups.log 2>&1)'
    cmd  = mdcmd + '\\\n\t>>' + logfile + '_mark_dups.log 2>&1' #)'

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile + '_mark_dups.log')
    os.system(cmd)

    bbicmd  = 'java -jar ' + settings.picard + ' BuildBamIndex \\\n'
    bbicmd += '\tINPUT=' + outdir + sample_id + '_dedup.bam \\\n'
    bbicmd += '\tOUTPUT=' + outdir + sample_id + '_dedup.bai \\\n'
    bbicmd += '\tTMP_DIR=' + settings.tmp_dir   

    #cmd  = 'time (' + bbicmd + '\\\n\t>>' + logfile + '_build_index.log 2>&1)'
    cmd  = bbicmd + '\\\n\t>>' + logfile + '_build_index.log 2>&1' #)'

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile + '_build_index.log')
    os.system(cmd)


start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '2_align/' + fam_id + '/'
outdir = settings.indirroot + '3_mark_dups/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir, '*.bam')

BAMs = glob.glob(indir + '*.bam')
BAMs.sort()

procs = []
status = []

for BAM in BAMs:
    proc = Process(target=mark_dups, args=(BAM, indir, outdir,))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
