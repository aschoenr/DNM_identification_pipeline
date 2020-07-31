#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions


def align(uBAM, indir, outdir, num_threads):

    sample_id = uBAM.split('_markilluminaadapters.unmapped.bam')[0].split(indir)[1]

    logfile = outdir + sample_id 

    stf_cmd  = 'java -Xmx8G -jar ' + settings.picard + ' SamToFastq \\\n'
    stf_cmd += '\tI=' + uBAM + ' \\\n'
    stf_cmd += '\tFASTQ=/dev/stdout \\\n'
    stf_cmd += '\tCLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \\\n'
    stf_cmd += '\tTMP_DIR=' + settings.tmp_dir + ' \\\n'
    stf_cmd += '\t2>>' + logfile + '_SamToFastq.log \\\n'

    bwa_cmd  = settings.bwa + ' mem -M -t ' + str(num_threads) + ' \\\n'
    bwa_cmd += '\t -p ' + settings.ref_genome + ' \\\n'
    bwa_cmd += '\t/dev/stdin \\\n'
    bwa_cmd += '\t2>>' + logfile + '_bwa.log \\\n'

    mba_cmd  = 'java -Xmx16G -jar ' + settings.picard + ' MergeBamAlignment \\\n'
    mba_cmd += '\tALIGNED_BAM=/dev/stdin \\\n'
    mba_cmd += '\tUNMAPPED_BAM=' + uBAM + ' \\\n'
    mba_cmd += '\tOUTPUT=' + outdir + sample_id + '.bam \\\n'
    mba_cmd += '\tR=' + settings.ref_genome + ' \\\n'
    mba_cmd += '\tCREATE_INDEX=true ADD_MATE_CIGAR=true \\\n'
    mba_cmd += '\tCLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \\\n'
    mba_cmd += '\tINCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \\\n'
    mba_cmd += '\tPRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \\\n'
    mba_cmd += '\tTMP_DIR=' + settings.tmp_dir + ' \\\n'
    mba_cmd += '\t2>>' + logfile + '_MergeBamAlignment.log '

    #cmd  = 'time (set -o pipefail && ' + stf_cmd + '\t| \\\n' + bwa_cmd + '\t| \\\n' + mba_cmd + ') '
    #cmd  = 'set -o pipefail && ' + stf_cmd + '\t| \\\n' + bwa_cmd + '\t| \\\n' + mba_cmd #+ ') '  
    cmd  = 'bash -c "set -o pipefail" && ' + stf_cmd + '\t| \\\n' + bwa_cmd + '\t| \\\n' + mba_cmd #+ ') ' 
    #cmd  = stf_cmd + '\t| \\\n' + bwa_cmd + '\t| \\\n' + mba_cmd # 2019-10-18 removed pipefail, wasn't working for some reason -Andrew

    print cmd
    print ''
    os.system('echo \'' + cmd + '\' > ' + outdir + sample_id + '_cmd.log')
    os.system(cmd)

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '1_ubam/' + fam_id + '/'
outdir = settings.indirroot + '2_align/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir, '*.bam')

uBAMs = glob.glob(indir + '*_markilluminaadapters.unmapped.bam')
uBAMs.sort()
num_threads = settings.total_cores / len(uBAMs)

procs = []
status = []

for uBAM in uBAMs:
    proc = Process(target=align, args=(uBAM, indir, outdir, num_threads,))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
