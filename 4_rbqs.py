#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions


def rbqs(BAM, indir, outdir, num_threads):

    sample_id = BAM.split('_dedup.bam')[0].split(indir)[1]
    logfile = outdir + sample_id

    brcmd  = settings.gatk + ' BaseRecalibrator \\\n'
    #brcmd += '\t--reference ' + settings.ref_g_uncomp + ' \\\n'
    brcmd += '\t--reference ' + settings.ref_genome + ' \\\n'
    brcmd += '\t--input ' + BAM + ' \\\n'
    #brcmd += '\t--TMP_DIR ' + settings.tmp_dir + ' \\\n'
    #brcmd += '\t--known-sites ' + settings.dbsnp + ' \\\n'
    brcmd += '\t--known-sites ' + settings.dbsnp_rbqs + ' \\\n'
    brcmd += '\t--known-sites ' + settings.mills_indels + ' \\\n' 
    brcmd += '\t--known-sites ' + settings.r1000G + ' \\\n' 
    brcmd += '\t--output ' + outdir + sample_id + '_recal_data.table \\\n'
    #brcmd += '\t-nct ' + str(num_threads) + ' ' #CHECK

    #cmd  = 'time (' + brcmd + '\t>>' + logfile + '_build_model.log 2>&1) '
    cmd  = brcmd + '\t>>' + logfile + '_build_model.log 2>&1 ' #) '


    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile + '_build_model.log')
    os.system(cmd)


    brcmd  = settings.gatk + ' ApplyBQSR \\\n'
    #brcmd += '\t--reference ' + settings.ref_g_uncomp + ' \\\n'
    brcmd += '\t--reference ' + settings.ref_genome + ' \\\n'
    brcmd += '\t--input ' + BAM + ' \\\n'
    brcmd += '\t-bqsr ' + outdir + sample_id + '_recal_data.table \\\n'
    brcmd += '\t--output ' + outdir + sample_id + '_recal_reads.bam \\\n'
    brcmd += '\t--static-quantized-quals 10 \\\n'
    brcmd += '\t--static-quantized-quals 20 \\\n'
    brcmd += '\t--static-quantized-quals 30 \\\n'
    brcmd += '\t--add-output-sam-program-record \\\n'
    brcmd += '\t--create-output-bam-md5 \\\n'
    brcmd += '\t--use-original-qualities \\\n'
    #brcmd += '\t-nct ' + str(num_threads) + ' ' #CHECK

    #cmd  = 'time (' + brcmd + '\t>>' + logfile + '_apply_BQSR.log 2>&1) '
    cmd  = brcmd + '\t>>' + logfile + '_apply_BQSR.log 2>&1 '#) '

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile + '_apply_BQSR.log')
    os.system(cmd)


start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '3_mark_dups/' + fam_id + '/'
outdir = settings.indirroot + '4_rbqs/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir, '*.bam')

BAMs = glob.glob(indir + '*.bam')
BAMs.sort()

procs = []
status = []
num_threads = settings.total_cores / len(BAMs)

for BAM in BAMs:
    proc = Process(target=rbqs, args=(BAM, indir, outdir, num_threads, ))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
