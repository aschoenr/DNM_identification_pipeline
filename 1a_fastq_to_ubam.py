#!/usr/bin/python -u
import glob
import gzip
import os
import datetime
import sys
import time
from multiprocessing import Process

import settings
import functions


def convert(sample_dir, outdir):

    curr_dir = sample_dir 
    sample_id = curr_dir.split('/')[-1]                         #SM
    log_file = outdir + sample_id + '_fastq_to_ubam.log'

    f = gzip.open(curr_dir + '/' + sample_id + '_R1.fastq.gz', 'r')
    temp = f.readline().strip()
    flowcell_id = temp.split(':')[2]                            #ID
    sequencing_unit = flowcell_id + '.' + temp.split(':')[-1]   #PU
    f.close()

    this_library = settings.library + sample_id
    
    in1 = curr_dir + '/' + sample_id + '_R1.fastq.gz'
    in2 = curr_dir + '/' + sample_id + '_R2.fastq.gz'

    logfile = outdir + sample_id + '_FastqToSam.log'

    #cmd  = 'time (java -Xmx8G -jar ' + settings.picard + ' FastqToSam \\\n'
    cmd  = 'java -Xmx8G -jar ' + settings.picard + ' FastqToSam \\\n'
    cmd += '\tTMP_DIR=' + settings.tmp_dir + ' \\\n'         
    cmd += '\tFASTQ=' + in1 + ' \\\n'
    cmd += '\tFASTQ2=' + in2 + ' \\\n'
    cmd += '\tOUTPUT=' + outdir + sample_id + '.unmapped.bam \\\n'
    cmd += '\tREAD_GROUP_NAME=' + flowcell_id + ' \\\n'
    cmd += '\tSAMPLE_NAME=' + sample_id + ' \\\n'
    cmd += '\tLIBRARY_NAME=' + settings.library + sample_id + ' \\\n'
    cmd += '\tPLATFORM_UNIT=' + sequencing_unit + ' \\\n'
    cmd += '\tPLATFORM=' + settings.platform + ' \\\n'
    cmd += '\tSEQUENCING_CENTER=BI \\\n'
    cmd += '\tRUN_DATE=' + str(datetime.datetime.now().isoformat()) + ' \\\n'
    cmd += '\t>>' + logfile + ' 2>&1' #)'

    print cmd
    print ''
    os.system('echo \'' + cmd + '\' > ' + logfile)
    os.system(cmd)

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir = settings.indirroot + '0_fastq_raw/' + fam_id + '/'
outdir = settings.indirroot + '1_ubam/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)

sample_dirs = glob.glob(indir + '*')
sample_dirs.sort()
num_threads = settings.total_cores / len(sample_dirs)

procs = []
status = []

for sample_dir in sample_dirs:
    proc = Process(target=convert, args=(sample_dir,outdir,))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
