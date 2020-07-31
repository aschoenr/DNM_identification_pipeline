#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import functions
import settings

def mark_adapters(uBAM):

    log_file = uBAM.split('.')[0] + '_markadapters.log'

    #cmd  = 'time (java -Xmx8G -jar ' + settings.picard + ' MarkIlluminaAdapters \\\n'
    cmd  = 'java -Xmx8G -jar ' + settings.picard + ' MarkIlluminaAdapters \\\n'
    cmd += '\tTMP_DIR=' + settings.tmp_dir + ' \\\n'
    #cmd += '\tNUM_PROCESSORS=' + str(num_threads) + ' \\\n'
    cmd += '\tI=' + uBAM + ' \\\n'
    cmd += '\tO=' + uBAM.split('.')[0] + '_markilluminaadapters.unmapped.bam \\\n'
    cmd += '\tM=' + uBAM.split('.')[0] + '_markilluminaadapters_metrics.txt \\\n'
    cmd += '\t>>' + log_file + ' 2>&1' #)'

    print cmd
    print ''
    os.system('echo \'' + cmd + '\n\'>' + log_file)
    os.system(cmd)

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir = settings.indirroot + '1_ubam/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(indir, '*markilluminaadapters*')

uBAMs = glob.glob(indir + '*.unmapped.bam')
for uBAM in uBAMs:
    if '_markilluminaadapters' in uBAM:
        uBAMs.remove(uBAM)
uBAMs.sort()

procs = []
status = []

for uBAM in uBAMs:
    proc = Process(target=mark_adapters, args=(uBAM,))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(indir, fam_id, procs, status)
elapsed_time = time.time() - start_time

print 'Elapsed time: ' + str(elapsed_time)
f = open(indir + 'time.log', 'a')
f.write(str(elapsed_time) + '\n')
f.close()
