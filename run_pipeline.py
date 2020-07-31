#!/usr/bin/python
import datetime
import sys
import os

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print '\n    USAGE: ./run_pipeline.py fam_id start_stage\n'
    sys.exit()

fam_id = sys.argv[1]

stages = ['1a_fastq_to_ubam.py', '1b_MarkIlluminaAdapters.py', 
          '2_align.py', '3_mark_duplicates.py', '4_rbqs.py'] #'5_dng.py']

start = 0
if len(sys.argv) == 3:
    start = int(sys.argv[2])


print datetime.datetime.now().isoformat()


#for stage in stages:
for i in range(start, len(stages)):
    stage = stages[i]
    cmd = './' + stage + ' ' + fam_id
    print cmd
    os.system(cmd)
    print datetime.datetime.now().isoformat()
