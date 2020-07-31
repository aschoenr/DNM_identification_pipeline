#!/usr/bin/python -u
import os 

data_dir        = '/home/aschoenr/work/*'
ignore          = '0_fastq_raw'
script_dir      = '/home/aschoenr/ngs/scriptsAndResults/scripts'
backup_location = '/s2/andrew_data/'


cmd = 'rsync -av --size-only --omit-dir-times --no-g --no-o --no-perms --rsh="ssh -c arcfour -o Compression=no" --exclude ' + ignore + ' ' + data_dir + ' ' + backup_location
print cmd
os.system(cmd) 

cmd = 'rsync -av --size-only --omit-dir-times --no-g --no-o --no-perms --rsh="ssh -c arcfour -o Compression=no" --exclude ' + ignore + ' ' + script_dir + ' ' + backup_location
print cmd
os.system(cmd)
