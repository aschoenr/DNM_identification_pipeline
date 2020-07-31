#!/usr/bin/python

import sys
import os

if len(sys.argv) != 2:
    print '\n    USAGE: ./analyze_resource_usage.py resource_usage_log\n'
    sys.exit()

f = open(sys.argv[1], 'r')
start_time = int(f.readline().split()[1])
f.close()

max_ut  = 0
max_vm  = 0
max_rm  = 0
max_tmp = 0


procs = {}
times = []
c = 0                       # current time period
ct = start_time             # current time
current_time = []           # proc stats in current time period

f = open(sys.argv[1], 'r')
for line in f:
    line = line.strip().split()
    
    line[0] = int(line[0])
    line[1] = int(line[1]) - start_time
    line[2] = float(line[2])
    line[3] = int(line[3])
    line[4] = int(line[4])
    line[5] = int(line[5])

    if line[0] in procs:
        procs[line[0]].append(line)
    else:
        procs[line[0]] = [line]

    if (line[1] - ct < 10):
        current_time.append(line)
    else:
        times.append(current_time)
        current_time = [line]
        #c += 1        

    ct = line[1]

times.append(current_time)

#total_usage = []
c = 0
f = open('./temp.txt', 'w')
for time in times:
    t_cpu = 0
    t_vm = 0
    t_rm = 0
    tmp = 0
    
    for proc in time:
        t_cpu += proc[2]
        t_vm  += proc[3]
        t_rm  += proc[4]
        if tmp < proc[5]:
            tmp = proc[5]
    
    #total_usage.append[t_cpu, t_vm, t_rm, tmp]
    f.write(str(c) + '\t' + str(t_cpu) + '\t' + str(t_vm/1024/1024))
    f.write('\t' + str(t_rm/1024/1024) + '\t' + str(tmp) + '\n')
    c += 1
f.close() 

out = sys.argv[1].split('.log')[0] + '_cpu.pdf'
print out

f = open('./temp.plt', 'w')

f.write('set terminal post color enhanced\n')
f.write('set output "' + out + '"\n')
f.write('set ylabel "CPU usage (%)"\n')
f.write('set xlabel "Minutes"\n')
f.write('set style line 11 lc rgb "red" lt 1 lw 3\n')
f.write('set style line 31 lc rgb "blue" lt 1 lw 3\n')
f.write('plot "./temp.txt" using 1:2 with lines ls 11 title ""\n')

out = sys.argv[1].split('.log')[0] + '_ram.pdf'
print out
f.write('set output "' + out + '"\n')
f.write('set ylabel "Memory usage (GB)"\n')
f.write('plot "./temp.txt" using 1:3 with lines ls 31 title "Virtual memory", ')
f.write('"./temp.txt" using 1:4 with lines ls 11 title "Real memory"\n')

out = sys.argv[1].split('.log')[0] + '_tmpdir.pdf'
print out
f.write('set output "' + out + '"\n')
f.write('set ylabel "Temp Directory usage (KB)"\n')
f.write('plot "./temp.txt" using 1:3 with lines ls 11 title ""\n')


f.close()

os.system('gnuplot ./temp.plt')
os.system('rm temp.txt')
os.system('rm temp.plt')





