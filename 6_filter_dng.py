#!/usr/bin/python -u
import glob
import sys
import os
import time
import copy
import cPickle
import multiprocessing
from scipy.stats import poisson
from intervaltree import Interval, IntervalTree

import functions
import settings


def power_set(s):
    result = [[]]
    for elem in s:
        result.extend([x + [elem] for x in result])
    return result


def in_lcr(chr, loc, lcrs):

    for lcr in lcrs:
        if chr in lcr:    
            if lcr[chr].overlaps(loc):    
                return True
    return False 


#def generate_genome_stats(fam_id, id):
def generate_genome_stats(fam_id, kid_id, outdir):

    logfile = outdir + kid_id + '_ggs.log'

    ggscmd  = 'java -jar ' + settings.picard + ' CollectAlignmentSummaryMetrics \\\n'
    ggscmd += '\tR=' + settings.ref_genome + ' \\\n'
    ggscmd += '\tI=' + settings.indirroot + '4_rbqs/' + fam_id + '/' + kid_id + '_recal_reads.bam \\\n'
    ggscmd += '\tO=' + outdir + kid_id + '_genome_stats.txt \\\n' 

    #cmd  = 'time (' + ggscmd + '\t>>' + logfile + ' 2>&1)'
    cmd  = ggscmd + '\t>>' + logfile + ' 2>&1 '#)'

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile)
    os.system(cmd)


def generate_vcfs(fam_id, p_id, outdir):

    logfile = outdir + p_id + '_gvcf.log'

    gvcfcmd  = settings.samtools  + ' mpileup \\\n'
    gvcfcmd += '\t-l ' + outdir + p_id + '.bed \\\n' 
    gvcfcmd += '\t-t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -d10000000 -O -s \\\n'
    gvcfcmd += '\t--output-QNAME -gf ' + settings.ref_genome + ' \\\n' 
    gvcfcmd += '\t' + settings.indirroot + '4_rbqs/' + fam_id + '/' + p_id + '_recal_reads.bam \\\n'
    gvcfcmd += '\t2>>' + logfile + ' | \\\n'
    gvcfcmd += settings.bcftools + ' view > ' + outdir + p_id + '.vcf \\\n'
    gvcfcmd += '\t2>>' + logfile

    #cmd = 'time (' + gvcfcmd + ')'
    cmd = gvcfcmd

    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + logfile)
    os.system(cmd)

def generate_vcfs_no_bed(fam_id, p_id, outdir, dng_results, mom, dad):

    for r in dng_results:
            
        chr = r.split()[4]
        pos = r.split()[6]

        cmd  = settings.samtools  + ' mpileup -r ' + chr + ':' + pos + ' ' 
        cmd += '-t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -d10000000 -O -s '
        cmd += '--output-QNAME -gf ' + settings.ref_genome + ' '
        cmd += settings.indirroot + '4_rbqs/' + fam_id + '/' + p_id + '_recal_reads.bam | '
        cmd += settings.bcftools + ' view -H >> ' + outdir + p_id + '.vcf'

        print cmd
        os.system(cmd) 

        cmd  = settings.samtools  + ' mpileup -r ' + chr + ':' + pos + ' '
        cmd += '-t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -d10000000 -O -s '
        cmd += '--output-QNAME -gf ' + settings.ref_genome + ' '
        cmd += settings.indirroot + '4_rbqs/' + fam_id + '/' + mom + '_recal_reads.bam | '
        cmd += settings.bcftools + ' view -H >> ' + outdir + p_id + '_' + mom + '.vcf'

        print cmd
        os.system(cmd)

        cmd  = settings.samtools  + ' mpileup -r ' + chr + ':' + pos + ' '
        cmd += '-t DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -d10000000 -O -s '
        cmd += '--output-QNAME -gf ' + settings.ref_genome + ' '
        cmd += settings.indirroot + '4_rbqs/' + fam_id + '/' + dad + '_recal_reads.bam | '
        cmd += settings.bcftools + ' view -H >> ' + outdir + p_id + '_' + dad + '.vcf'


def load_vcfs(vcf_snp_results, vcf_indel_results, outdir, out_file=None):
    print '\nduplicate mpileup results'

    if out_file is not None:
        out_file.write('\nduplicate mpileup results')

    for p in vcf_snp_results:
        c = 0
        f = open(outdir + p + '.vcf', 'r')
        for line in f:
            if line[0] != '#':
                if 'INDEL' in line:
                    if vcf_snp_results[p][line.split()[0]].get(line.split()[1]) is not None:
                        c += 1 
                    if vcf_indel_results[p][line.split()[0]].get(line.split()[1]) is not None:
                        print '\nWTF\n' + p + ' ' + line.split()[0] + ' ' + line.split()[1] + '\n'
                    vcf_indel_results[p][line.split()[0]][line.split()[1]] = line.strip()            
                else:
                    if vcf_indel_results[p][line.split()[0]].get(line.split()[1]) is not None:
                        c += 1
                    if vcf_snp_results[p][line.split()[0]].get(line.split()[1]) is not None:
                        print '\nWTF\n' + p + ' ' + line.split()[0] + ' ' + line.split()[1] + '\n'
                    vcf_snp_results[p][line.split()[0]][line.split()[1]] = line.strip()
        f.close()
        print p + ' ' + str(c)
        if out_file is not None:
            out_file.write( p + ' ' + str(c))

def final_filter(p_id, dng_results, vcf_snp_results, vcf_indel_results, mom, dad, outdir, result_counts):
    
    f_snp   = open(outdir + p_id + '_snps.dnm', 'w')
    f_indel = open(outdir + p_id + '_indels.dnm', 'w')
    
    for r in dng_results[p_id]:
    
        indel = r.startswith('DENOVO-INDEL')

        cols = r.strip().split()
        chr  = cols[4]
        pos  = cols[6]
        #ref  = cols[8]

        full_vcf_data = ''
        if indel:
            full_vcf_data = vcf_indel_results[p_id][chr].get(pos,'')
            if full_vcf_data == '':
                full_vcf_data = vcf_snp_results[p_id][chr].get(pos,'')
        else:
            full_vcf_data = vcf_snp_results[p_id][chr].get(pos,'')
            if full_vcf_data == '':
                full_vcf_data = vcf_indel_results[p_id][chr].get(pos,'')

        vcf_data = full_vcf_data.strip().split('\t')
        data = vcf_data[-1].split(':')
        alts = vcf_data[4].split(',')
        ADFs = data[3].split(',')
        ADRs = data[4].split(',')

    
        alts_copy = copy.deepcopy(alts)
        good = False
        for i in range(len(alts_copy)):
            if int(ADFs[i+1]) > 1 and int(ADRs[i+1]) > 1:
                good = True
            else:
                alts.remove(alts_copy[i])
                
        if not good:
            result_counts[p_id][3] += 1
            continue
    
        if indel:
            mom_full_vcf_data = vcf_indel_results[mom][chr].get(pos,'')
            if mom_full_vcf_data == '':
                mom_full_vcf_data = vcf_snp_results[mom][chr].get(pos,'')
        else:
            mom_full_vcf_data = vcf_snp_results[mom][chr].get(pos,'')
            if mom_full_vcf_data == '':
                mom_full_vcf_data = vcf_indel_results[mom][chr].get(pos,'')

        mom_vcf_data = mom_full_vcf_data.strip().split('\t')
        mom_data = mom_vcf_data[-1].split(':')
        mom_alts = mom_vcf_data[4].split(',')
        mom_ADs  = mom_data[5].split(',')
        mom_DP   = mom_data[1]
    
        for i in range(len(mom_alts)):
            if mom_alts[i] not in alts:
                continue
            if float(mom_ADs[i+1])/float(mom_DP) > settings.parental_VAF_cutoff: #0.1:
                alts.remove(mom_alts[i])

        if len(alts) == 0:
            result_counts[p_id][4] += 1
            continue

        if indel:
            dad_full_vcf_data = vcf_indel_results[dad][chr].get(pos,'')
            if dad_full_vcf_data == '':
                dad_full_vcf_data = vcf_snp_results[dad][chr].get(pos,'')
        else:
            dad_full_vcf_data = vcf_snp_results[dad][chr].get(pos,'')
            if dad_full_vcf_data == '':
                dad_full_vcf_data = vcf_indel_results[dad][chr].get(pos,'')

        dad_vcf_data = dad_full_vcf_data.strip().split('\t')
        dad_data = dad_vcf_data[-1].split(':')
        dad_alts = dad_vcf_data[4].split(',')
        dad_ADs  = dad_data[5].split(',')
        dad_DP   = dad_data[1]
    
        for i in range(len(dad_alts)):
            if dad_alts[i] not in alts:
                continue
            if float(dad_ADs[i+1])/float(dad_DP) > settings.parental_VAF_cutoff: #0.1:
                alts.remove(dad_alts[i])

        if len(alts) == 0:
            result_counts[p_id][4] += 1
        else:
            if indel:
                f_indel.write(r.strip() + '\n')
                f_indel.write(full_vcf_data.strip() + '\n')
                f_indel.write(mom_full_vcf_data.strip() + '\n')
                f_indel.write(dad_full_vcf_data.strip() + '\n')
            else:
                f_snp.write(r.strip() + '\n')
                f_snp.write(full_vcf_data.strip() + '\n')
                f_snp.write(mom_full_vcf_data.strip() + '\n')
                f_snp.write(dad_full_vcf_data.strip() + '\n')

            result_counts[p_id][5] += 1

    f_snp.close()
    f_indel.close()

def create_interval_tree(data_file, output_file=None):
    f = open(data_file, 'r')
    int_tree = {}
    for i in range(1,23):
        chr = 'chr' + str(i)
        int_tree[chr] = IntervalTree()

    int_tree['chrX'] = IntervalTree()
    int_tree['chrY'] = IntervalTree()

    for line in f:
        line = line.strip().split()
        if line[0] in int_tree:
            int_tree[line[0]].addi(int(line[1]), int(line[2])+1, (int(line[1]), int(line[2])))
        #else:
        #    print '\nWTF\n' + str(line) + '\n'
        #    sys.exit()
    f.close()

    if output_file is not None:
        f = open(output_file, 'wb')
        cPickle.dump(int_tree, f, -1)
        f.close()

    return int_tree

def load_interval_tree(filename, l=None):

    print '\tloading ' + filename + '...'

    f = open(filename, 'rb')
    int_tree = cPickle.load(f)

    if l is not None:
        l.append(int_tree)
    
    f.close()
    return int_tree



start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot  + '5_dng/' + fam_id + '/'
outdir = settings.indirroot + '6_filtered_DNG/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)

summary_file = open(outdir + 'summary_' + str(time.time()) + '.txt', 'w')

summary_file.write('parental_VAF_cutoff = ' + str(settings.parental_VAF_cutoff) + '\n')
print 'parental_VAF_cutoff = ' + str(settings.parental_VAF_cutoff)

summary_file.write('depth cutff = ' + str(settings.depth_cutoff) + '\n\n')
print 'depth cutff = ' + str(settings.depth_cutoff) + '\n'

print 'loading low complexity regions...'
lcrs = []
'''
#lcrs.append(create_interval_tree(settings.segmental_dups, settings.seg_dups_tree))
#lcrs.append(create_interval_tree(settings.simple_repeats, settings.simp_reps_tree))
#lcrs.append(create_interval_tree(settings.repeatmasker, settings.remask_tree))
lcrs.append(load_interval_tree(settings.seg_dups_tree))
lcrs.append(load_interval_tree(settings.simp_reps_tree))
#lcrs.append(load_interval_tree(settings.remask_tree))
'''
manager = multiprocessing.Manager()
s_list = manager.list()
ps = []
summary_file.write('low complexity region files used for filtering:\n')
ps.append(multiprocessing.Process(target=load_interval_tree,args=(settings.seg_dups_tree,s_list,)))
summary_file.write('\t' + settings.seg_dups_tree + '\n')
ps.append(multiprocessing.Process(target=load_interval_tree,args=(settings.simp_reps_tree,s_list,)))
summary_file.write('\t' + settings.simp_reps_tree + '\n')
#ps.append(multiprocessing.Process(target=load_interval_tree,args=(settings.remask_tree,s_list,)))
#summary_file.write('\t' + settings.remask_tree + '\n')

for p in ps:
    p.start()

for p in ps:
    p.join()

for t in s_list:
    lcrs.append(t)

print str(len(lcrs)) + ' low complexity region files loaded into interval trees'


print '\nloading pedigree information'
    # load pedigree file to get family member IDs

ped_file = settings.ped_dir + fam_id + '.ped'
f = open(ped_file, 'r')
dad = ''
mom = ''
kids = []
dng_results = {}
vcf_snp_results = {}
vcf_indel_results = {}
coverage_cutoff = {}
f = open(ped_file, 'r')
for line in f:
    line = line.strip().split()
    dad = line[2]
    mom = line[3]
    kids.append(line[1])
    dng_results[line[1]] = []
    vcf_snp_results[line[1]] = {}
    vcf_indel_results[line[1]] = {}
    coverage_cutoff[line[1]] = 0
f.close()

### FOR NEWLY FORMATTED PED FILES WITH ENTRIES FOR MOM & DAD
kids.remove(mom)
kids.remove(dad)
dng_results.pop(mom, None)
dng_results.pop(dad, None)
vcf_snp_results.pop(mom, None)
vcf_snp_results.pop(dad, None)
vcf_indel_results.pop(mom, None)
vcf_indel_results.pop(dad, None)
coverage_cutoff.pop(mom, None)
coverage_cutoff.pop(dad, None)
####
#print kids
#sys.exit()  

vcf_snp_results[mom] = {}
vcf_snp_results[dad] = {}
vcf_indel_results[mom] = {}
vcf_indel_results[dad] = {}

    # generate genome summay statistics to find mean coverage
# UNCOMMENT
#'''
print 'generating genome summary statistics..'

procs = []
for k in kids:
    proc = multiprocessing.Process(target=generate_genome_stats, args=(fam_id, k, outdir,))
    procs.append(proc)
    #status.append(False)
    proc.start()

for proc in procs:
    proc.join()
#'''

    # load mean coverage from above

print '\nid\tmean_read_count\tread_cutoff'
summary_file.write('\nid\tmean_read_count\tread_cutoff\n')

for k in kids:
    f = open(outdir + k + '_genome_stats.txt', 'r')
    data = f.readlines()
    f.close()
    
    while len(data[-1].strip()) == 0:
        del data[-1]

    mean = float(data[-1].split()[7])/settings.ref_size #PF_ALIGNED_BASES
    coverage_cutoff[k] = poisson.ppf(settings.depth_cutoff, mean)
    print str(k) + '\t' + str(mean) + '\t' + str(coverage_cutoff[k])
    summary_file.write(str(k) + '\t' + str(mean) + '\t' + str(coverage_cutoff[k]) + '\n')
print ''
summary_file.write('\n')


chromosomes = []
for i in range(1,23):
    chr = 'chr' + str(i)
    chromosomes.append(chr)
    for k in kids:
        vcf_snp_results[k][chr] = {}
        vcf_indel_results[k][chr] = {}
    vcf_snp_results[mom][chr] = {}
    vcf_indel_results[mom][chr] = {}
    vcf_snp_results[dad][chr] = {}
    vcf_indel_results[dad][chr] = {}

chromosomes.extend(['chrX', 'chrY'])

result_counts = {}

for k in kids:
    vcf_snp_results[k]['chrX'] = {}
    vcf_snp_results[k]['chrY'] = {}
    vcf_indel_results[k]['chrX'] = {}
    vcf_indel_results[k]['chrY'] = {}
    result_counts[k] = [0,0,0,0,0,0]
vcf_snp_results[mom]['chrX'] = {}
vcf_snp_results[dad]['chrX'] = {}
vcf_snp_results[mom]['chrY'] = {}
vcf_snp_results[dad]['chrY'] = {}
vcf_indel_results[mom]['chrX'] = {}
vcf_indel_results[dad]['chrX'] = {}
vcf_indel_results[mom]['chrY'] = {}
vcf_indel_results[dad]['chrY'] = {}

#lcr_count = 0
#high_read_count = 0

for chr in chromosomes:
    print "filtering " + chr + '...'
    f = open(indir + fam_id + '_' + chr + '.dnm', 'r') 
    z = 0
    for line in f:
        id = line.split()[2]
        if id in kids:
            result_counts[id][0] += 1
            read_count = int(line.strip().split()[33])
            if read_count <= coverage_cutoff[id]:                   # filter high read counts
                if not in_lcr(chr, int(line.split()[6]), lcrs):     # filter low complexity region
                    dng_results[line.split()[2]].append(line.strip())
                else:
                    result_counts[id][2] += 1
            else:
                result_counts[id][1] += 1
    f.close()
print''

for k in kids:
    print k + '\t' + str(result_counts[k])

    # create BED files

mf = open(outdir + mom + '.bed', 'w')
df = open(outdir + dad + '.bed', 'w')

for k in kids:
    f = open(outdir + k + '.bed', 'w')
    for r in dng_results[k]:
        chr  = r.split()[4]
        pos2 = r.split()[6]
        pos  = str(int(pos2) -1)
        l = chr + '\t' + pos + '\t' + pos2 + '\n'
        f.write(l)
        mf.write(l)
        df.write(l)
    f.close()
mf.close()
df.close()

fam = [mom,dad]
fam.extend(kids)

#'''
    # generate VCFs
procs = []
for person in fam:
    proc = multiprocessing.Process(target=generate_vcfs, args=(fam_id, person, outdir,))
    procs.append(proc)
    #status.append(False)
    proc.start()

for proc in procs:
    proc.join()   
#'''

load_vcfs(vcf_snp_results, vcf_indel_results, outdir)

summary_file.write("id\ttotal_dnm_calls\thigh_read_count\tlow_complexity\tbidirectional\tparent_allele\n")

print "\nid\ttotal_dnm_calls\thigh_read_count\tlow_complexity\tbidirectional\tparent_allele"

for k in kids:
    final_filter(k, dng_results, vcf_snp_results, vcf_indel_results, mom, dad, outdir, result_counts)
    out = k
    for i in result_counts[k]:
        out += '\t' + str(i)
    print out
    summary_file.write(out + '\n')

summary_file.write('\n')
print''
kids.sort()

dng_results = {}
for k in kids:
    dng_results[k] = set([])
    f1 = open(outdir + k + '_snps.dnm', 'r')
    for line in f1:
        if not line.startswith('DENOVO'):
            line = line.strip().split('\t')
            dng_results[k].add(line[0] + '_' + line[1])
    f1.close()

    f1 = open(outdir + k + '_indels.dnm', 'r')
    for line in f1:
        if not line.startswith('DENOVO'):
            line = line.strip().split('\t')
            dng_results[k].add(line[0] + '_' + line[1])
    f1.close()



#non_unique = iopen(outdir + 'non_unique_DNMs.txt', 'w')

kids_powerset = power_set(kids)
kids_powerset.remove([])
for s in kids_powerset:
    r = set([])
    r = r | dng_results[s[0]]
    if len(s) > 1:
        for k in s[1:]:
            r = r & dng_results[k]
    print str(s) + '\t' + str(len(r))
    if len(s) > 1 and len(r) > 0:
        l = list(r)
        l.sort()
        summary_file.write(str(s) + '\t' + str(len(r)) + '\n' + str(l) + '\n\n')


elapsed_time = str(time.time() - start_time)

print "\nElapsed time: " + elapsed_time
summary_file.write('\nElapsed time: ' + elapsed_time + '\n')
summary_file.close()

time_file = open(outdir + 'time.log', 'w')
time_file.write(str(elapsed_time) + '\n')
time_file.close()
