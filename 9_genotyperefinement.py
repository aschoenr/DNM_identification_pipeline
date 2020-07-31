#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions


def gtr(mom, dad, child, indir, outdir):
    
    #vcfcmd = 'cat ' + indir  + fam_id + '_recalibrated_variants.vcf | \\\n'
    vcfcmd = settings.vcf_subset + ' -c \\\n'
    vcfcmd += '\t' + dad + ',' + mom + ',' + child + ' \\\n'
    vcfcmd += '\t' + indir  + fam_id + '_recalibrated_variants.vcf \\\n'
    vcfcmd += '\t> ' + outdir + child + '.vcf \\\n'

    #cmd  = 'time (' + vcfcmd + '\t2>>' + outdir + child + '_vcf-subset.log) '
    cmd  = vcfcmd + '\t2>>' + outdir + child + '_vcf-subset.log '#) '
    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + outdir + child + '_vcf-subset.log')
    os.system(cmd)

    ped_file  = open(settings.ped_dir + fam_id + '.ped', 'r')
    temp_ped  = outdir + child + '.ped'
    temp_file = open(temp_ped, 'w')

    for line in ped_file:
        cols = line.split()
        if cols[1] in [mom, dad, child]:
            temp_file.write(line)
    ped_file.close()
    temp_file.close()



    '''
    java -jar GenomeAnalysisToolkit.jar 
       -R human_g1k_v37_decoy.fasta 
       -T CalculateGenotypePosteriors 
       --supporting 1000G_phase3_v4_20130502.sites.vcf 
       -ped trio.ped 
       -V recalibratedVariants.vcf 
       -o recalibratedVariants.postCGP.vcf
    '''

    cgpcmd  = settings.gatk  + ' CalculateGenotypePosteriors \\\n'
    cgpcmd += '\t-V ' + outdir + child + '.vcf \\\n'
    cgpcmd += '\t-O ' + outdir + child + '.postCGP.vcf \\\n'
    cgpcmd += '\t-R ' + settings.ref_genome + ' \\\n'
    cgpcmd += '\t--supporting ' + settings.r1000G_BO + ' \\\n'
    cgpcmd += '\t-ped ' + temp_ped + ' \\\n'

    #cmd  = 'time (' + cgpcmd + ' \t>>' + outdir + child + '_cgp.log 2>&1) '
    cmd  = cgpcmd + ' \t>>' + outdir + child + '_cgp.log 2>&1 '#) '
    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + outdir + child + '_cgp.log')
    os.system(cmd)


    '''
    java -jar $GATKjar 
        -T VariantFiltration 
        -R $bundlePath/b37/human_g1k_v37_decoy.fasta 
        -V recalibratedVariants.postCGP.vcf 
        -G_filter "GQ < 20.0" 
        -G_filterName lowGQ 
        -o recalibratedVariants.postCGP.Gfiltered.vcf
    '''

    vrcmd  = settings.gatk  + ' VariantFiltration \\\n'
    vrcmd += '\t-R ' + settings.ref_genome + ' \\\n'
    vrcmd += '\t-V ' + outdir + child + '.postCGP.vcf \\\n'
    vrcmd += '\t--G-filter "GQ < 20.0" \\\n'
    vrcmd += '\t--G-filter-name lowGQ \\\n'
    vrcmd += '\t-O ' + outdir + child + '.postCGP.Gfiltered.vcf \\\n'

    #cmd  = 'time (' + vrcmd + ' \t>>' + outdir + child + '_vf.log 2>&1) '
    cmd  = vrcmd + ' \t>>' + outdir + child + '_vf.log 2>&1 '#) '
    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + outdir + child + '_vf.log')
    os.system(cmd)


    '''
    java -jar $GATKjar 
        -T VariantAnnotator 
        -R $bundlePath/b37/human_g1k_v37_decoy.fasta 
        -V recalibratedVariants.postCGP.Gfiltered.vcf 
        -A PossibleDeNovo 
        -ped trio.ped 
        -o recalibratedVariants.postCGP.Gfiltered.deNovos.vcf
    '''

    #dnmcmd  = settings.gatk  + ' VariantAnnotator \\\n'
    dnmcmd  = 'java -jar ' + settings.gatk3  + ' -T VariantAnnotator \\\n'
    dnmcmd += '\t-R ' + settings.ref_genome + ' \\\n'
    dnmcmd += '\t-V ' + outdir + child + '.postCGP.Gfiltered.vcf \\\n'
    dnmcmd += '\t-A PossibleDeNovo \\\n'
    #dnmcmd += '\t--ped ' + settings.ped_dir + fam_id + '.ped \\\n'
    dnmcmd += '\t-ped ' + temp_ped + ' \\\n'
    #dnmcmd += '\t-O ' + outdir + fam_id + '_recalibratedVariants.postCGP.Gfiltered.DNMs.vcf \\\n\t'
    dnmcmd += '\t-o ' + outdir + child + '.postCGP.Gfiltered.DNMs.vcf \\\n'

    #cmd  = 'time (' + dnmcmd + ' \t>>' + outdir + child + '_dnm.log 2>&1) '
    cmd  = dnmcmd + ' \t>>' + outdir + child + '_dnm.log 2>&1 '#) '
    print cmd + '\n'
    os.system('echo \'' + cmd + '\n\'>' + outdir + child + '_dnm.log')
    os.system(cmd)

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '8_rvqs/' + fam_id + '/'
outdir = settings.indirroot + '9_genotyperefinement/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)

mom, dad, daughters, sons = functions.load_ped(settings.ped_dir + fam_id + '.ped')
kids = daughters + sons

procs = []
status = []

for kid in kids:
    proc = Process(target=gtr, args=(mom, dad, kid, indir, outdir))
    procs.append(proc)
    status.append(False)
    proc.start()

functions.log_resource_use(outdir, fam_id, procs, status)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()
