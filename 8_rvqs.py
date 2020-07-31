#!/usr/bin/python -u
import glob
import os
import sys
import time
from multiprocessing import Process

import settings
import functions

start_time = time.time()
print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
print sys.argv
functions.check_args(sys.argv)
fam_id = sys.argv[1]

indir  = settings.indirroot + '7_jointgenotype/' + fam_id + '/'
outdir = settings.indirroot + '8_rvqs/' + fam_id + '/'

functions.check_indir(indir)
functions.check_outdir(outdir)


logfile = outdir + fam_id + '_CV.log'

''' INCORRECT TOOL, SEE BELOW
java -jar GenomeAnalysisTK.jar \
    -T CombineVariants \
    -R reference.fasta \
    --variant input1.vcf \
    --variant input2.vcf \
    -o output.vcf \
    -genotypeMergeOptions PRIORITIZE
    -priority chr1,chr2,....
'''

vcfs = glob.glob(indir + '*.vcf')
vcfs.sort()

cvcmd  = 'java -jar ' + settings.gatk3 + ' \\\n '
cvcmd += '\t-T CombineVariants \\\n '
cvcmd += '\t-R ' + settings.ref_genome + ' \\\n '

for vcf in vcfs:
    v = vcf.split('.')[0].split('/')[-1]
    cvcmd += '\t--variant:' + v + ' ' + vcf + ' \\\n '
cvcmd += '\t-o ' + outdir + fam_id + '.vcf \\\n '
cvcmd += '\t-genotypeMergeOptions PRIORITIZE \\\n '
cvcmd += '\t-priority '
for vcf in vcfs:
    cvcmd += vcf.split('.')[0].split('/')[-1] + ','
cvcmd  = cvcmd[0:-1]
cvcmd += ' \\\n'

#cmd  = 'time (' + cvcmd + '\t>>' + logfile + ' 2>&1) '
cmd  = cvcmd + '\t>>' + logfile + ' 2>&1 '#) '

print cmd + '\n'
os.system('echo \'' + cmd + '\n\'>' + logfile)
os.system(cmd)
#sys.exit()

'''
java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
    -R reference.fasta \
    -V input1.vcf \
    -V input2.vcf \
    -out output.vcf \
    -assumeSorted


vcfs = glob.glob(indir + '*.vcf')
vcfs.sort()

cvcmd  = 'java -cp ' + settings.gatk3 + ' org.broadinstitute.gatk.tools.CatVariants \\\n '
cvcmd += '\t-R ' + settings.ref_genome + ' \\\n '

for vcf in vcfs:
    cvcmd += '\t-V ' + vcf + ' \\\n '
cvcmd += '\t--outputFile ' + outdir + fam_id + '.vcf \\\n '
cvcmd += '\t-assumeSorted \\\n'

cmd  = 'time (' + cvcmd + '\t>>' + logfile + ' 2>&1) '

print cmd + '\n'
os.system('echo \'' + cmd + '\n\'>' + logfile)
os.system(cmd)
'''

#def rvqs(chr):
 
'''
java -jar GenomeAnalysisTK.jar \ 
    -T VariantRecalibrator \ 
    -R reference.fa \ 
    -input raw_variants.vcf \ 
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf \ 
    -resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf \ 
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf \ 
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf \ 
    -an DP \ 
    -an QD \ 
    -an FS \ 
    -an SOR \ 
    -an MQ \
    -an MQRankSum \ 
    -an ReadPosRankSum \ 
    -an InbreedingCoeff \
    -mode SNP \ 
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ 
    -recalFile recalibrate_SNP.recal \ 
    -tranchesFile recalibrate_SNP.tranches \ 
    -rscriptFile recalibrate_SNP_plots.R  
'''

logfile = outdir + fam_id + '_VR_SNP.log'

vrcmd  = settings.gatk + ' VariantRecalibrator \\\n'
vrcmd += '\t-R ' + settings.ref_genome + ' \\\n '
vrcmd += '\t-V ' + outdir + fam_id + '.vcf \\\n'
#vrcmd += '-nt ' + str(settings.total_cores) + ' \\\n'
vrcmd += '\t--resource hapmap,known=false,training=true,truth=true,prior=15.0:'
vrcmd += settings.hapmap + ' \\\n'
vrcmd += '\t--resource omni,known=false,training=true,truth=true,prior=12.0:' 
vrcmd += settings.omni + ' \\\n'
vrcmd += '\t--resource 1000G,known=false,training=true,truth=false,prior=10.0:' 
vrcmd += settings.r1000G + ' \\\n'
vrcmd += '\t--resource dbsnp,known=true,training=false,truth=false,prior=2.0:' 
vrcmd += settings.dbsnp_rvqs + ' \\\n'
vrcmd += '\t-an DP \\\n'
vrcmd += '\t-an QD \\\n' 
vrcmd += '\t-an FS \\\n' 
vrcmd += '\t-an SOR \\\n'
vrcmd += '\t-an MQ \\\n'
vrcmd += '\t-an MQRankSum \\\n' 
vrcmd += '\t-an ReadPosRankSum \\\n' 
#vrcmd += -an InbreedingCoeff
vrcmd += '\t-mode SNP \\\n' 
vrcmd += '\t--max-gaussians 4 \\\n'
vrcmd += '\t-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n' 
vrcmd += '\t-O ' + outdir + fam_id + '_recalibrate_SNP.recal \\\n'
vrcmd += '\t--tranches-file ' + outdir + fam_id + '_recalibrate_SNP.tranches \\\n'
vrcmd += '\t--rscript-file ' + outdir + fam_id + '_recalibrate_SNP_plots.R \\\n'
#cmd  = 'time (' + vrcmd + '\t>>' + logfile + ' 2>&1) '
cmd  = vrcmd + '\t>>' + logfile + ' 2>&1 '#) '


print cmd + '\n'
os.system('echo \'' + cmd + '\n\'>' + logfile)
os.system(cmd)


'''
java -jar GenomeAnalysisTK.jar \ 
    -T ApplyRecalibration \ 
    -R reference.fa \ 
    -input raw_variants.vcf \ 
    -mode SNP \ 
    --ts_filter_level 99.0 \ 
    -recalFile recalibrate_SNP.recal \ 
    -tranchesFile recalibrate_SNP.tranches \ 
    -o recalibrated_snps_raw_indels.vcf 
'''

logfile = outdir + fam_id + '_AR_SNP.log'

arcmd  = settings.gatk + ' ApplyVQSR \\\n'
arcmd += '\t-R ' + settings.ref_genome + ' \\\n '
#arcmd += '\t--nt ' + str(settings.total_cores) + ' \\\n'
arcmd += '\t-V ' + outdir + fam_id + '.vcf \\\n'
arcmd += '\t-mode SNP \\\n'
arcmd += '\t--ts-filter-level 99.5 \\\n'
arcmd += '\t--recal-file ' + outdir + fam_id + '_recalibrate_SNP.recal \\\n'
arcmd += '\t--tranches-file ' + outdir + fam_id + '_recalibrate_SNP.tranches \\\n'
arcmd += '\t-O ' + outdir + fam_id + '_recalibrated_snps_raw_indels.vcf \\\n'
#cmd  = 'time (' + arcmd + '\t>>' + logfile + ' 2>&1) '
cmd  = arcmd + '\t>>' + logfile + ' 2>&1 '#) '

print cmd + '\n'
os.system('echo \'' + cmd + '\n\'>' + logfile)
os.system(cmd)


'''
java -jar GenomeAnalysisTK.jar \ 
    -T VariantRecalibrator \ 
    -R reference.fa \ 
    -input recalibrated_snps_raw_indels.vcf \ 
    -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.b37.vcf  \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.b37.vcf \
    -an QD \
    -an DP \ 
    -an FS \ 
    -an SOR \ 
    -an MQRankSum \ 
    -an ReadPosRankSum \ 
    -an InbreedingCoeff
    -mode INDEL \ 
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ 
    --maxGaussians 4 \ 
    -recalFile recalibrate_INDEL.recal \ 
    -tranchesFile recalibrate_INDEL.tranches \ 
    -rscriptFile recalibrate_INDEL_plots.R 
'''

logfile = outdir + fam_id + '_VR_INDEL.log'

vrcmd  = settings.gatk + ' VariantRecalibrator \\\n'
vrcmd += '\t-R ' + settings.ref_genome + ' \\\n '
vrcmd += '\t-V ' + outdir + fam_id + '_recalibrated_snps_raw_indels.vcf \\\n'
#vrcmd += '-nt ' + str(settings.total_cores) + ' \\\n'
#vrcmd += '\t--resource hapmap,known=false,training=true,truth=true,prior=15.0:'
#vrcmd += settings.hapmap + ' \\\n'
#vrcmd += '\t--resource omni,known=false,training=true,truth=true,prior=12.0:'
#vrcmd += settings.omni + ' \\\n'
vrcmd += '\t--resource mills,known=false,training=true,truth=true,prior=12.0:'
vrcmd += settings.mills_indels + ' \\\n'
vrcmd += '\t--resource dbsnp,known=true,training=false,truth=false,prior=2.0:'
vrcmd += settings.dbsnp_rvqs + ' \\\n'
vrcmd += '\t-an DP \\\n'
vrcmd += '\t-an QD \\\n'
vrcmd += '\t-an FS \\\n'
vrcmd += '\t-an SOR \\\n'
vrcmd += '\t-an MQ \\\n'
vrcmd += '\t-an MQRankSum \\\n'
vrcmd += '\t-an ReadPosRankSum \\\n'
#vrcmd += -an InbreedingCoeff
vrcmd += '\t-mode INDEL \\\n'
vrcmd += '\t-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n'
vrcmd += '\t--max-gaussians 4 \\\n'
vrcmd += '\t-O ' + outdir + fam_id + '_recalibrate_INDEL.recal \\\n'
vrcmd += '\t--tranches-file ' + outdir + fam_id + '_recalibrate_INDEL.tranches \\\n'
vrcmd += '\t--rscript-file ' + outdir + fam_id + '_recalibrate_INDEL_plots.R \\\n'
#cmd  = 'time (' + vrcmd + '\t>>' + logfile + ' 2>&1) '
cmd  = vrcmd + '\t>>' + logfile + ' 2>&1 '#) '

print cmd + '\n'
os.system('echo \'' + cmd + '\n\'>' + logfile)
os.system(cmd)


'''
java -jar GenomeAnalysisTK.jar \ 
    -T ApplyRecalibration \ 
    -R reference.fa \ 
    -input recalibrated_snps_raw_indels.vcf \ 
    -mode INDEL \ 
    --ts_filter_level 99.0 \ 
    -recalFile recalibrate_INDEL.recal \ 
    -tranchesFile recalibrate_INDEL.tranches \ 
    -o recalibrated_variants.vcf 
'''

logfile = outdir + fam_id + '_AR_INDEL.log'

arcmd  = settings.gatk + ' ApplyVQSR \\\n'
arcmd += '\t-R ' + settings.ref_genome + ' \\\n '
#arcmd += '\t--nt ' + str(settings.total_cores) + ' \\\n'
arcmd += '\t-V ' + outdir + fam_id + '_recalibrated_snps_raw_indels.vcf \\\n'
arcmd += '\t-mode INDEL \\\n'
arcmd += '\t--ts-filter-level 99.0 \\\n'
arcmd += '\t--recal-file ' + outdir + fam_id + '_recalibrate_INDEL.recal \\\n'
arcmd += '\t--tranches-file ' + outdir + fam_id + '_recalibrate_INDEL.tranches \\\n'
arcmd += '\t-O ' + outdir + fam_id + '_recalibrated_variants.vcf \\\n'
#cmd  = 'time (' + arcmd + '\t>>' + logfile + ' 2>&1) '
cmd  = arcmd + '\t>>' + logfile + ' 2>&1 '#) '

print cmd + '\n'
os.system('echo \'' + cmd + '\n\'>' + logfile)
os.system(cmd)

elapsed_time = time.time() - start_time
print 'Elapsed time: ' + str(elapsed_time)
f = open(outdir + 'time.log', 'w')
f.write(str(elapsed_time) + '\n')
f.close()

