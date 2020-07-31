total_cores    = 30
num_hc_procs   = 10     #haplotypecaller

#DIRECTORIES
app_dir        = '/dnm_pipeline/apps/'
ref_genome_dir = '/mnt/nfs1/danielle/hg38/'
indirroot      = '/mnt/nfs1/danielle/data/'
tmp_dir        = '/mnt/nfs1/danielle/temp/'
ped_dir        = '/mnt/nfs1/danielle/data/peds/'
hg38_lcr_dir   = '/mnt/nfs1/danielle/hg38_lcr/'


#APPS
picard         = app_dir + 'picard.jar'
bwa            = app_dir + 'bwa-0.7.17/bwa'
gatk           = app_dir + 'gatk-4.0.11.0/gatk'
gatk3          = app_dir + 'gatk-3.8-1-0/GenomeAnalysisTK.jar'
samtools       = 'samtools'
bcftools       = 'bcftools'
dng            = app_dir + 'denovogear/bin/dng'
vcf_subset     = 'vcf-subset'


#REFERENCE GENOME
ref_genome     = ref_genome_dir + 'Homo_sapiens_assembly38.fasta' #.gz'
ref_size       = 3209286105.0

##VQSR
hapmap         = ref_genome_dir + 'hapmap_3.3.hg38.vcf.gz' 
omni           = ref_genome_dir + '1000G_omni2.5.hg38.vcf.gz'
r1000G         = ref_genome_dir + '1000G_phase1.snps.high_confidence.hg38.vcf.gz'
r1000G_BO      = ref_genome_dir + '1000G_phase1.snps.high_confidence.hg38_BIALLELIC_ONLY.vcf.gz'
dbsnp_rbqs     = ref_genome_dir + 'dbsnp_146.hg38.vcf.gz'
dbsnp_rvqs     = ref_genome_dir + 'dbsnp_146.hg38.vcf.gz'
mills_indels   = ref_genome_dir + 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'

segmental_dups = hg38_lcr_dir + 'hg38.segmental_dups.bed'
seg_dups_tree  = hg38_lcr_dir + 'hg38.segmental_dups_interval_tree.bin'
simple_repeats = hg38_lcr_dir + 'hg38.simple_repeats.bed'
simp_reps_tree = hg38_lcr_dir + 'hg38.simple_repeats_interval_tree.bin'
repeatmasker   = hg38_lcr_dir + 'hg38.repeatmasker.bed'
remask_tree    = hg38_lcr_dir + 'hg38.repeatmasker_interval_tree.bin'
lcr_tree_file  = hg38_lcr_dir + 'hg38.interval_tree.bin'


# 1a_fastq_to_ubam.py
platform       = 'ILLUMINA'                                       
library        = 'WAKE_'

#FILTERING
parental_VAF_cutoff = 0.05
depth_cutoff        = 0.9999

segmental_dups      = hg38_lcr_dir + 'hg38.segmental_dups.bed'
seg_dups_tree       = hg38_lcr_dir + 'hg38.segmental_dups_interval_tree.bin'
simple_repeats      = hg38_lcr_dir + 'hg38.simple_repeats.bed'
simp_reps_tree      = hg38_lcr_dir + 'hg38.simple_repeats_interval_tree.bin'
repeatmasker        = hg38_lcr_dir + 'hg38.repeatmasker.bed'
remask_tree         = hg38_lcr_dir + 'hg38.repeatmasker_interval_tree.bin'
lcr_tree_file       = hg38_lcr_dir + 'hg38.interval_tree.bin'
