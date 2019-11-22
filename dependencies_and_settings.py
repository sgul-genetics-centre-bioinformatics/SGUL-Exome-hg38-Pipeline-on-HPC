#version 1.0

###################################
## SGUL Genetics Centre, 2019
## Dionysios Grigoriadis

## This python script provides the main dependencies and settings for running the
## hg38_WholeGenomeGermlineMultiSample.py pyton script and the
## hg38_WholeGenomeGermlineMultiSample.wdl workflow subsequently.

## The user can use this script to specify inputs and settings of the workflow.
###################################

###################################
#Main directories
maindir = "/storage/root/homes/dgrigoriadis/NGS/Exomes/hg38/"
drive = "/storage/root/homes/dgrigoriadis/resources/"
drive2 = drive + "Genome_reference_files/hg38/"
###################################

###################################
#Software
bwa = drive + "bwa/0.7.16a/bwa"
gatk = drive + "gatk/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar"
samtools = drive + "samtools/1.9/samtools"
sambamba = drive + "anaconda3/bin/sambamba"
verifybamid = drive + "/storage/root/homes/dgrigoriadis/resources/verifybamid/1.1.3/verifyBamID/bin/verifyBamID"
fastp = drive + "anaconda3/bin/fastp"
fastqc = drive + "anaconda3/bin/fastqc"
multiqc = drive + "anaconda3/bin/multiqc"
###################################

###################################
#References
contamination_sites_ud = drive2 + "Homo_sapiens_assembly38.contam.exome_calling_regions.v1.UD"
contamination_sites_bed = drive2 + "Homo_sapiens_assembly38.contam.exome_calling_regions.v1.bed"
contamination_sites_mu = drive2 + "Homo_sapiens_assembly38.contam.exome_calling_regions.v1.mu"
interval_list = drive2 + "exome_calling_regions.v1.interval_list"
ref_dict = drive2 + "Homo_sapiens_assembly38.dict"
ref_fasta = drive2 + "Homo_sapiens_assembly38.fasta"
ref_fasta_index = drive2 + "Homo_sapiens_assembly38.fasta.fai"
ref_alt = drive2 + "Homo_sapiens_assembly38.fasta.64.alt"
ref_sa = drive2 + "Homo_sapiens_assembly38.fasta.64.sa"
ref_amb = drive2 + "Homo_sapiens_assembly38.fasta.64.amb"
ref_bwt = drive2 + "Homo_sapiens_assembly38.fasta.64.bwt"
ref_ann = drive2 + "Homo_sapiens_assembly38.fasta.64.ann"
ref_pac = drive2 + "Homo_sapiens_assembly38.fasta.64.pac"

known_indels_sites_vcfs = [drive2 + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
                           drive2 + "Homo_sapiens_assembly38.known_indels.vcf.gz"]

known_indels_sites_indices = [drive2 + "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
                           drive2 + "Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"]

dbsnp_vcf = drive2 + "Homo_sapiens_assembly38.dbsnp138.vcf"
dbsnp_vcf_index = drive2 + "Homo_sapiens_assembly38.dbsnp138.vcf.idx"
###################################

###################################
#Settings
break_bands_at_multiples_of = 1000000
haplotype_scatter_count = 50
read_length = 150
###################################
