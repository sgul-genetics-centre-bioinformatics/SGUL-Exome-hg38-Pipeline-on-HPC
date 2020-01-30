#version 1.0

###################################
## SGUL Genetics Centre, 2019
## Dionysios Grigoriadis

## This python script formats the directories and creates the inputs for running the
## hg38_WholeGenomeGermlineMultiSample.wdl workflow.
###################################
print("\n")
print("\n")
print("###################################################################################\n"
      "###########  Welcome to the PexOmes Pipelines by SGUL Genetics Centre  ############\n"
      "#                                                                                 #\n"
      "#This pipeline takes your raw sequencing data (fastq or bam file formats) and     #\n"
      "#performs Reads Pre-processing (trimming) - Various QC steps - Alignment -        #\n"
      "#Deduplication - BQSR and Variant Calling by implementing GATK4 Best Practices    #\n"
      "#for Exome Sequencing Processing                                                  #\n"
      "#---------------------------------------------------------------------------------#\n"
      "#Developers: Dionysios Grigoriadis, Dr Alan Pittman                               #\n"
      "#                                                                                 #\n"
      "#####     St George's University of London Genetics Centre, October 2019      #####\n"
      "###################################################################################\n")


#Parsing input
from optparse import OptionParser
import time
#User Input
#User Input
parser = OptionParser()
parser.add_option("-p", "--Project_name", dest="projectname",
                  help="Project name. Should be the name of the directory "
                       "that you keep the files with your BAM or FASTQ files e.g. apittman-tiny.")
# parser.add_option("-ws", "--Workflow_Start", dest="workflowstart",
#                  help="Optional: Workflow start: Where do you want to start your workflow from? Options:\n"
#                       "If you type alignment -> the workflow will run: alignment,markduplicates,baserecalibrator,variantcalling\n"
#                       "If you type markduplicates -> the workflow will run: markduplicates,baserecalibrator,variantcalling\n"
#                       "If you type baserecalibrator -> the workflow will run: baserecalibrator,variantcalling\n"
#                       "If you type variantcalling -> the workflow will run: variantcalling\n",
#                   default="alignment")
(options, args) = parser.parse_args()

#Import Necessary Packages
print("--INFO: Loading Packages:")
import os
import sys
import time
import glob
import datetime
import subprocess
import json
print("--INFO: Done!\n\n")
print("--INFO: Loading required packages:")
time.sleep(1)

#Loading Dependencies and Utilities
print("--INFO: Loading the dependencies_and_settings, tasks and utils:")
from dependencies_and_settings import *
from utils import *
from tasks import *
print("--INFO: Done!\n\n")
print("--INFO: Loading required packages:")
time.sleep(1)

#Processing and Cleaning Input
myproject = options.projectname
# workstart = options.workflowstart
#myproject = "apittman-test-samples-tiny3"
#workstart = "alignment"
mydir = maindir
now = datetime.datetime.now()

myproject = myproject.strip()
#workstart = workstart.strip()
#workstart = workstart.lower()

mydir = mydir.strip()

print("--INFO: Activating conda R-environment")
p = subprocess.Popen(["conda","activate","r-environment"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

print("--INFO: Making everything fine with the project directories:")
if mydir.endswith("/"):
    pass
else:
    mydir = mydir + "/"

#Formating directories
aligned_dir = mydir + "Aligned/"
fastqc_dir = mydir + "ReadsQC/"
fastq_raw_dir = mydir + "FASTQ_raw/"
gvcf_dir = mydir + "gVCF/"
runs_dir = mydir + "Runs/"
ubams_dir = mydir + "UBAMs/"
unaligned_dir = mydir + "Unaligned/"
vcf_dir = mydir + "VCF/"
tmp_dir = mydir + "tmp/"

#Make sure that the directories are correct
time.sleep(1)
while True:
    filelist = [aligned_dir, fastqc_dir, fastq_raw_dir, gvcf_dir, runs_dir, ubams_dir, unaligned_dir, vcf_dir]
    if all([os.path.isdir(f) for f in filelist]):
        print("        Base Directories Checked and OK\n")
        break
    else:
        print("        Base Directories Problem. Check spelling\n")
        break



print("\n")
time.sleep(1)
print("--INFO: Creating project-specific directories:")
#Creating project directories
aligned_project = aligned_dir + myproject
fastqc_project = fastqc_dir + myproject
fastq_raw_project = fastq_raw_dir + myproject
gvcf_project = gvcf_dir + myproject
runs_project = runs_dir + myproject
ubams_project = ubams_dir + myproject
unaligned_project = unaligned_dir + myproject
vcf_project = vcf_dir + myproject
time.sleep(1)
print("        Done! Example: "+aligned_project)
print("\n")

filelist = [aligned_project, fastqc_project, fastq_raw_project, gvcf_project,
            runs_project, ubams_project, unaligned_project, vcf_project]
for f in filelist:
    if not os.path.exists(f):
        os.makedirs(f)

#Defining project specific variables
initial_bams = glob.glob(ubams_project+"/*/*.bam")
initial_fqs = glob.glob(ubams_project+"/*/*.fq.gz") + \
              glob.glob(ubams_project+"/*/*.fq") + \
              glob.glob(ubams_project+"/*/*.fastq.gz") + \
              glob.glob(ubams_project+"/*/*.fastq")

if len(initial_bams)>=1 and len(initial_fqs)==0:
    print("--INFO: BAM files detected\n")
    dfiles = "bam"
    sample_names = [x.split("/")[-2] for x in initial_bams]
    indict = {}
    for sid in sample_names:
        indict[sid] = glob.glob(ubams_project+"/"+sid+"/"+"*R1*")+glob.glob(ubams_project+"/"+sid+"/"+"*.bam")
        if len(indict[sid])!=1:
            print("--ERROR Inputing Fastq files. Please check the R1 and R2 names.")
            raise()
    time.sleep(1)
    print("--INFO: Run settings:")
    print("        Time: " + str(now))
    print("        Project specified by the user: " + myproject)
    print("        Working directory: " + maindir)
    print("        Input BAM files: " + ", ".join(initial_bams))
    print("        Input sample_names: " + ", ".join(sample_names))
    print("        Number of BAM files: " + str(len(initial_bams)) + " files")
    print("\n")



if len(initial_bams)==0 and len(initial_fqs)>=1:
    print('--INFO: Fastq files detected\n')
    sample_names = list(set([x.split("/")[-2] for x in initial_fqs]))
    if len(sample_names)!=(len(initial_fqs)*2):
        print('--ERROR: Detected Less Fastqs than expected\n')
        raise()
    dfiles = "fastq"
    indict = {}
    for sid in sample_names:
        indict[sid] = glob.glob(ubams_project+"/"+sid+"/"+"*R1*")+glob.glob(ubams_project+"/"+sid+"/"+"*R2*")
        if len(indict[sid])!=2:
            print("--ERROR Inputing Fastq files. Please check the R1 and R2 names.")
            raise()
    time.sleep(1)
    print("--INFO: Run settings:")
    print("        Time: " + str(now))
    print("        Project specified by the user: " + myproject)
    print("        Working directory: " + maindir)
    print("        Input BAM files: " + ", ".join(initial_fqs))
    print("        Input sample_names: " + ", ".join(sample_names))
    print("        Number of BAM files: " + str(len(initial_fqs)) + " files")
    print("\n")



if len(initial_bams)>=1 and len(initial_fqs)>=1:
    print('--ERROR: Both BAMS and FastQ files detected.')
    raise()



for f in filelist:
    filelist2 = [f + "/" + x for x in sample_names]
    for f2 in filelist2:
        if not os.path.exists(f2):
            os.makedirs(f2)


time.sleep(1)
######################################################################################################
######################################################################################################
#Running the pipeline for each of the samples:
print("--INFO: Submitting the pipeline jobscripts: ")

#------------------------------------------------------------------------------------#
#GetBwaVersion
print("--INFO: Getting the bwa version: ")
getbwaversion = GetBwaVersion(bwa, path=runs_project+"/")
print(getbwaversion)
bv_file = runs_project+"/"+"GetBwaVersion/bwa_version.txt"
time.sleep(1)

while True:
    if hasRQJob(getbwaversion)==True:
        time.sleep(3)
        pass
    elif hasRQJob(getbwaversion)==False:
        time.sleep(1)
        break

with open(bv_file) as infile:
    bwa_version = infile.readlines()

bwa_version = bwa_version[0].strip('\n')
print("--INFO: Getting the bwa version: ")
print("        Bwa version captured: "+str(bwa_version)+"\n\n")
#------------------------------------------------------------------------------------#

print("--INFO: Starting Iterations:\n")

for samp_id in list(indict.keys()):

    sins = indict[samp_id]

    print("--INFO: Sample Name: "+samp_id)

    samHeader = "'@RG\\tID:" + samp_id + "\\tSM:" + samp_id + "\\tLB:" + samp_id + "\\tPL:ILLUMINA'"

    aligned_sample = aligned_project + "/" + samp_id + "/"
    fastqc_sample = fastqc_project + "/" + samp_id + "/"
    fastq_raw_sample = fastq_raw_project + "/" + samp_id + "/"
    gvcf_sample = gvcf_project + "/" + samp_id + "/"
    runs_sample = runs_project + "/" + samp_id + "/"
    ubams_sample = ubams_project + "/" + samp_id + "/"
    unaligned_sample = unaligned_project + "/" + samp_id + "/"
    vcf_sample = vcf_project + "/" + samp_id + "/"

    sam2indexed = SamToFastqAndBwaMemAndMba_FromBAM2(sample=samp_id,
                                       input_bam=sins[0],
                                       sam_header=samHeader,
                                       threads=8,
                                       output_fq=fastq_raw_sample+samp_id+".raw.fastq",
                                       trimmed_fq_basename=fastq_raw_sample+samp_id+".trimmed",
                                       report_basename =  fastq_raw_sample+samp_id,
                                       report_dir = fastq_raw_sample,
                                       output_aligned_sam = samp_id + ".aligned.unmerged.sam",
                                       output_bam_basename = aligned_sample+samp_id,
                                       path = runs_sample,
                                       jdep = check_last_job(getbwaversion),
                                       gatk=gatk, bwa=bwa, fastqc=fastqc, fastp=fastp, sambamba=sambamba,
                                       ref_fasta=ref_fasta, task="SamToFastqAndBwaMemAndMba_FromBAM2",
                                       cpu=8, mem="14gb",wtime="30:00:00")

    time.sleep(0.2)

    markdups = MarkDuplicates(sample=samp_id,
                              input_bam=aligned_sample+samp_id+".sorted.bam",
                              output_bam_basename=aligned_sample+samp_id+".unique",
                              output_prefix = aligned_sample+samp_id,
                              metrics_filename=aligned_sample+samp_id+".duplicate_metrics",
                              path = runs_sample,
                              jdep = check_last_job(sam2indexed),
                              task="MarkDuplicates",
                              cpu=1,mem="30gb",wtime="200:00:00")

    time.sleep(0.2)

    sortdup = SortSam(sample=samp_id,
                      input_bam=aligned_sample+samp_id+".unique.bam",
                      output_prefix = aligned_sample+samp_id,
                      output_bam_basename=aligned_sample+samp_id+".unique.sorted",
                      path = runs_sample,
                      jdep = check_last_job(markdups), gatk=gatk,
                      cpu=1,mem="5gb",wtime="50:00:00")

    time.sleep(0.2)

    baserec = BaseRecalibrator(sample=samp_id,
                                 input_bam=aligned_sample+samp_id+".unique.sorted.bam",
                                 recalibration_report_filename=aligned_sample+samp_id+".recal_data.csv",
                                 path=runs_sample,
                                 jdep = check_last_job(sortdup),
                                 task="BaseRecalibrator",
                                 dbsnp_vcf=dbsnp_vcf,
                                 known_indels_sites_vcfs=known_indels_sites_vcfs,
                                 ref_fasta=ref_fasta,
                                 gatk=gatk,
                                 cpu=1, mem="8gb",wtime="200:00:00")

    time.sleep(0.2)

    ancovar = AnalyzeCovariates(sample=samp_id,
                                recalibration_report_filename=aligned_sample+samp_id+".recal_data.csv",
                                plot_file=aligned_sample+samp_id+".recal_AnalyzeCovariates.pdf",
                                path = runs_sample,
                                jdep = check_last_job(baserec),
                                cpu=1, mem="4gb", wtime="5:00:00")

    time.sleep(0.2)

    applybqsr = ApplyBQSR(sample=samp_id,
                          input_bam=aligned_sample+samp_id+".unique.sorted.bam",
                          output_bam_basename=aligned_sample+samp_id+"_sorted_unique_recalibrated",
                          recalibration_report=aligned_sample+samp_id+".recal_data.csv",
                          output_prefix=aligned_sample+samp_id,
                          path = runs_sample,
                          jdep = check_last_job(baserec),
                          cpu=1,mem="5gb",wtime="100:00:00")

    time.sleep(0.2)

    wgs_metrics = CollectWgsMetrics(sample=samp_id,
                                    input_bam=aligned_sample+samp_id+"_sorted_unique_recalibrated.bam",
                                    metrics_filename=aligned_sample+samp_id+".wgs_metrics",
                                    path = runs_sample,
                                    jdep = check_last_job(applybqsr),
                                    cpu=1, mem="3gb", wtime="5:00:00")

    time.sleep(0.2)

    collect_agg = CollectAggregationMetrics(sample=samp_id,
                                            input_bam=aligned_sample+samp_id+"_sorted_unique_recalibrated.bam",
                                            output_bam_prefix=aligned_sample+samp_id,
                                            path=runs_sample,
                                            jdep=check_last_job(applybqsr),
                                            cpu=1, mem="7gb", wtime="20:00:00")

    time.sleep(0.2)

    hapcall = HaplotypeCaller(sample=samp_id,
                              input_bam=aligned_sample+samp_id+"_sorted_unique_recalibrated.bam",
                              base_file_name=gvcf_sample+samp_id,
                              bamout_base_name=aligned_sample+samp_id,
                              path=runs_sample,
                              jdep=check_last_job(applybqsr),
                              cpu=2, mem="20gb", wtime="200:00:00")

    time.sleep(0.2)

    valvcf = ValidateVCF(sample=samp_id,
                         input_vcf=gvcf_sample+samp_id+"g.vcf.gz",
                         path=runs_sample,
                         jdep=check_last_job(hapcall),
                         cpu=1, mem="4gb",wtime="20:00:00")

print("Done!")
