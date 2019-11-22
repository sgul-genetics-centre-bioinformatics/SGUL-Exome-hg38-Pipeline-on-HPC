from dependencies_and_settings import *
from utils import *
import os
import sys
import time
import glob
import datetime
import subprocess
import json


def InputBAMtoUnmappedBAM(sample, input_bam, output_unmapped_bam, path, gatk=gatk, cpu=1,mem="14gb",
                          wtime="15:00:00",task="InputBAMtoUnmappedBAM"):

    path = prepare_submission(path=path,task=task)

    scommand = "-Xmx8G -jar " + gatk +\
        "RevertSam " +\
        "--I "+input_bam +" "\
        "--O "+output_unmapped_bam+".bam "+\
        "--SANITIZE true "\
        "--MAX_DISCARD_FRACTION 0.005 "\
        "--ATTRIBUTE_TO_CLEAR XT "\
        "--ATTRIBUTE_TO_CLEAR XN "\
        "--ATTRIBUTE_TO_CLEAR AS "\
        "--ATTRIBUTE_TO_CLEAR OC "\
        "--ATTRIBUTE_TO_CLEAR OP "\
        "--SORT_ORDER queryname "\
        "--RESTORE_ORIGINAL_QUALITIES true "\
        "--REMOVE_DUPLICATE_INFORMATION true "\
        "--REMOVE_ALIGNMENT_INFORMATION true "\

    jobid = torque_submit(scommand, sample=sample, task=task, to_wait_id="",
                          wtime=wtime, cpu=cpu, mem=mem, cwd=path)

    return(jobid)



def GetBwaVersion(bwa, path, task="GetBwaVersion"):

    path = prepare_submission(path=path,task=task)

    scommand = bwa + " 2>&1 | grep -e '^Version' | sed 's/Version: //' > bwa_version.txt; sleep 10"

    jobid = torque_submit(scommand, sample="bwa", task=task, to_wait_id="",
                          wtime="00:05:00", nodes=1, cpu=1, mem="512mb", cwd=path)

    return(jobid)

def SamToFastqAndBwaMemAndMba_FromBAM1(sample, input_bam, bwa_commandline, bwa_version, output_fq, trimmed_fq_basename,
                               output_aligned_sam, output_bam_basename, report_basename, report_dir, path, jdep="",
                               gatk=gatk, bwa=bwa, fastqc=fastqc, fastp=fastp, ref_fasta=ref_fasta, task="SamToFastqAndBwaMemAndMba_FromBAM1",
                               cpu=8, mem="14gb",wtime="30:00:00"):

    path = prepare_submission(path=path,task=task)

    scommand = "bash_ref_fasta="+ref_fasta+";\n\n"+\
        "java -Xms5000m -jar "+gatk+\
        " SamToFastq "\
        "--INPUT "+input_bam+" "+\
        "--FASTQ "+output_fq+" "+\
        "--INTERLEAVE true --NON_PF true;\n\n"+\
        ""\
        ""+\
        fastqc+" "+output_fq+" -o "+report_dir+";\n\n"\
        ""\
        ""+\
        fastp+" -i "+output_fq+" --interleaved_in"+" "+\
        "-o "+trimmed_fq_basename+".R1.fastq.gz "+\
        "-O "+trimmed_fq_basename+".R2.fastq.gz "+\
        "--html "+report_basename+".after_fastp.html "+\
        "--json " + report_basename + ".after_fastp.json " +\
        "-w 10;\n\n"+\
        ""\
        ""+\
        fastqc+" "+\
        trimmed_fq_basename+".R1.fastq.gz "+ \
        trimmed_fq_basename+".R2.fastq.gz "+ \
        "-o "+report_dir+";\n\n"\
        ""\
        ""+\
        bwa + " " + bwa_commandline + " -t "+str(threads)+" "+\
        trimmed_fq_basename+".R1.fastq.gz "+\
        trimmed_fq_basename+".R2.fastq.gz "+\
        "> "+output_aligned_sam+" "\
        "2> >(tee "+output_bam_basename+".bwa.stderr.log >&2);\n\n"+\
        ""\
        ""+\
        "java -Xms3000m -jar "+gatk+" "+\
        "MergeBamAlignment "\
        "--VALIDATION_STRINGENCY SILENT "\
        "--EXPECTED_ORIENTATIONS FR "\
        "--ATTRIBUTES_TO_RETAIN X0 "\
        "--ATTRIBUTES_TO_REMOVE NM "\
        "--ATTRIBUTES_TO_REMOVE MD "\
        "--ALIGNED_BAM "+output_aligned_sam+" "+\
        "--UNMAPPED_BAM "+input_bam+" "+\
        "--OUTPUT "+output_bam_basename+".bam "+\
        "--REFERENCE_SEQUENCE "+ref_fasta+" "+\
        "--PAIRED_RUN true "\
        "--SORT_ORDER unsorted "\
        "--IS_BISULFITE_SEQUENCE false "\
        "--ALIGNED_READS_ONLY false "\
        "--CLIP_ADAPTERS false "\
        "--MAX_RECORDS_IN_RAM 2000000 "\
        "--ADD_MATE_CIGAR true "\
        "--MAX_INSERTIONS_OR_DELETIONS -1 "\
        "--PRIMARY_ALIGNMENT_STRATEGY MostDistant "\
        '--PROGRAM_RECORD_ID "bwamem" '\
        '--PROGRAM_GROUP_VERSION "'+bwa_version+'" '\
        '--PROGRAM_GROUP_COMMAND_LINE "'+bwa_commandline+' " '\
        '--PROGRAM_GROUP_NAME "bwamem" '\
        '--UNMAPPED_READ_STRATEGY COPY_TO_TAG '\
        '--ALIGNER_PROPER_PAIR_FLAGS true '\
        '--UNMAP_CONTAMINANT_READS true '\
        '--ADD_PG_TAG_TO_READS false;\n\n'\
        ""\
        'if test -f '+output_bam_basename+'.bam; then rm '+trimmed_fq_basename+".R*.fastq.gz "+output_fq+" "+output_aligned_sam+" fi;\n\n"\
        ""\
        'grep -m1 "read .* ALT contigs" '+output_bam_basename+'.bwa.stderr.log | grep -v "read 0 ALT contigs"'

    jobid = torque_submit(command=scommand, sample=sample, task=task,
                          cwd=path, to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)


def SamToFastqAndBwaMemAndMba_FromBAM2(sample, input_bam, sam_header, threads, output_fq, trimmed_fq_basename,
                               output_aligned_sam, output_bam_basename, report_basename, report_dir, path, jdep="",
                               gatk=gatk, bwa=bwa, fastqc=fastqc, fastp=fastp, sambamba=sambamba, ref_fasta=ref_fasta, task="SamToFastqAndBwaMemAndMba_FromBAM2",
                               cpu=8, mem="14gb",wtime="30:00:00"):

    path = prepare_submission(path=path,task=task)

    scommand = "bash_ref_fasta="+ref_fasta+"; "+\
        "java -Xms5000m -jar "+gatk+\
        " SamToFastq "\
        "--INPUT "+input_bam+" "+\
        "--FASTQ "+output_fq+" "+\
        "--INTERLEAVE true --NON_PF true;\n\n"+\
        ""\
        ""+\
        fastqc+" "+output_fq+" -o "+report_dir+";\n\n"\
        ""\
        ""+\
        fastp+" -i "+output_fq+" --interleaved_in"+" "+\
        "-o "+trimmed_fq_basename+".R1.fastq.gz "+\
        "-O "+trimmed_fq_basename+".R2.fastq.gz "+\
        "--html "+report_basename+".after_fastp.html "+\
        "--json " + report_basename + ".after_fastp.json " +\
        "-w 10;\n\n"+\
        ""\
        ""+\
        fastqc+" "+\
        trimmed_fq_basename+".R1.fastq.gz "+ \
        trimmed_fq_basename+".R2.fastq.gz "+ \
        "-o "+report_dir+";\n\n"\
        ""\
        ""+\
        bwa + " mem " +\
        ref_fasta + " " +\
        trimmed_fq_basename+".R1.fastq.gz "+\
        trimmed_fq_basename+".R2.fastq.gz "+\
        "-t "+str(threads)+" "+\
        "-R "+sam_header+" "+\
        "> "+output_aligned_sam+" "\
        "2> >(tee "+output_bam_basename+".bwa.stderr.log >&2);\n\n"+\
        ""\
        ""+\
        'echo "sambamba1";\n'+\
        sambamba + " view -S -f bam -t " + str(threads) + " " + output_aligned_sam + " > " + output_bam_basename+".bam;\n\n"+\
        'echo "sambamba2";\n'+\
        sambamba + " sort -m 5000000000 -t " + str(threads) + " " + output_bam_basename+".bam " + "-o " + output_bam_basename+".sorted.bam;\n\n"+\
        'echo "sambamba3";\n'+\
        sambamba + " index -t " + str(threads) + " " + output_bam_basename+".sorted.bam;\n\n"+\
        ""\
        'echo "remove";\n'+\
        'if test -f '+output_bam_basename+'.sorted.bam; then rm '+trimmed_fq_basename+".R*.fastq.gz "+output_fq+" "+output_aligned_sam+" "+output_bam_basename+".bam "+"\n\n"\
        'fi;\n\n'\
        ""\
        'echo "grep";\n'+\
        'grep -m1 "read .* ALT contigs" '+output_bam_basename+'.bwa.stderr.log | grep -v "read 0 ALT contigs"'

    jobid = torque_submit(command=scommand, sample=sample, task=task,
                          cwd=path, to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)

def SamToFastqAndBwaMemAndMba_FromFASTA_1(sample, input_fq, threads, sam_header, trimmed_fq_basename,
                               output_aligned_sam, output_bam_basename, report_basename, report_dir, path, jdep="",
                               bwa=bwa, fastqc=fastqc, fastp=fastp, sambamba=sambamba, ref_fasta=ref_fasta, task="SamToFastqAndBwaMemAndMba_FromFASTA_1",
                               cpu=8, mem="14gb",wtime="30:00:00"):

    path = prepare_submission(path=path,task=task)

    scommand = fastp+" -i "+input_fq+" --interleaved_in"+" "+\
        "-o "+trimmed_fq_basename+".R1.fastq.gz "+\
        "-O "+trimmed_fq_basename+".R2.fastq.gz "+\
        "--html "+report_basename+".after_fastp.html "+\
        "--json " + report_basename + ".after_fastp.json " +\
        "-w 10;\n\n"+\
        ""\
        ""+\
        fastqc+" "+\
        trimmed_fq_basename+".R1.fastq.gz "+ \
        trimmed_fq_basename+".R2.fastq.gz "+ \
        "-o "+report_dir+";\n\n"\
        ""\
        ""+\
        bwa + " mem " +\
        ref_fasta + " " +\
        trimmed_fq_basename+".R1.fastq.gz "+\
        trimmed_fq_basename+".R2.fastq.gz "+\
        "-t "+str(threads)+" "+\
        "-R "+sam_header+" "+\
        "> "+output_aligned_sam+" "\
        "2> >(tee "+output_bam_basename+".bwa.stderr.log >&2);\n\n"+\
        ""\
        ""+\
        sambamba + " view -S -f bam -t " + str(threads) + " " + output_aligned_sam + " > " + output_bam_basename+".bam;\n\n"+\
        sambamba + " sort -m 5000000000 -t " + str(threads) + " " + output_bam_basename+".bam " + "-o " + output_bam_basename+".sorted.bam;\n\n"+\
        sambamba + " index -t " + str(threads) + " " + output_bam_basename+".sorted.bam;\n\n"+\
        ""\
        'if test -f '+output_bam_basename+'.bam; then rm '+trimmed_fq_basename+".R*.fastq.gz "+output_fq+" "+output_aligned_sam+" fi;\n\n"\
        ""\
        'grep -m1 "read .* ALT contigs" '+output_bam_basename+'.bwa.stderr.log | grep -v "read 0 ALT contigs"'

    jobid = torque_submit(command=scommand, sample=sample, task=task,
                          cwd=path, to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)


def SamToFastqAndBwaMemAndMba_FromFASTA_2(sample, input_fq1, input_fq2, threads, sam_header, trimmed_fq_basename,
                               output_aligned_sam, output_bam_basename, report_basename, report_dir, path, jdep="",
                               bwa=bwa, fastqc=fastqc, fastp=fastp, sambamba=sambamba, ref_fasta=ref_fasta, task="SamToFastqAndBwaMemAndMba_FromFASTA_2",
                               cpu=8, mem="14gb",wtime="30:00:00"):

    path = prepare_submission(path=path,task=task)

    scommand = fastp+" -i "+input_fq1+" -I "+input_fq2+\
        "-o "+trimmed_fq_basename+".R1.fastq.gz "+\
        "-O "+trimmed_fq_basename+".R2.fastq.gz "+\
        "--html "+report_basename+".after_fastp.html "+\
        "--json " + report_basename + ".after_fastp.json " +\
        "-w 10;\n\n"+\
        ""\
        ""+\
        fastqc+" "+\
        trimmed_fq_basename+".R1.fastq.gz "+ \
        trimmed_fq_basename+".R2.fastq.gz "+ \
        "-o "+report_dir+";\n\n"\
        ""\
        ""+\
        bwa + " mem " +\
        ref_fasta + " " +\
        trimmed_fq_basename+".R1.fastq.gz "+\
        trimmed_fq_basename+".R2.fastq.gz "+\
        "-t "+str(threads)+" "+\
        "-R "+sam_header+" "+\
        "> "+output_aligned_sam+" "\
        "2> >(tee "+output_bam_basename+".bwa.stderr.log >&2);\n\n"+\
        ""\
        ""+\
        'echo "sambamba1";\n'+\
        sambamba + " view -S -f bam -t " + str(threads) + " " + output_aligned_sam + " > " + output_bam_basename+".bam;\n\n"+\
        'echo "sambamba2";\n'+\
        sambamba + " sort -m 5000000000 -t " + str(threads) + " " + output_bam_basename+".bam " + "-o " + output_bam_basename+".sorted.bam;\n\n"+\
        'echo "sambamba3";\n'+\
        sambamba + " index -t " + str(threads) + " " + output_bam_basename+".sorted.bam;\n\n"+\
        ""\
        'echo "remove";\n'+\
        'if test -f '+output_bam_basename+'.bam; then rm '+trimmed_fq_basename+".R*.fastq.gz "+output_fq+" "+output_aligned_sam+" fi;\n\n"\
        ""\
        'echo "grep";\n'+\
        'grep -m1 "read .* ALT contigs" '+output_bam_basename+'.bwa.stderr.log | grep -v "read 0 ALT contigs"'

    jobid = torque_submit(command=scommand, sample=sample, task=task,
                          cwd=path, to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)

def SortSam(sample, input_bam, output_bam_basename, path, output_prefix,
            jdep, gatk=gatk,  task="SortSam", cpu=1,mem="5gb",wtime="10:00:00"):

    path = prepare_submission(path=path,task=task)

    scommand = "java -jar "+gatk+\
        " SortSam "\
        "--INPUT "+input_bam+" "+\
        "--OUTPUT "+output_bam_basename+".bam "+\
        '--SORT_ORDER "coordinate" '\
        "--CREATE_INDEX true "\
        "--CREATE_MD5_FILE true "\
        "--MAX_RECORDS_IN_RAM 300000;\n\n"\
        ""\
        "" \
        'if test -f ' + output_bam_basename + '.bam; then rm ' + output_prefix+".unique.bam"+"\n"\
        'fi;\n\n'\

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)


def MarkDuplicates(sample, input_bam, output_bam_basename, metrics_filename, output_prefix,
                   path, jdep, task="MarkDuplicates", gatk=gatk,cpu=1,mem="8gb",wtime="15:00:00"):

    path = prepare_submission(path=path,task=task)

    scommand = "java -Xms6g -jar "+gatk+ \
        " MarkDuplicates " \
        "--INPUT " + input_bam + " " +\
        "--OUTPUT " + output_bam_basename+".bam "+\
        "--METRICS_FILE " + metrics_filename + " "+\
        "--VALIDATION_STRINGENCY SILENT " + \
        "--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 "\
        '--CLEAR_DT "false" '\
        '--ADD_PG_TAG_TO_READS false;\n\n'\
        ''\
        ''\
        'rm '+output_prefix+".bam "+output_prefix+".sorted.bam"

    #        '--ASSUME_SORT_ORDER "queryname" '\

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)


def BaseRecalibrator(sample, input_bam, recalibration_report_filename, path, jdep,
                     task="BaseRecalibrator", cpu=1, mem="6gb",
                     wtime="15:00:00",dbsnp_vcf=dbsnp_vcf, known_indels_sites_vcfs=known_indels_sites_vcfs,
                     ref_fasta=ref_fasta, gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "java -jar "+gatk+\
        " BaseRecalibrator "+\
        "-R "+ref_fasta+" "\
        "-I "+input_bam+" "\
        "--use-original-qualities "+\
        "-O "+recalibration_report_filename+" "\
        "--known-sites "+dbsnp_vcf+" "+\
        "--known-sites "+" -known-sites ".join([x for x in known_indels_sites_vcfs])

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)


    return(jobid)


def ApplyBQSR(sample, input_bam, output_bam_basename, recalibration_report, output_prefix,
              path, jdep, cpu=1, mem="5gb",
              wtime="10:00:00", ref_fasta=ref_fasta, task="ApplyBQSR", gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "java -XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps " \
               "-XX:+PrintGCDetails -Xloggc:gc_log.log " \
               "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -jar "+gatk+" "+\
        "ApplyBQSR "+\
        "--create-output-bam-md5 " \
        "--add-output-sam-program-record " \
        "-R "+ref_fasta+" "\
        "-I "+input_bam+" "\
        "-O "+output_bam_basename+".bam "\
        "-bqsr "+recalibration_report+" "\
        "--static-quantized-quals 10 " \
        "--static-quantized-quals 20 " \
        "--static-quantized-quals 30 " \
        "--use-original-qualities"+";\n\n"\
        ''\
        ''\
        'rm '+output_prefix+".bam "+output_prefix+".unique.sorted.*"

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path, to_wait_id=jdep, cpu=cpu,
                          mem=mem, wtime=wtime)

    return(jobid)


def CollectQualityYieldMetrics(sample, input_bam, metrics_filename, path,
                               jdep, cpu=1, mem="2gb", task="CollectQualityYieldMetrics", wtime="5:00:00", gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "java -jar "+gatk+" "+ \
               "CollectQualityYieldMetrics " \
               "--INPUT "+input_bam+" " \
               "--OUTPUT "+metrics_filename+" " \
               "--OQ true"

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path, to_wait_id=jdep, cpu=cpu,
                          mem=mem, wtime=wtime)

    return(jobid)


def AnalyzeCovariates(sample, recalibration_report_filename, plot_file, path, jdep, task="AnalyzeCovariates",
                      cpu=1, mem="4gb", wtime="5:00:00", gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "conda activate r-environment; " \
               "librjava -jar "+gatk+" "+ \
               "AnalyzeCovariates " \
               "-bqsr "+recalibration_report_filename+" " \
               "-plots "+plot_file

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path, to_wait_id=jdep, cpu=cpu,
                          mem=mem, wtime=wtime)

    return(jobid)


def CollectWgsMetrics(sample, input_bam, metrics_filename, path, jdep, cpu=1, mem="3gb", wtime="5:00:00",
                      read_length=read_length, task="CollectWgsMetrics", interval_list=interval_list, ref_fasta=ref_fasta, gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "java -Xms2000m -jar "+gatk+" "+ \
               "CollectWgsMetrics " \
               "--INPUT "+input_bam+" " \
               "--INTERVALS "+interval_list+" "\
               "--OUTPUT "+metrics_filename+" "\
               "--REFERENCE_SEQUENCE "+ref_fasta+" "\
               "--READ_LENGTH "+str(read_length)+" "\
               "--VALIDATION_STRINGENCY SILENT " \
               "--INCLUDE_BQ_HISTOGRAM true " \
               "--USE_FAST_ALGORITHM true"


    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path, to_wait_id=jdep, cpu=cpu,
                          mem=mem, wtime=wtime)

    return(jobid)

def CollectAggregationMetrics(sample, input_bam, output_bam_prefix, path, jdep, cpu=1, mem="7gb", wtime="20:00:00",
                      task="CollectAggregationMetrics", ref_fasta=ref_fasta, gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "touch "+output_bam_prefix+".gc_bias.detail_metrics "+\
      output_bam_prefix+".gc_bias.pdf "+\
      output_bam_prefix+".gc_bias.summary_metrics "+\
      output_bam_prefix+".insert_size_metrics "+\
      output_bam_prefix+".insert_size_histogram.pdf; " \
      "java -Xms5000m -jar "+gatk+" "+\
      "CollectMultipleMetrics " \
      "--INPUT "+input_bam+" " \
      "--REFERENCE_SEQUENCE "+ref_fasta+" "\
      "--OUTPUT "+output_bam_prefix+" "\
      '--ASSUME_SORTED true ' \
      '--PROGRAM "null" ' \
      '--PROGRAM "CollectAlignmentSummaryMetrics" ' \
      '--PROGRAM "CollectInsertSizeMetrics" ' \
      '--PROGRAM "CollectSequencingArtifactMetrics" ' \
      '--PROGRAM "CollectGcBiasMetrics" ' \
      '--PROGRAM "QualityScoreDistribution" ' \
      '--PROGRAM "CollectGcBiasMetrics" ' \
      '--METRIC_ACCUMULATION_LEVEL "null" ' \
      '--METRIC_ACCUMULATION_LEVEL "SAMPLE" ' \
      '--METRIC_ACCUMULATION_LEVEL "LIBRARY"' \

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)

def CollectVariantCallingMetrics(sample, input_vcf, metrics_basename, path, jdep,
                                 cpu=1, mem="4gb", wtime="20:00:00", task="CollectVariantCallingMetrics",
                                 interval_list=interval_list, dbsnp_vcf=dbsnp_vcf,
                                 ref_dict=ref_dict, gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "java -Xms2000m -jar "+gatk+" " \
               "CollectVariantCallingMetrics " \
               "--INPUT "+ input_vcf+" " \
               "--OUTPUT "+metrics_basename+" "\
               "--TARGET_INTERVALS "+interval_list+" "\
               "--SEQUENCE_DICTIONARY "+ref_dict+" "\
               "--DBSNP "+dbsnp_vcf+" "\
               "--GVCF_INPUT true"

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)


def ValidateVCF(sample, input_vcf, path, jdep, cpu=1, mem="4gb",
                wtime="20:00:00", interval_list=interval_list, task="ValidateVCF",
                dbsnp_vcf=dbsnp_vcf, ref_fasta=ref_fasta, gatk=gatk):

    path = prepare_submission(path=path,task=task)

    scommand = "java -Xms3000m -jar "+gatk+" " \
           "ValidateVariants " \
           "-V "+input_vcf+" " \
           "-R "+ref_fasta+" " \
           "-L "+interval_list+" " \
           "--validation-type-to-exclude ALLELES " \
           "--dbsnp "+dbsnp_vcf

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)


def HaplotypeCaller(sample, input_bam, base_file_name, bamout_base_name, path, jdep,
                    cpu=2, mem="20gb", wtime="100:00:00", task="HaplotypeCaller", ref_fasta=ref_fasta, ref_dict=ref_dict,
                    gatk=gatk, interval_list=interval_list):

    path = prepare_submission(path=path,task=task)

    scommand = "java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -jar " + gatk + " " \
               "HaplotypeCaller " \
               "-R "+ref_fasta+" "\
               "-I "+input_bam+" "\
               "-L "+interval_list+" "\
               "-O "+base_file_name+"g.vcf.gz "\
               "-G StandardAnnotation " \
               "-G StandardHCAnnotation " \
               "-G AS_StandardAnnotation " \
               "-new-qual " \
               "-GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 " \
               "-GQB 60 -GQB 70 -GQB 80 -GQB 90 " \
               "-ERC GVCF " \
               "-bamout "+bamout_base_name+".bamout.bam"

    jobid = torque_submit(command=scommand, sample=sample, task=task, cwd=path,
                          to_wait_id=jdep, cpu=cpu, mem=mem, wtime=wtime)

    return(jobid)





