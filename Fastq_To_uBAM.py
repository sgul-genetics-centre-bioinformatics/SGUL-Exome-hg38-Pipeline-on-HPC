#!/usr/bin/env python

#SGUL GENETICS CENTRE EXOME PIPELINE MAY 2018
#Alan Pittman

#standard pipeline steps:
###############################################################################################

# run BWA command
# run samtools sam to bam conversion command
# run samtools bam sort command 
# run samtools index
# run picard mark PCR Duplicates
# run samtools index again
# run GATK BaseRecalibrator #### need to use big VCF file
# run GATK AnalyzeCovariates ### work in progress need some R libraries installed 
# run GATK ApplyBQSR
	
###############################################################################################

import os
import sys
import subprocess
import csv
import time

#software and resources:
###############################################################################################

DRIVE = "/storage/root/homes/apittman/resources/"
java = DRIVE + "java/jre1.8.0_171/bin/java"
picard = DRIVE + "picard-2.815/picard.jar"

###############################################################################################

def display(message):
    print(message)
	
def ScriptWriter(step):
	for command in step:
		jobScript.write(str(command) + " ")
		
	jobScript.write(" \n")
	jobScript.write(" ")
	jobScript.write(" \n")
	

myProject = sys.argv
del myProject[0]

myProject = str(myProject)
myProject = myProject.lstrip('[')
myProject = myProject.lstrip("'")
myProject = myProject.rstrip(']')
myProject = myProject.rstrip("'")    

print("###################################################################################\n"
      "###########  Welcome to the PexOmes Pipelines by SGUL Genetics Centre  ############\n"
      "#                                                                                 #\n"
      "#  This scripy takes your raw sequencing data in fastq format and converts it to  #\n"
	  "#                                  uBAM file                                      #\n"
      "#---------------------------------------------------------------------------------#\n"
      "#Developers: Dionysios Grigoriadis, Dr Alan Pittman                               #\n"
      "#                                                                                 #\n"
      "#####     St George's University of London Genetics Centre, October 2019      #####\n"
      "###################################################################################\n")

print("your project is :")
print(myProject)

print("\n")

dirpath = os.getcwd()
projectDataDir=dirpath + "/" + "Unaligned" + "/" + myProject + "/"
print(os.listdir(projectDataDir))

#generate list of all the samples in the project directory 
mySamples = os.listdir(projectDataDir)

#for each sample in the project folder, identify the fastq files (R1 and R2) for BWA alignment
for sample in mySamples:

	sampleDirectory = projectDataDir + sample + "/"

	myinputfiles = os.listdir(sampleDirectory) # assuming we have two fastq files per sample
	
	inputFASTQ1 = myinputfiles[0]
	inputFASTQ2 = myinputfiles[1]
	
	print("your input fastq files are:")
	
	#L1_2.fq.gz # check naming format of fastq files 
	print(inputFASTQ1)
	print(inputFASTQ2)
		
	#make analysis output directory for alignment
	
	uBAMdir = "/storage/root/homes/apittman/exomes_hg38/UBAMs/"
	
	projectFolder = uBAMdir + myProject 
	makeDirectoryProject = ['mkdir', projectFolder] 
	subprocess.call(makeDirectoryProject)
	
	sampleFolder = uBAMdir + myProject + "/" + sample
	makeDirectory = ['mkdir', sampleFolder] 
	subprocess.call(makeDirectory)
		
	############# SAMPLE SPECIFIC VARIABLES ###################################

	path_inputFASTQ1 = sampleDirectory + inputFASTQ1
	path_inputFASTQ2 = sampleDirectory + inputFASTQ2
	
	F1 = "F1=" + path_inputFASTQ1 
	F2 = "F2=" + path_inputFASTQ2
	
	SAMPLE_NAME = "SM=" + sample
	
	outputBAM = sampleFolder + "/" + sample + ".bam" # output bam file
	
	picardOuBAM = "O=" + outputBAM
	
	########## Job submission at top of script ########################

	interpreter = "#!/bin/bash"
	nodesandcores = "#PBS -l nodes=1:ppn=10"
	walltime = "#PBS -l walltime=36:00:00"
	pbs = "#PBS -S /bin/bash"
	abe = "#PBS -m abe"
	email = "#PBS -M apittman@sgul.ac.uk"
		
	FastQtouBAM = [java, 
	'-jar', 
    picard, 
    'FastqToSam',
    F1,
    F2,
    SAMPLE_NAME,	
	picardOuBAM]
	
#	java -jar picard.jar FastqToSam \
#      F1=file_1.fastq \
#      O=fastq_to_bam.bam \
#      SM=for_tool_testing 
		
	################ Now writing the job scripts ######################
	
	ScriptDirectory = dirpath + "/" + "Jobs" + "/" + myProject + "/"
		
	if not os.path.exists(ScriptDirectory):
		os.mkdir(ScriptDirectory)
		print("Directory " , ScriptDirectory ,  " Created ")
	else:    
		print("Directory " , ScriptDirectory ,  " already exists")
	
	jobScript = open( ScriptDirectory + sample + ".sh","w") # add path here
	
	
	jobScript.write(str(interpreter))
	jobScript.write(" \n")
	jobScript.write(str(nodesandcores))
	jobScript.write(" \n")
	jobScript.write(str(walltime))
	jobScript.write(" \n")
	jobScript.write(str(pbs))
	jobScript.write(" \n")
	jobScript.write(str(abe))
	jobScript.write(" \n")
	jobScript.write(str(email))
	jobScript.write(" \n")
	jobScript.write(" ")
	jobScript.write(" \n")
	
	ScriptWriter(FastQtouBAM)

	jobScript.write(" \n")
	jobScript.write("exit")
	jobScript.write(" \n")
		
	jobScript.close()	

	################ Submit Jobs to Q #####################################
	
	job = ScriptDirectory + sample + ".sh"
	
	submission = ['qsub', job]
	
	subprocess.call(submission)
	
	time.sleep(3)

