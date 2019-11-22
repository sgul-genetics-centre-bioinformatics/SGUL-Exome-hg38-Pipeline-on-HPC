# SGUL-Exome-hg38-Pipeline-on-HPC

This is SGUL Genetic Centre's pipeline for pre-processing and variant calling of raw Exome Sequencing data on TORQUE/PBS HPC scheduling system.

This is a python pipeline equivalent of the GATK Best practices for Exome sequencing found here:https://github.com/gatk-workflows/gatk4-exome-analysis-pipeline/tree/master/tasks#

Joint Genotyping step is not available at the moment. Under construction :)

Before running the pipeline ensure you have the suitable directories tree by using this script:
dirs_create.sh


* PexOmes_raw2gvcf.py
Main executable python module
run: python PexOmes_raw2gvcf.py -p dios_project

help: python PexOmes_raw2gvcf.py --help

* dependencies_and_settings.py
Resources and software dependencies

* tasks.py
Steps of the pipeline. Each function create and runs a jobscipt and returns the id of the job

* utils.py	
Utilities - Usefull functions for the pipeline
