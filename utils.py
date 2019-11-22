import os
import sys
import time
import glob
import datetime
import subprocess
import json


def torque_qcomm(command, sample, task, to_wait_id="",
                  wtime="24:00:00", nodes=1, cpu=1, mem="1gb", cwd="./"):
    """Function which takes a pipeline command, runtime arguments, jobname, job dependencies
    as input and returns the job submission script"""
    out_filename = cwd + "/" + task + ".jobscript"
    with open(out_filename, 'w+') as out_file:
        out_file.write('#!/bin/bash')
        out_file.write('\n')
        out_file.write('#PBS -V')
        out_file.write('\n')
        out_file.write('#PBS -N '+str(sample)+"."+str(task))
        out_file.write('\n')
        out_file.write('#PBS -V')
        out_file.write('\n')
        if to_wait_id!="":
            out_file.write('#PBS -W '+"depend=afterok:"+to_wait_id)
            out_file.write('\n')
        out_file.write("#PBS -lwalltime="+str(wtime)+",nodes="+str(nodes)+":ppn="+str(cpu)+",mem="+str(mem))
        out_file.write('\n')
        out_file.write("#PBS -d "+str(cwd))
        out_file.write('\n')
        out_file.write("#PBS -e " + str(cwd + "/" + task + ".stderr"))
        out_file.write('\n')
        out_file.write("#PBS -o " + str(cwd + "/" + task + ".stdout"))
        out_file.write('\n')
        out_file.write('\n')
        out_file.write(command)
        out_file.write('\n')
        out_file.write('sleep 5')
        out_file.write('\n')
        out_file.write("exit")
    out_file.close()

    return(out_file)


def torque_submit(command, sample, task, to_wait_id="",
                 wtime="24:00:00", nodes=1, cpu=1, mem="1gb", cwd="./"):
    """Function which takes a pipeline command, runtime arguments, jobname, job dependencies
    as input, submits the job to the cluster and returns the job id"""
    out_filename = cwd + "/" + task + ".jobscript"
    with open(out_filename, 'w+') as out_file:
        out_file.write('#!/bin/bash')
        out_file.write('\n')
        out_file.write('#PBS -V')
        out_file.write('\n')
        out_file.write('#PBS -N ' + str(sample) + "." + str(task))
        out_file.write('\n')
        if to_wait_id != "":
            out_file.write('#PBS -W ' + "depend=afterok:" + to_wait_id)
            out_file.write('\n')
        out_file.write(
            "#PBS -lwalltime=" + str(wtime) + ",nodes=" + str(nodes) + ":ppn=" + str(cpu) + ",mem=" + str(mem))
        out_file.write('\n')
        out_file.write("#PBS -d " + str(cwd))
        out_file.write('\n')
        out_file.write("#PBS -e " + str(cwd + "/" + task + ".stderr"))
        out_file.write('\n')
        out_file.write("#PBS -o " + str(cwd + "/" + task + ".stdout"))
        out_file.write('\n')
        out_file.write('\n')
        out_file.write(command)
        out_file.write('\n')
        out_file.write('sleep 5')
        out_file.write('\n')
        out_file.write("exit")
    out_file.close()

    p = subprocess.Popen(["qsub",out_filename], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    jout = out.decode("utf-8").strip("\n")
    jerr = err.decode("utf-8").strip("\n")

    if jerr != "":
        print(jerr)
        raise()

    return (jout)

def prepare_submission(path,task):

    if (path.endswith("/")):
        pass
    else:
        path = path + "/"

    if not os.path.exists(path):
        os.makedirs(path)

    if not os.path.exists(path+task):
        os.makedirs(path+task)

    return(path+task)

def hasRQJob(job):
    p = subprocess.Popen(["qstat", job], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    status = out.decode("utf-8").strip("\n").split()[-2] in ('Q','R')
    return(status)

def check_last_job(job):
    p = subprocess.Popen(["qstat", job], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    status = out.decode("utf-8").strip("\n").split()
    if job in status and status[-2] in ('H','R','Q'):
        return(job)
    else:
        return("")