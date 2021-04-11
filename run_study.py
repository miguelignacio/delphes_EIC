#!/usr/bin/env python
#SBATCH -J EICSTUDY    # job name
#SBATCH -o logs/eicstudy-%A_%a.out
#SBATCH -e logs/eicstudy-%A_%a.out
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -p htc             # queue (partition) -- batch, parallel, etc.  parallel-medium
#SBATCH -t 12:00:00        # run time (hh:mm:ss)
#SBATCH -D .               # Directory where executable will be run
#SBATCH --mail-user=ssekula@smu.edu
#SBATCH --mail-type=fail    # email me when the job finishes


# This script was written originally for use on a SLURM-based batch system,
# like the one available at SMU ("ManeFrame II"). It can be run on the command-line
# alone; if the script doesn't detect the requested SLURM environment variables,
# it will ask you to specify them. For instance, to run the first variation in a
# study,
#
# SLURM_ARRAY_TASK_ID=0 ./run_study.py -t ./bfield_study_template.tcl -p "['PARAM_BFIELD']" -v "[[1.0,1.5,2.0,2.5,3.0,3.5]]" -n bfield
#
# The above will manually select the first variation in the array passed to "-v" and run it.
#
# General idea:
#
# A template script is executed by the Delphes simulation executable. Before execution,
# placeholder parameters are filled in by the script and a copy of the template is 
# generated (with the parameter values set) for actual execution. For example:
#
# SLURM_ARRAY_TASK_ID=0 ./run_study.py -t ./bfield_study_template.tcl -p "['PARAM_BFIELD']" -v "[[2.0]]" -n bfield
# 
# will open ./bfield_study_template.tcl, find all instanced of PARAM_BFIELD, and replace them with the number "2.0".
# This allows command-line variation of one or more parameters to perform a study of changes to simulation.
#

import subprocess
import math
import os
import sys
import shutil
import glob
import re
import ast

import argparse

global tclparams, pythiaparams
tclparams = {}
pythiaparams = {}


def LoadParamsTcl(tclfile):
    tcl = open(tclfile)

    param_pattern = re.compile("^set\s+(.*?)\s+(.*)")
    for line in tcl:
        line = line.rstrip('\n')
        match = param_pattern.match(line)
        if match != None:
            tclparams[match.group(1)] = str(match.group(2))
    
    print("Customization TCL Parameters are:")
    print(tclparams)

def LoadParamsPythia(cmndfile):
    cmnd = open(cmndfile)

    param_pattern = re.compile("^\!\s*(.*?)\s*=\s*(.*)")
    for line in cmnd:
        line = line.rstrip('\n')

        if line.find('END HEADER') != -1:
            break

        match = param_pattern.match(line)
        if match != None:
            pythiaparams[match.group(1)] = str(match.group(2))
    
    print("Customization Pythia8 Parameters are:")
    print(pythiaparams)


def WriteTclFile(filename):
    global args
    tclfile = open(filename, 'w')

    for param in args.params:
        value = args.params[param]
        if (param not in tclparams.keys()) and (param not in pythiaparams.keys()):
            print("WARNING: you tried to set %s, which is not an available parameter!" % (param))
            continue
        tclparams[param] = str(value)

    for param in tclparams:
        line = "set %s %s\n" % (param, tclparams[param])
        tclfile.write(line)
    tclfile.close()
    
def WritePythiaFile(output_filename):
    global args

    for param in args.params:
        value = args.params[param]
        if (param not in pythiaparams.keys()) and (param not in tclparams.keys()):
            print("WARNING: you tried to set %s, which is not an available parameter!" % (param))
            continue
        pythiaparams[param] = str(value)


    with open(args.commands, "rt") as input_template:
        with open(output_filename, "wt") as output_template:
            for line in input_template:
                for param in pythiaparams:
                    value = str(pythiaparams[param])
                    line = line.replace(param, value)
                    
                output_template.write(line)


def TransferTclFile(tclfile, taskdir):
    with open(tclfile, "rt") as input_template:
        with open("%s/%s" % (taskdir, tclfile), "wt") as output_template:
            for line in input_template:
                line = line.replace("customizations.tcl", "%s/customizations.tcl" % (taskdir))
                    
                output_template.write(line)



parser = argparse.ArgumentParser()

parser.add_argument("-n", "--name", type=str,
                    help="name for this study set (e.g. bfield)")
parser.add_argument("-t", "--template", type=str,
                    help="template TCL file for study")
parser.add_argument("-c", "--commands", type=str,
                    help="template command file for study")
parser.add_argument("-p", "--params", type=ast.literal_eval, default={},
                    help="environment variables to set for study")
#parser.add_argument("-v", "--values", type=ast.literal_eval, default=[],
#                    help="parameter values for study")
parser.add_argument("-f", "--force", default=False, action='store_true',
                    help="force-overwrite existing output")

global args
args = parser.parse_args()

# Create the task superdirectory

if not os.path.exists(args.name):
    try:
        os.makedirs(args.name)
    except OSError:
        print("%s already exists... continuing..." % (args.name))


SLURM_ARRAY_TASK_ID="0"

try:
    SLURM_ARRAY_TASK_ID=os.environ["SLURM_ARRAY_TASK_ID"]
except:
    print("Please set the SLURM_ARRAY_TASK_ID environment variable to a number (e.g. 0) before running this script.")
    sys.exit()
    pass

print("Task ID requested: %d" % (int(SLURM_ARRAY_TASK_ID)))



value_index = int(SLURM_ARRAY_TASK_ID)


# Load the available parameters and their default values from
# customizations.tcl and from the header of the PYTHIA8 command
# file

LoadParamsTcl("customizations.tcl")
LoadParamsPythia(args.commands)


# Handle random number seed issues and other defaults

random_seed = 0
if "PARAM_RANDOM_SEED" not in args.params:
    random_seed = abs(hash(args.name)) % (10 ** 8) + value_index
    pythiaparams["PARAM_RANDOM_SEED"] = random_seed

# Execute the study

taskdir="%s/%d" % (args.name, value_index)
tclfile = "%s/customizations.tcl" % (taskdir)
cmndfile = "%s/%s" % (taskdir, args.commands)

if os.path.exists(taskdir) and not args.force:
    print("Skipping this task directory --- it already exists. Cleanup before overwriting!")
    print(taskdir)
else:
    if not os.path.exists(taskdir):
        os.makedirs(taskdir)

    # Replace parameter placeholders with values for this study
    WriteTclFile(tclfile)
    WritePythiaFile(cmndfile)

    # Copy files to the task directory before execution
    copy_files = [args.template]
    for a_file in copy_files:
        subprocess.call("cp %s %s" % (a_file, taskdir), shell=True)
        #TransferTclFile(a_file, taskdir)
    # Write the random number seed to disk
    rndm_file = open(taskdir+"/random_seed.dat", "w")
    rndm_file.write(str(random_seed))
    rndm_file.close()

    # Execute the study
    subprocess.call("DelphesPythia8 {0[taskdir]}/{0[template]} {0[taskdir]}/{0[commands]} {0[taskdir]}/out.root".format({'taskdir': taskdir, 'template': args.template, 'commands': args.commands}), shell=True)
