#!/usr/bin/env python
#SBATCH -J EICSKIM    # job name
#SBATCH -o logs/eicskim-%A_%a.out
#SBATCH -e logs/eicskim-%A_%a.out
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=6G
#SBATCH -p htc             # queue (partition) -- batch, parallel, etc.  parallel-medium
#SBATCH -t 00:05:00        # run time (hh:mm:ss)
#SBATCH -D .               # Directory where executable will be run
#SBATCH --mail-user=ssekula@smu.edu
#SBATCH --mail-type=fail    # email me when the job finishes


# This script was written originally for use on a SLURM-based batch system,
# like the one available at SMU ("ManeFrame II"). It can be run on the command-line
# alone; if the script doesn't detect the requested SLURM environment variables,
# it will ask you to specify them. For instance, to run the first variation in a
# study,
#
# SLURM_ARRAY_TASK_ID=0 ./run_skim.py -i <INPUT DIRECTORY> -o <OUTPUT DIRECTORY> -f out.root -c "Jet.Flavor==4"
#
# will skim the <INPUT DIRECTORY>/0/out.root file into <OUTPUT DIRECTORY>/0/out.root file, using the cut
# Jet.Flavor==4. Any event containing such a jet will be retained in the skim.

import subprocess
import math
import os
import sys
import shutil
import glob
import re
import ast

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", type=str,
                    help="directory holding all input ROOT files for skimming")
parser.add_argument("-o", "--output", type=str,
                    help="directory holding all input ROOT files for skimming")
parser.add_argument("-r", "--rootfile", type=str,
                    help="Name of the ROOT file in each subdirectory of the input directory")
parser.add_argument("-c", "--cuts", type=str,
                    help="ROOT selection string-style cuts")
parser.add_argument("-f", "--force", default=False, action='store_true',
                    help="force-overwrite existing output")

args = parser.parse_args()



# Create the task superdirectory

if not os.path.exists(args.output):
    try:
        os.makedirs(args.output)
    except OSError:
        print("%s already exists... continuing..." % (args.output))


SLURM_ARRAY_TASK_ID="0"

try:
    SLURM_ARRAY_TASK_ID=os.environ["SLURM_ARRAY_TASK_ID"]
except:
    print("Please set the SLURM_ARRAY_TASK_ID environment variable to a number (e.g. 0) before running this script.")
    sys.exit()
    pass

print("Task ID requested: %d" % (int(SLURM_ARRAY_TASK_ID)))



value_index = int(SLURM_ARRAY_TASK_ID)


# Execute the skim

taskdir="%s/%d" % (args.output, value_index)
inputfile="%s/%d/%s" % (args.input, value_index, args.rootfile)
outputfile="%s/%d/%s" % (args.output, value_index, args.rootfile)

if os.path.exists(taskdir) and not args.force:
    print("Skipping this task directory --- it already exists. Cleanup before overwriting!")
    print(taskdir)
else:
    if not os.path.exists(taskdir):
        os.makedirs(taskdir)

    # Execute the study
    subprocess.call("root -q -l -b ./DelphesSkim.C'+(\"{0[input]}\",\"{0[output]}\",\"{0[cuts]}\")'".format({'input': inputfile, 'output': outputfile, 'cuts': args.cuts}), shell=True)
