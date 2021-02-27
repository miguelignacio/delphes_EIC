#!/usr/bin/env python
#SBATCH -J EICSTUDY    # job name
#SBATCH -o logs/simple-%A_%a.out
#SBATCH -e logs/simple-%A_%a.out
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
# SLURM_ARRAY_TASK_ID=0 ./simpleanalysis.py --input <DIRECTORY CONTAINING FILES> --name <OUTPUT DIRECTORY>
#
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


parser = argparse.ArgumentParser()

parser.add_argument("-n", "--name", type=str,
                    help="name for this study set (e.g. bfield)")
parser.add_argument("-i", "--input", type=str,
                    help="input directory of Delphes ROOT files")
parser.add_argument("-m", "--modules", type=str,
                    help="analysis module list")
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



# Load all the ROOT files to Process

root_files = glob.glob(args.input + "/*/*.root")
root_files = sorted(root_files)

if len(root_files) == 0:
    print("No appropriate ROOT files exist in the input directory")
    sys.exit()
else:
    print("Input directory contains %d files for processing" % (len(root_files)))


#print(root_files)


# Use the task ID to load the necessary file from the array
fileNumber = int(SLURM_ARRAY_TASK_ID)

# Search out the file number in this list
root_file = None
for test_file in root_files:
    if test_file.find('%d/out.root' % fileNumber) != -1:
        root_file = test_file
        break

if root_file == None:
    print("Unable to find an input file consistent with the task ID = %d" %
          (int(SLURM_ARRAY_TASK_ID)))
    sys.exit()

# Get the full path to the input file
root_file = os.path.abspath(root_file)
print("Processing %s" % (root_file))


# Execute the study

taskdir="%s/%d" % (args.name, fileNumber)

if os.path.exists(taskdir) and not args.force:
    print("Skipping this task directory --- it already exists. Cleanup before overwriting!")
    print(taskdir)
else:
    if not os.path.exists(taskdir):
        os.makedirs(taskdir)

    # Copy or Link needed files to working directory
    subprocess.call("cp -a mva_taggers %s/" % (taskdir), shell=True);
    # Execute the study
    subprocess.call('cd {0[taskdir]}; SimpleAnalysis.exe --input_dir {0[root_file]} --output_file out.root --module_sequence "{0[modules]}"'.format({'taskdir': taskdir, 'root_file': root_file, 'modules': args.modules}), shell=True)
