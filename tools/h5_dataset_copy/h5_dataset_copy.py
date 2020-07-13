import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil
import os

class bcolors :
    """color and font style definitions for changing output appearance"""
    # Reset (user after applying a color to return to normal coloring)
    ENDC   ='\033[0m'    

    # Regular Colors
    BLACK  ='\033[0;30m' 
    RED    ='\033[0;31m' 
    GREEN  ='\033[0;32m' 
    YELLOW ='\033[0;33m' 
    BLUE   ='\033[0;34m' 
    PURPLE ='\033[0;35m' 
    CYAN   ='\033[0;36m' 
    WHITE  ='\033[0;37m' 

    # Text Style
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def red(text) :
    return bcolors.RED+text+bcolors.ENDC

def green(text) :
    return bcolors.GREEN+text+bcolors.ENDC

def blue(text) :
    return bcolors.BLUE+text+bcolors.ENDC

def yellow(text) :
    return bcolors.YELLOW+text+bcolors.ENDC

# import h5 I/O routines
try :
    import h5py
    h5py_module_loaded = True
except ImportError :
    print(red('Could not import h5py module. This is required for handling .h5 files.'))
    exit(0)

# Start the timer
start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for copying a dataset (e.g. PartData) from a .h5 state file to another file. Supply original datafile, target data file and dataset name.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-s', '--source', type=str, help='Source .h5 file containing the dataset that is to be copied.')
parser.add_argument('-t', '--target', type=str, help='Target .h5 file into which the dataset is to be copied.')
parser.add_argument('-d', '--dataset', type=str, help='Name of the h5 dataset in the source .h5 file.')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
print('='*132)
print("Running with the following command line options")
for arg in list(args.__dict__) :
    print(arg.ljust(15)+" = [ "+str(getattr(args,arg))+" ]")
print('='*132)
print()

# Check arguments
if not (args.source and os.path.exists(args.source)):
    parser.error('No source file given or file does not exist, add -s or --source followed by the path to the source file')
if not (args.target and os.path.exists(args.target)):
    parser.error('No target file given or file does not exist, add -t or --target followed by the path to the source file')
if not (args.dataset):
    parser.error('No dataset given, add -d or --dataset followed by the name of the dataset in the source file')



# Open h5 file and read container info
# --------------------------------------------
#     r       : Readonly, file must exist
#     r+      : Read/write, file must exist
#     w       : Create file, truncate if exists
#     w- or x : Create file, fail if exists
#     a       : Read/write if exists, create otherwise (default
# --------------------------------------------
# 1. Open .h5 source file
try:
    f1 = h5py.File(args.source,'r')
except Exception as e:
    print(e)
    print("Could not open the source file [%s]" % (args.source))
    exit(0)

# Usage:
# -------------------
# available keys         : print("Keys: %s" % f1.keys())                                # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
# first key in list      : a_group_key = list(f1.keys())[0]                             # yields 'DG_Solution'
# available attributes   : print('\n'.join(x for x in f1.attrs))                        # yields 'File_Type\n File_Version\n MeshFile'
# get specific attribute : file_version  = f1.attrs.get('File_Version', default=-1.)[0] # yields, e.g., 1.5
# -------------------

# 2. Read the dataset from the hdf5 file
try :
    dset1 = f1[args.dataset][:]
    dtype1 = f1[args.dataset].dtype
except :
    print(e)
    print("Could not read dataset [%] from the source file [%s]" % (args.dataset,args.source))
    exit(0)

print('Source file'.ljust(len(args.source))," | dataset(shape) data type")
print(132*'-')
#print("".ljust(-len(args.source)),args.source," | %s%s %s" % (args.dataset, str(dset1.shape),dtype1))
print(args.source.ljust(len(args.source))," | %s%s %s" % (args.dataset, str(dset1.shape),dtype1))




# 3. Open .h5 target file
try:
    f2 = h5py.File(args.target,'r+')
except Exception as e:
    print(e)
    print("Could not open the target file [%s]" % (args.target))
    exit(0)

# 0. remove dataset container in target file (if it exists)
try:
    del f2[args.dataset]
except Exception as e:
    #print(e)
    #print("Dataset does not exist in target file (%s). Copying." % e)
    pass


f2.create_dataset(args.dataset, data=dset1)

f1.close()
f2.close()
print("Done.")
