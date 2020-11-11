import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil

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
    print(red('Could not import h5py module. This is required for analyse functions.'))
    exit(0)

# Start the timer
start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for merging the PartState container of multiple .h5 files into a single file.\nSupply a single state file or a group of state files by using the wildcard "*", e.g. MySimulation_000.000000* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Files (.h5) that are to be merged together.', nargs='+')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
print('='*132)
print("Running with the following command line options")
for arg in list(args.__dict__) :
    print(arg.ljust(15)+" = [ "+str(getattr(args,arg))+" ]")
print('='*132)

# Loop over all files and identify the latet simulation time
maxtime=-1.0
NbrOfFiles = len(args.files)
files = []
for statefile in args.files :
    #print(statefile)
    pat = r'^.*\_(.*)\..*$'
    match = re.search(pat, statefile)
    timestr = match.group(1)
    try :
        time = float(timestr)
        maxtime=max(time,maxtime)
        if time == maxtime :
            maxtimestr = timestr
            maxtimeFile = statefile
            newFile = re.sub(timestr+'.h5', '', statefile)+maxtimestr+'_merged.h5'
        files.append(statefile)
    except :
        print("not considering "+statefile)

print("t_max     : %s" % maxtime)
print("t_max_str : %s" % timestr)
print("newfile   : %s" % newFile)
print()

# Get maximum number of characters in h5 file names
max_length=0
for statefile in files :
    max_length = max(max_length,len(statefile))

print(132*"-")
s="Example.h5"
print("".ljust(max_length-len(s)),s," | PartData(dim1, dim2)")
print(132*"-")
n = 0
for statefile in files :
    n+=1
    # Open h5 file and read container info
    # --------------------------------------------
    #     r       : Readonly, file must exist
    #     r+      : Read/write, file must exist
    #     w       : Create file, truncate if exists
    #     w- or x : Create file, fail if exists
    #     a       : Read/write if exists, create otherwise (default
    # --------------------------------------------
    # When sorting is used, the sorted array is written to the original .h5 file with a new name
    f1 = h5py.File(statefile,'r+')

    file_version  = f1.attrs.get('File_Version', default=-1.)[0]
    
    # Usage:
    # -------------------
    # available keys         : print("Keys: %s" % f1.keys())                                # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
    # first key in list      : a_group_key = list(f1.keys())[0]                             # yields 'DG_Solution'
    # available attributes   : print('\n'.join(x for x in f1.attrs))                        # yields 'File_Type\n File_Version\n MeshFile'
    # get specific attribute : file_version  = f1.attrs.get('File_Version', default=-1.)[0] # yields, e.g., 1.5
    # -------------------
    data_set= 'PartData'
    
    # 1.1.1   Read the dataset from the hdf5 file
    b1 = f1[data_set][:]
    print("".ljust(max_length-len(statefile)),statefile," | PartData%s" % str(b1.shape))

    # Save old file
    if n > 1 :
        # Compare shape of the dataset of both files, throw error if they do not conincide
        #if b1.shape[0] != b2.shape[0] : # e.g.: b1.shape = (48, 1, 1, 32)
            #s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\nThe dimensions dim1 = %s and dim2 = %s must be equal!\n\nAborted!" % (statefile,statefile_old,b1.shape,b2.shape,b1.shape[0],b2.shape[0])
            #print(s)
            #exit(1)
        # Concatenate files depending on the file version (beginning with v1.5 PartData's dimensions are switched)
        if file_version < 1.5 :
            # Concatenate columns
            b1_merged = np.concatenate((b1_merged, b1), axis=1)
        else :
            # Concatenate rows
            b1_merged = np.concatenate((b1_merged, b1), axis=0)
    else :
        b1_merged = b1
    b2 = b1
    statefile_old = statefile
    f1.close()
print(132*"-")
print("Files have been merged into %s | %s%s" % (newFile,data_set,b1_merged.shape))



# Copy old file and modify PartState in the new file
shutil.copyfile(maxtimeFile, newFile)
# Write new .h5 file
f1 = h5py.File(newFile,'r+')
del f1[data_set]
# file: Create new dataset
dset = f1.create_dataset(data_set, shape=b1_merged.shape, dtype=np.float64)
# write as C-continuous array via np.ascontiguousarray()
dset.write_direct(np.ascontiguousarray(b1_merged))
f1.close()
