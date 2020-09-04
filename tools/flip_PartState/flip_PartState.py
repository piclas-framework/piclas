import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil

def ChangeFileVersion(statefile) :
    # Open h5 file and read container info
    # --------------------------------------------
    #     r       : Readonly, file must exist
    #     r+      : Read/write, file must exist
    #     w       : Create file, truncate if exists
    #     w- or x : Create file, fail if exists
    #     a       : Read/write if exists, create otherwise (default
    # --------------------------------------------
    # When sorting is used, the sorted array is written to the original .h5 file with a new name
    # 1. Open .h5 data file
    f1 = h5py.File(statefile,'r+')

    # 2. Get curret file version
    file_version = f1.attrs.get('File_Version', default=-1.)
    if type(file_version) is type([]) :
        file_version = file_version[0]

    # 3. Change the file version to 1.4 or 1.5 depending on the original file version
    if file_version > 0. :
        if file_version < 1.5 :
            file_version_new=1.5
        else :
            file_version_new=1.4
        f1.attrs.modify('File_Version',file_version_new)
        print( )
        print("Changed file version from %s to %s" % (file_version,file_version_new))

    # 4. Close .h5 data file
    f1.close()

def FlipDataset(statefile,data_set) :
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

    # Usage:
    # -------------------
    # available keys         : print("Keys: %s" % f1.keys())                                # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
    # first key in list      : a_group_key = list(f1.keys())[0]                             # yields 'DG_Solution'
    # available attributes   : print('\n'.join(x for x in f1.attrs))                        # yields 'File_Type\n File_Version\n MeshFile'
    # get specific attribute : file_version  = f1.attrs.get('File_Version', default=-1.)[0] # yields, e.g., 1.5
    # -------------------

    # 1.1.1   Read the dataset from the hdf5 file
    try :
        b1 = f1[data_set][:]
    except :
        #print('Dataset %s does not exist' % data_set)
        return

    dataType = f1[data_set].dtype

    if args.debug:
        print( )
        print(yellow("    size :         "+ str(f1[data_set].size)))
        print(yellow("    shape :        "+ str(f1[data_set].shape)))
        print(yellow("    dtype :        "+ str(f1[data_set].dtype)))
        print(yellow("    chunks :       "+ str(f1[data_set].chunks)))
        print(yellow("    compression :  "+ str(f1[data_set].compression)))
        print(yellow("    shuffle :      "+ str(f1[data_set].shuffle)))

    old_shape=b1.shape

    # Switch dimenions and store in file
    # -----------------------------------------
    # 0. remove original PartData container
    del f1[data_set]

    # 1. Swtich dimension of PartData
    b1 = np.swapaxes(b1,0,1)
    new_shape=b1.shape
    print("".ljust(max_length-len(statefile)),statefile," | %s%s => %s%s" % (data_set, str(old_shape), data_set, str(new_shape)))

    # 2. Create new dataset 'dset'
    #dset = f1.create_dataset(data_set, shape=b1.shape, dtype=np.float64)
    dset = f1.create_dataset(data_set, shape=b1.shape, dtype=getattr(np, str(dataType)))

    # 3. Write as C-continuous array via np.ascontiguousarray()
    if not b1.any() :
        print(" %s has dimension %s. Skipping" % (data_set,b1.shape))
        pass
    else :
        dset.write_direct(np.ascontiguousarray(b1))

    # 4. Close .h5 data file
    f1.close()

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
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for switching the dimensions of the PartState container of multiple .h5 files.\nSupply a single state file or a group of state files by using the wildcard "*", e.g. MySimulation_000.000000* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Files (.h5) that are to be merged together.', nargs='+')
parser.add_argument('-d', '--debug', action='store_true', help='Print additional imformation regarding the files onto screen.')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
print('='*132)
print("Running with the following command line options")
for arg in list(args.__dict__) :
    print(arg.ljust(15)+" = [ "+str(getattr(args,arg))+" ]")
print('='*132)
print()

# Get maximum number of characters in h5 file names
max_length=0
for statefile in args.files :
    max_length = max(max_length,len(statefile))

print(132*"-")
s="Example.h5"
print("".ljust(max_length-len(s)),s," | PartData(dim1, dim2)")
print(132*"-")
for statefile in args.files :
    # Flip PartData
    FlipDataset(statefile,'PartData')

    # Check for additional arrays: 'CloneVibQuantData', 'VibQuantData' and 'CloneData'
    FlipDataset(statefile,'CloneVibQuantData')
    FlipDataset(statefile,'VibQuantData')
    FlipDataset(statefile,'CloneData')

    # Change file version
    ChangeFileVersion(statefile)

print(132*"-")

