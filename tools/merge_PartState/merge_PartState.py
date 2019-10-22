import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil

# import h5 I/O routines
try :
    import h5py
    h5py_module_loaded = True
except ImportError :
    #raise ImportError('Could not import h5py module. This is required for anaylze functions.')
    print(tools.red('Could not import h5py module. This is required for anaylze functions.'))
    exit(0)


start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for merging the PartState container of multiple .h5 files into a single file.\nSupply a single state file or a group of state files by using the wildcard "*", e.g. MySimulation_000.000000* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
#parser.add_argument('-c', '--carryon', action='store_true', help='''Continue build/run process. 
  #--carryon         : build non-existing binary-combinations and run all examples for thoses builds
  #--carryon --run   : run all failed examples''')
#parser.add_argument('-e', '--exe', help='Path to executable of code that should be tested.')
#parser.add_argument('-d', '--debug', type=int, default=0, help='Debug level.')
#parser.add_argument('-j', '--buildprocs', type=int, default=0, help='Number of processors used for compiling (make -j XXX).')
#parser.add_argument('-b', '--basedir', help='Path to basedir of code that should be tested (contains CMakeLists.txt).')
#parser.add_argument('-y', '--dummy', action='store_true',help='Use dummy_basedir and dummy_checks for fast testing on dummy code.')
#parser.add_argument('-r', '--run', action='store_true' ,help='Run all binaries for all examples with all run-combinations for all existing binaries.')
#parser.add_argument('-s', '--save', action='store_true',help='Do not remove output directories buildsXXXX in output_dir after successful run.')
#parser.add_argument('-t', '--compiletype', help='Override all CMAKE_BUILD_TYPE settings by ignoring the value set in builds.ini (e.g. DEBUG or RELEASE).')
#parser.add_argument('-a', '--hlrs', action='store_true', help='Run on with aprun (hlrs system).')
#parser.add_argument('-z', '--rc', dest='referencescopy', help='Create/Replace reference files that are required for analysis. After running the program, the output files are stored in the check-/example-directory.', action='store_true')
#parser.set_defaults(referencescopy=False)
parser.add_argument('files', type=str, help='Files (.h5) that are to be merged together.', nargs='+')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
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
    
    # Usage:
    # -------------------
    # available keys   : print("Keys: %s" % f1.keys())         # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
    # first key in list: a_group_key = list(f1.keys())[0]      # yields 'DG_Solution'
    # -------------------
    data_set= 'PartData'
    
    # 1.1.1   Read the dataset from the hdf5 file
    b1 = f1[data_set][:]
    print(statefile,b1.shape)

    # Save old file
    if n > 1 :
        # Compare shape of the dataset of both files, throw error if they do not conincide
        if b1.shape[0] != b2.shape[0] : # e.g.: b1.shape = (48, 1, 1, 32)
            s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\nThe dimensions dim1 = %s and dim2 = %s must be equal!\n\nAborted!" % (statefile,statefile_old,b1.shape,b2.shape,b1.shape[0],b2.shape[0])
            print(s)
            exit(1)
        b1_merged = np.concatenate((b1_merged, b1), axis=1)
    else :
        b1_merged = b1
    b2 = b1
    statefile_old = statefile
    f1.close()



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
