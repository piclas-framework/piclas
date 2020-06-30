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

def GetDataSets(statefile) :
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
    data_set = 'SurfaceData2'

    datasets = list(f1.keys())
    
    f1.close()

    return datasets
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
parser = argparse.ArgumentParser(description='description:\n\
  Tool for merging multiple .h5 PICLas state files (MyProject_State_000.000000*.h5) containing different hdf5 containers into a single file.\n\
  Supply a single state file or a group of state files by using the wildcard "*", e.g. MyProject_000.00* for a list of file names.\n\n\
  Note the different modes that will be selected automatically depending on the name of the h5 container:\n\n\
    * \'%s\', \'%s\', \'%s\' and \'%s\' are added together\n\
    * \'%s\' container will be concatenated (not added together). To use this feature, supply [-p], [--particles]\n\
    * \'%s\' container in *_TimeAvg_*.h5 files will be averaged over the number of supplied files ([-a], [--average] is automatically turned on)\n\n\
  Note that different types of files, e.g. "_DSMCSurfState_" and "_State_", are not allowed to be mixed. Choose one type!' % (green("SurfaceData"), green("ElemData"), green("DG_Source"), green("DG_Solution"), green("PartData"), green("DG_Solution"))
,formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('files', type=str, help='files (.h5) that are to be merged together.', nargs='+')
parser.add_argument('-a', '--average', action='store_true', help='average the values in the containers by dividing the sum of the arrays by the number of arrays (default=OFF)')
parser.add_argument('-p', '--particles', action='store_true', help='Concatenate \'PartData\' containers (default=OFF)')

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
mintime=1.0e99
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
        mintime=min(time,mintime)
        if time == maxtime :
            maxtimestr = timestr
            maxtimeFile = statefile
            newFile = re.sub(timestr+'.h5', '', statefile)+maxtimestr+'_merged.h5'
        if time == mintime :
            mintimestr = timestr
        files.append(statefile)
    except :
        print("not considering "+statefile)
print("t_min     : %s" % mintime)
print("t_min_str : %s" % mintimestr)

print("t_max     : %s" % maxtime)
print("t_max_str : %s" % maxtimestr)

dt = maxtime-mintime
print("delta t   : %s" % dt)
print("newfile   : %s" % newFile)
print()

# Get maximum number of characters in h5 file names
max_length=0
AlreadyFound = False
for statefile in files :
    if '_TimeAvg_' in statefile and not AlreadyFound :
        print(green("Detected '_TimeAvg_' state file. Switching [-a], [--average] on.\nThis means that the arrays will be averaged after being summed up (by dividing by the total number of files)."))
        args.average = True
        AlreadyFound = True
    max_length = max(max_length,len(statefile))
print(132*"-")
datasets = GetDataSets(files[0])
print("Found the following data sets:\n",datasets)

print(132*"-")
s="Example.h5"
print("".ljust(max_length-len(s)),s," | PartData(dim1, dim2)")
print(132*"-")
n = 0

merged = False
foundPartData = False
foundDG_Solution = False
foundDG_Source = False
foundElemData = False
foundSurfaceData = False
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


    # --------------------------------------------------------------------------------------------------------------------------------------
    # 'PartData' container (concatenate)
    # --------------------------------------------------------------------------------------------------------------------------------------
    if args.particles :
        merged = True
        foundPartData = True
        file_version  = f1.attrs.get('File_Version', default=-1.)[0]
        
        # Usage:
        # -------------------
        # available keys         : print("Keys: %s" % f1.keys())                                # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
        # first key in list      : a_group_key = list(f1.keys())[0]                             # yields 'DG_Solution'
        # available attributes   : print('\n'.join(x for x in f1.attrs))                        # yields 'File_Type\n File_Version\n MeshFile'
        # get specific attribute : file_version  = f1.attrs.get('File_Version', default=-1.)[0] # yields, e.g., 1.5
        # -------------------
        data_set = 'PartData'

        if not data_set in list(f1.keys()) :
            print(red("ERROR: Dataset ['%s'] not found in file. Stop." % data_set))
            exit(0)
        
        # 1.1.1   Read the dataset from the hdf5 file
        PartData = f1[data_set][:]
        print("".ljust(max_length-len(statefile)),statefile," | %s%s" % (data_set, str(PartData.shape)))

        # Save old file
        if n > 1 :
            # Concatenate files depending on the file version (beginning with v1.5 PartData's dimensions are switched)
            if file_version < 1.5 :
                # Compare shape of the dataset of both files, throw error if they do not coincide
                if PartData.shape[0] != PartData_old_shape[0] : # e.g.: PartData.shape = (48, 32)
                    s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\nThe dimensions dim1 = %s and dim2 = %s must be equal!\n\nAborted!"\
                            % (statefile,statefile_old,PartData.shape,PartData_old_shape,PartData.shape[0],PartData_old_shape[0])
                    print(red(s))
                    exit(1)

                # Concatenate columns
                PartData_merged = np.concatenate((PartData_merged, PartData), axis=1)
            else :
                # Compare shape of the dataset of both files, throw error if they do not coincide
                if PartData.shape[1] != PartData_old_shape[1] : # e.g.: PartData.shape = (48, 32)
                    s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\nThe dimensions dim1 = %s and dim2 = %s must be equal!\n\nAborted!"\
                            % (statefile,statefile_old,PartData.shape,PartData_old_shape,PartData.shape[1],PartData_old_shape[1])
                    print(red(s))
                    exit(1)
                # Concatenate rows
                PartData_merged = np.concatenate((PartData_merged, PartData), axis=0)
        else :
            PartData_merged = PartData


        # store the old shape for checking with next state file
        PartData_old_shape = PartData.shape

    # --------------------------------------------------------------------------------------------------------------------------------------
    # 'DG_Solution' container (sum and maybe average if it is a *_TimeAvg_*.h5 file)
    # --------------------------------------------------------------------------------------------------------------------------------------
    data_set = 'DG_Solution'
    if data_set in list(f1.keys()) :
        merged = True
        foundDG_Solution = True
        
        # 1.1.1   Read the dataset from the hdf5 file
        DG_Solution = f1[data_set][:]
        print("".ljust(max_length-len(statefile)),statefile," | %s%s" % (data_set, str(DG_Solution.shape)))

        # Save old file
        if n > 1 :
            # Compare shape of the dataset of both files, throw error if they do not conincide
            if DG_Solution.shape != DG_Solution_old_shape : # e.g.: DG_Solution.shape = (48, 1, 1, 32)
                s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\n\nAborted!"\
                        % (statefile, statefile_old, DG_Solution.shape, DG_Solution_old_shape)
                print(red(s))
                exit(1)
        
            # Add ne current array stored in DG_Solution to b
            DG_Solution_merged = DG_Solution_merged + DG_Solution
        else :
            # Create empty array when processing the first file
            DG_Solution_merged = DG_Solution

        # store the old shape for checking with next state file
        DG_Solution_old_shape = DG_Solution.shape

    # --------------------------------------------------------------------------------------------------------------------------------------
    # 'DG_Source' container (sum and maybe average if it is a *_TimeAvg_*.h5 file)
    # --------------------------------------------------------------------------------------------------------------------------------------
    data_set = 'DG_Source'
    if data_set in list(f1.keys()) :
        merged = True
        foundDG_Source = True
        
        # 1.1.1   Read the dataset from the hdf5 file
        DG_Source = f1[data_set][:]
        print("".ljust(max_length-len(statefile)),statefile," | %s%s" % (data_set, str(DG_Source.shape)))

        # Save old file
        if n > 1 :
            # Compare shape of the dataset of both files, throw error if they do not conincide
            if DG_Source.shape != DG_Source_old_shape : # e.g.: DG_Source.shape = (48, 1, 1, 32)
                s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\n\nAborted!"\
                        % (statefile, statefile_old, DG_Source.shape, DG_Source_old_shape)
                print(red(s))
                exit(1)
        
            # Add ne current array stored in DG_Source to b
            DG_Source_merged = DG_Source_merged + DG_Source
        else :
            # Create empty array when processing the first file
            DG_Source_merged = DG_Source

        # store the old shape for checking with next state file
        DG_Source_old_shape = DG_Source.shape

    # --------------------------------------------------------------------------------------------------------------------------------------
    # 'ElemData' container (sum and maybe average if it is a *_TimeAvg_*.h5 file)
    # --------------------------------------------------------------------------------------------------------------------------------------
    data_set = 'ElemData'
    if data_set in list(f1.keys()) :
        merged = True
        foundElemData = True
        
        # 1.1.1   Read the dataset from the hdf5 file
        ElemData = f1[data_set][:]
        print("".ljust(max_length-len(statefile)),statefile," | %s%s" % (data_set, str(ElemData.shape)))

        # Save old file
        if n > 1 :
            # Compare shape of the dataset of both files, throw error if they do not conincide
            if ElemData.shape != ElemData_old_shape : # e.g.: ElemData.shape = (48, 1, 1, 32)
                s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\n\nAborted!"\
                        % (statefile, statefile_old, ElemData.shape, ElemData_old_shape)
                print(red(s))
                exit(1)
        
            # Add ne current array stored in ElemData to b
            ElemData_merged = ElemData_merged + ElemData
        else :
            # Create empty array when processing the first file
            ElemData_merged = ElemData

        # store the old shape for checking with next state file
        ElemData_old_shape = ElemData.shape

    # --------------------------------------------------------------------------------------------------------------------------------------
    # 'SurfaceData' container (sum and maybe average if it is a *_TimeAvg_*.h5 file)
    # --------------------------------------------------------------------------------------------------------------------------------------
    data_set = 'SurfaceData'
    if data_set in list(f1.keys()) :
        merged = True
        foundSurfaceData = True
        
        # 1.1.1   Read the dataset from the hdf5 file
        SurfaceData = f1[data_set][:]
        print("".ljust(max_length-len(statefile)),statefile," | %s%s" % (data_set, str(SurfaceData.shape)))

        # Save old file
        if n > 1 :
            # Compare shape of the dataset of both files, throw error if they do not conincide
            if SurfaceData.shape != SurfaceData_old_shape : # e.g.: SurfaceData.shape = (48, 1, 1, 32)
                s="\nDatasets are not compatible due to different shapes: Files [%s] and [%s] have shapes %s and %s\n\nAborted!"\
                        % (statefile, statefile_old, SurfaceData.shape, SurfaceData_old_shape)
                print(red(s))
                exit(1)
        
            # Add ne current array stored in SurfaceData to b
            SurfaceData_merged = SurfaceData_merged + SurfaceData
        else :
            # Create empty array when processing the first file
            SurfaceData_merged = SurfaceData

        # store the old shape for checking with next state file
        SurfaceData_old_shape = SurfaceData.shape




    statefile_old = statefile # store name when checking the dimensions in the next file
    f1.close()

# Check if anything was merged
if merged :
    print(132*"-")
    print("Finished reading data from %s files." % n)
    n = float(n)
    
    # Copy old file and modify PartState in the new file
    shutil.copyfile(maxtimeFile, newFile)
    # Write new .h5 file
    f1 = h5py.File(newFile,'r+')
    
    
    if foundPartData:
        data_set = 'PartData'
        del f1[data_set]
        # file: Create new dataset
        PartDataSet = f1.create_dataset(data_set, shape=PartData_merged.shape, dtype=np.float64)
        # Check if average is required
        if args.average:
            print(green("Averaged 'PartData' by dividing by %s" % n))
            PartData_merged = PartData_merged / n
        # write as C-continuous array via np.ascontiguousarray()
        PartDataSet.write_direct(np.ascontiguousarray(PartData_merged))
        print("Files have been merged into %s | %s%s" % (newFile,data_set,PartData_merged.shape))
    if foundDG_Solution:
        data_set = 'DG_Solution'
        del f1[data_set]
        # file: Create new dataset
        DG_SolutionSet = f1.create_dataset(data_set, shape=DG_Solution_merged.shape, dtype=np.float64)
        # Check if average is required
        if args.average:
            print(green("Averaged 'DG_Solution' by dividing by %s" % n))
            DG_Solution_merged = DG_Solution_merged / n
        # write as C-continuous array via np.ascontiguousarray()
        DG_SolutionSet.write_direct(np.ascontiguousarray(DG_Solution_merged))
        print("Files have been merged into %s | %s%s" % (newFile,data_set,DG_Solution_merged.shape))
    if foundDG_Source:
        data_set = 'DG_Source'
        del f1[data_set]
        # file: Create new dataset
        DG_SourceSet = f1.create_dataset(data_set, shape=DG_Source_merged.shape, dtype=np.float64)
        # Check if average is required
        if args.average:
            print(green("Averaged 'DG_Source' by dividing by %s" % n))
            DG_Source_merged = DG_Source_merged / n
        # write as C-continuous array via np.ascontiguousarray()
        DG_SourceSet.write_direct(np.ascontiguousarray(DG_Source_merged))
        print("Files have been merged into %s | %s%s" % (newFile,data_set,DG_Source_merged.shape))
    if foundElemData:
        data_set = 'ElemData'
        del f1[data_set]
        # file: Create new dataset
        ElemDataSet = f1.create_dataset(data_set, shape=ElemData_merged.shape, dtype=np.float64)
        # Check if average is required
        if args.average:
            print(green("Averaged 'ElemData' by dividing by %s" % n))
            ElemData_merged = ElemData_merged / n
        # write as C-continuous array via np.ascontiguousarray()
        ElemDataSet.write_direct(np.ascontiguousarray(ElemData_merged))
        print("Files have been merged into %s | %s%s" % (newFile,data_set,ElemData_merged.shape))
    if foundSurfaceData:
        data_set = 'SurfaceData'
        del f1[data_set]
        # file: Create new dataset
        SurfaceDataSet = f1.create_dataset(data_set, shape=SurfaceData_merged.shape, dtype=np.float64)
        # Check if average is required
        if args.average:
            print(green("Averaged 'SurfaceData' by dividing by %s" % n))
            SurfaceData_merged = SurfaceData_merged / n
        # write as C-continuous array via np.ascontiguousarray()
        SurfaceDataSet.write_direct(np.ascontiguousarray(SurfaceData_merged))
        print("Files have been merged into %s | %s%s" % (newFile,data_set,SurfaceData_merged.shape))
    

    print(132*"-")
    
    
    
    
    
    f1.close()
