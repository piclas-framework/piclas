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
    #raise ImportError('Could not import h5py module. This is required for analyse functions.')
    print(tools.red('Could not import h5py module. This is required for analyse functions.'))
    exit(0)


start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for clipping a bounding box of the PartState by supplying a single state file.\nNote that a new file is created and its size may be very large as compared with the original file.\nUsage example: python visuBox_PartState.py 2Dplasma_detail_State_000.00000144247600000.h5 -x 0.4 0.45 -y 0.4 0.45 -z 0.0 1.00\n\nARGUMENTS:', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-x', '--xdir', type=str, default=['-1.0e99', '1.0e99'], help='x- y- and z-coordinates.', nargs='+')
parser.add_argument('-y', '--ydir', type=str, default=['-1.0e99', '1.0e99'], help='x- y- and z-coordinates.', nargs='+')
parser.add_argument('-z', '--zdir', type=str, default=['-1.0e99', '1.0e99'], help='x- y- and z-coordinates.', nargs='+')
parser.add_argument('statefile', type=str, help='File (.h5) that contains a PartState container.')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
print("Running with the following command line options")
for arg in list(args.__dict__) :
    print(arg.ljust(15)+" = [ "+str(getattr(args,arg))+" ]")
print('='*132)

#if len(args.xyz) != 6:
    #print("supply exactly 6 coordinates for xyz")
    #exit(0)

xmin = float(args.xdir[0])
xmax = float(args.xdir[1])
ymin = float(args.ydir[0])
ymax = float(args.ydir[1])
zmin = float(args.zdir[0])
zmax = float(args.zdir[1])

print("x-dir: %s, %s" % (xmin,xmax))
print("y-dir: %s, %s" % (ymin,ymax))
print("z-dir: %s, %s" % (zmin,zmax))

print(args.statefile)
pat = r'^.*\_(.*)\..*$'
match = re.search(pat, args.statefile)
timestr = match.group(1)
try :
    time = float(timestr)
    newFile = re.sub(timestr+'.h5', '', args.statefile)+timestr+'_visuBox.h5'
except :
    print("something wrong with "+args.statefile)
    exit(0)

print("t_max_str : %s" % timestr)
print("newfile   : %s" % newFile)
# Open h5 file and read container info
# --------------------------------------------
#     r       : Readonly, file must exist
#     r+      : Read/write, file must exist
#     w       : Create file, truncate if exists
#     w- or x : Create file, fail if exists
#     a       : Read/write if exists, create otherwise (default
# --------------------------------------------
# When sorting is used, the sorted array is written to the original .h5 file with a new name
f1 = h5py.File(args.statefile,'r+')

# Usage:
# -------------------
# available keys   : print("Keys: %s" % f1.keys())         # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
# first key in list: a_group_key = list(f1.keys())[0]      # yields 'DG_Solution'
# -------------------
data_set1= 'PartData'
#data_set2= 'DG_Solution'

# 1.1.1   Read the dataset from the hdf5 file
b1 = f1[data_set1][:]

print(132*"-")
print("Original dataset")
print(b1)
print(b1.shape)
print(b1.shape[1])

# Find points outside of box and mark with NaN
for column in range(b1.shape[1]):
    # Check x-direction
    if b1[0,column] < xmin:
        b1[:,column] = np.nan
        continue
    if b1[0,column] > xmax:
        b1[:,column] = np.nan
        continue
    
    # Check y-direction
    if b1[1,column] < ymin:
        b1[:,column] = np.nan
        continue
    if b1[1,column] > ymax:
        b1[:,column] = np.nan
        continue
    
    # Check z-direction
    if b1[2,column] < zmin:
        b1[:,column] = np.nan
        continue
    if b1[2,column] > zmax:
        b1[:,column] = np.nan
        continue

# Delete all columns with NaN
b1 = b1[:,~np.all(np.isnan(b1), axis=0)]

print(132*"-")
print("New dataset")
print(b1)
if not b1.any() : # The truth value of an array with more than one element is ambiguous.
    print("Resulting array is empty!")
    exit(0)
print(b1.shape)
print(132*"-")
# Copy old file and modify PartState in the new file
shutil.copyfile(args.statefile, newFile)

# Write new .h5 file
f2 = h5py.File(newFile,'r+')

# Delete 'PartData' in new file in order to replace it with the reduced array
del f2[data_set1]

#for key, prm in f1.attrs.items() :
    #f2.attrs[key]=prm

# close old file
f1.close()

# file: Create new dataset
dset = f2.create_dataset(data_set1, shape=b1.shape, dtype=np.float64)
# write as C-continuous array via np.ascontiguousarray()
dset.write_direct(np.ascontiguousarray(b1))

# file: Create new dataset
#dset = f2.create_dataset(data_set2, shape=b2.shape, dtype=np.float64)
# write as C-continuous array via np.ascontiguousarray()
#dset.write_direct(np.ascontiguousarray(b2))

f2.close()
