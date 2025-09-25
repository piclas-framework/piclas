import time
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, draw, show
import math

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
    print(red('Could not import h5py module. This is required for reading .h5 files. Exit.'))
    exit(0)


def NxM(num):
    """convert num to N x M matrix

    :num: TODO
    :returns: TODO

    """
    base = int(round(math.sqrt(num)))
    remainder = int(num-float(base*base))

    col = base
    row = base + min(1,max(0,remainder))

    return row, col

def createPlot(title,NbrOfAttributes,VarNames,NbrOfRP,t,data):
    global fig
    global ax
    global numbers
    # Calculate the number of required plots (assume that the first attribute is the time)
    NbrOfPlots = NbrOfAttributes-1
    # Get the number of rows and columns
    row, col = NxM(NbrOfPlots)
    # Check if the figure has already been created
    if not plt.fignum_exists(1):
        fig, ax = plt.subplots(nrows=row, ncols=col)
        ls = 'solid'
    else:
        ls = 'dashed'
    # Check if the reference solution is plotted or not
    if 'ref' in title:
        s = 'ref'
    else:
        s = 'not ref'
    # Get the number of figures
    numbers = plt.get_fignums()
    # Add title to the figure
    fig.suptitle('%s' % title, fontsize=16)
    # Initialise number of variables
    iVar=0
    # Loop over rows
    for row in ax:
        # Loop over columns
        for col in row:
            # Count number of variables
            iVar += 1
            # Initialize legend container for each sub-plot
            legend = []
            # Stop is all RPs have been plotted
            if iVar > NbrOfPlots:
                break
            # Loop over the RPs
            for iRP in range(NbrOfRP):
                colplt, = col.plot(t, data[:,iRP,iVar], linestyle=ls, label='RP$_{%s}$ (%s)' % (iRP, s) )
                legend.append(colplt)
            # Display title
            col.set_title(f'%s' % VarNames[iVar-1])
            # Display legend
            col.legend(handles=legend)

    #plt.show()
    draw()

    return fig, ax


def ReadFileAndCreatePlot(statefile):
    """123

    :statefile: TODO
    :returns: TODO

    """
    global fig
    global ax
    global numbers
    start = time.time()

    # Open h5 file and read container info
    # --------------------------------------------
    #     r       : Readonly, file must exist
    #     r+      : Read/write, file must exist
    #     w       : Create file, truncate if exists
    #     w- or x : Create file, fail if exists
    #     a       : Read/write if exists, create otherwise (default
    # --------------------------------------------
    # When sorting isused, the sorted array is written to the original .h5 file with a new name
    f1 = h5py.File(statefile,'r+')
    file_stats = os.stat(statefile)
    if file_stats.st_size > 1024 * 1024:
        if file_stats.st_size > 1024 * 1024 * 1024:
            filesize = '%.2f GB' % round(file_stats.st_size / (1024 * 1024 * 1024), 2)
        else:
            filesize = '%.2f MB' % round(file_stats.st_size / (1024 * 1024), 2)
    else:
        filesize = '%.2f kB' % round(file_stats.st_size / (1024), 2)

    # Usage:
    # -------------------
    # available keys         : print("Keys: %s" % f1.keys())                                # yields, e.g., <KeysViewHDF5 ['DG_Solution', 'PartData']>
    # first key in list      : a_group_key = list(f1.keys())[0]                             # yields 'DG_Solution'
    # available attributes   : print('\n'.join(x for x in f1.attrs))                        # yields 'File_Type\n File_Version\n MeshFile'
    # get specific attribute : file_version  = f1.attrs.get('File_Version', default=-1.)[0] # yields, e.g., 1.5
    # -------------------

    # 1   Read the dataset from the hdf5 file
    dataset='RP_Data'
    try :
        data = f1[dataset][:]
    except :
        print('Dataset %s does not exist' % dataset)
        exit(0)

    dataType = f1[dataset].dtype

    t = data[:,0,0]

    debug=True
    debug=False
    if debug:
        print()
        attributes = ', '.join(x for x in f1.attrs)
        print(yellow("Attributes: %s" % attributes))
        VarNames = f1.attrs.get('VarNames', default="")
        print(yellow("    VarNames: %s" % VarNames))
        print()
        print(yellow("Keys: %s" % f1.keys()))
        print(yellow("    size :         "+ str(f1[dataset].size)))
        print(yellow("    shape :        "+ str(f1[dataset].shape)))
        print(yellow("    dtype :        "+ str(f1[dataset].dtype)))
        print(yellow("    chunks :       "+ str(f1[dataset].chunks)))
        print(yellow("    compression :  "+ str(f1[dataset].compression)))
        print(yellow("    shuffle :      "+ str(f1[dataset].shuffle)))

    VarNames = f1.attrs.get('VarNames', default="")
    NbrOfPoints     = f1[dataset].shape[0]
    NbrOfRP         = f1[dataset].shape[1]
    NbrOfAttributes = f1[dataset].shape[2]
    if debug:
        print("NbrOfPoints     = %s" % (NbrOfPoints))
        print("NbrOfRP         = %s" % (NbrOfRP))
        print("NbrOfAttributes = %s (including time)" % (NbrOfAttributes))

    fig, ax = createPlot(statefile, NbrOfAttributes, [bytesqequence.decode('ascii') for bytesqequence in VarNames], NbrOfRP, t, data)

    f1.close()

    end = time.time()
    elapsed = end-start
    rounded = "{:.2f}".format(round(elapsed, 2))
    s = "----- Required time for %s: %s sec [%s] for %s" % (statefile,rounded,str(datetime.timedelta(seconds=int(elapsed))),filesize)
    print(s)



# ============================================================
# Program starts here
# ============================================================

statefiles=['twt_RP_000.0000000002000000_reference.h5', 'twt_RP_000.00000000020000000.h5']
for statefile in statefiles:
    ReadFileAndCreatePlot(statefile)

# at the end call show to ensure window won't close.
show()

