import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil
import os.path
import configparser

# Bind raw_input to input in Python 2
try:
    input = raw_input
    import ConfigParser as ConfPars
except NameError:
    import configparser as ConfPars
    pass

def CreateConfig(parameterFile):
    if os.path.exists(parameterFile):
        with open(parameterFile, "r") as f:
            lines = f.readlines()
            
        config = {}
        for line in lines:
            line=line.replace(" ", "").strip("\n")
            if not line.startswith("!"):
                idx = line.find('!')
                if idx > -1:
                    line=line[:idx]
                if line:
                    keyval=line.split("=")
                    key=keyval[0]
                    val=keyval[1]

                    # check if the same parameter name (e.g. 'BoundaryName') occurs more than once in the list and
                    # move multiple occurances to a separate key/value where the value is a list of all occurances
                    # this must be done, because dicts cannot have the same key name more than once (it is a dictionary)
                    found, number = isKeyOf(config,key)
                    if found :
                        if type(config[key]) == tuple:
                            config[key] = config[key] + tuple([val])
                        else:
                            config[key] = tuple([config[key]]) + tuple([val])
                    else :
                        config[key] = val

        with open(parameterFile+".new", "w") as output_new: 
            output_new.write("[Section1]\n")
            for x in config.keys():
                if type(config[x]) == tuple:
                    #print(x.ljust(40) +" = " + ", ".join(str(y) for y in config[x]))
                    output_new.write(x.ljust(40) +" = " + ", ".join(str(y) for y in config[x])+"\n")
                else:
                    #print(x.ljust(40) +" = " + config[x].ljust(25))
                    output_new.write(x.ljust(40) +" = " + config[x].ljust(25)+"\n")
        exit_code=0
    else:
        exit_code=-1

    return exit_code

def isKeyOf(a,key_IN) :
    """Check if the dictionary 'a' contains a key 'key_IN'"""
    found = False
    number = 0
    for key in a.keys() :
        if key == key_IN :
            number += 1
            found = True
    return found, number

def ReadConfig(parameterFile):
    # Create config file that is conform with the format of ConfigParser()
    exit_code = CreateConfig(parameterFile)

    # Read config with ConfigParser()
    if exit_code == 0:
        config =  ConfPars.ConfigParser()
        config.read(parameterFile+".new")
    else:
        config=None

    return config


def GetInfoFromDataset(statefile,data_set,NbrOfFiles,species_info_read) :
    global InitialDataRead
    global SpeciesIndex
    global MPFIndex
    global Charge
    global PreviousTime
    global data
    global FileCount

    #if not InitialDataRead :
        #return

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
    
    # 1   Read the dataset from the hdf5 file
    try :
        b1 = f1[data_set][:]
    except :
        print('Dataset %s does not exist' % data_set)
        return

    dataType = f1[data_set].dtype

    if args.debug:
        print( )
        print(yellow("Keys: %s" % f1.keys()))
        print(yellow("    size :         "+ str(f1[data_set].size)))
        print(yellow("    shape :        "+ str(f1[data_set].shape)))
        print(yellow("    dtype :        "+ str(f1[data_set].dtype)))
        print(yellow("    chunks :       "+ str(f1[data_set].chunks)))
        print(yellow("    compression :  "+ str(f1[data_set].compression)))
        print(yellow("    shuffle :      "+ str(f1[data_set].shuffle)))


    Time = f1.attrs.get('Time', default=-1.)
    Time = Time[0]

    # 2   Read the attributes and determine the species charges from user input
    if InitialDataRead :
        PreviousTime = 0.0
        data = np.empty([NbrOfFiles, 2], dtype=float)
        FileCount=0
        InitialDataRead = False

        # Attributes
        Attributes=[]
        for x in f1.attrs :
            Attributes.append("%s" % x)


        # VarNamesParticles
        f1VarNamesParticles = f1.attrs.get('VarNamesParticles', default=-1.)
        s = ' '.join('%s' % x for x in f1VarNamesParticles)
        l = s.split("'")[1::2];

        VarNamesParticles = []
        for x in l :
            VarNamesParticles.append("%s" % x)
        SpeciesIndex = VarNamesParticles.index('Species')
        MPFIndex = VarNamesParticles.index('MacroParticleFactor')

        # Only request user input for species if this info could not be read automatically from parameter.ini
        if not species_info_read:
            Species=[]
            for i in range(len(b1)):
                SpecID = b1[i][SpeciesIndex]
                if SpecID not in Species :
                    Species.append(SpecID)
            Species = [x+1 for x in range(int(max(Species)))]

            #for i in range(len(b1)):
                #for j in range(len(b1[i])):
                    #print(b1[i][j])

            Charge={}
            for iSpec in Species :
                txt = input(yellow("    Enter the charge for species %s in multiples of the elementary charge e=1.602176634e-19 : " % iSpec))
                Charge[iSpec] = float(txt)*1.602176634e-19
            print(yellow("    The charges of the species are: %s" % Charge))


    # 2   Sum the chares in PartData and divide by the time difference to the current
    TotalCharge = 0.
    for i in range(len(b1)):
        SpecID = int(b1[i][SpeciesIndex])
        if SpecID not in Charge :
            txt = input(yellow("    Enter the charge for species %s in multiples of the elementary charge e=1.602176634e-19 : " % SpecID))
            Charge[SpecID] = float(txt)*1.602176634e-19
            print(yellow("    The charges of the species are: %s" % Charge))
        TotalCharge += Charge[SpecID] * b1[i][MPFIndex]

    dt=Time-PreviousTime
    Current = TotalCharge / dt
    PreviousTime = Time
    data[FileCount,0] = Time
    data[FileCount,1] = Current
    FileCount += 1

    # 3   Close .h5 data file
    f1.close()

    print("%s  dt=%.2E  charge=%.2E  current=%.2E"% (statefile, dt, TotalCharge, Current))

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
    print(red('Could not import h5py module. This is required for anaylze functions.'))
    exit(0)

# Start the timer
start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for adding up the charges of the particles in the PartState container of multiple .h5 files.\nSupply a single state file or a group of state files by using the wildcard "*", e.g. MySimulation_000.000000* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
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
#print()

# Get maximum number of characters in h5 file names
max_length=0
for statefile in args.files :
    max_length = max(max_length,len(statefile))

parameterFile="parameter.ini"
config = ReadConfig(parameterFile)

if config :
    species_info_read=True

    nSpecies = config.get("Section1","Part-nSpecies")

    # Set species ID's
    Species=[]
    Species = [x+1 for x in range(int(max(nSpecies)))]

    # Set charge of species
    Charge={}
    for iSpec in Species :
        Charge[iSpec] = float(config.get("Section1","Part-Species%d-ChargeIC" % iSpec))
    print(yellow("    Automatically fetched species ID's and charges from %s" % parameterFile))
    print(yellow("    The charges of the species are: %s" % Charge))
else:
    species_info_read=False

InitialDataRead = True
NbrOfFiles = len(args.files)
#print(NbrOfFiles)
for statefile in args.files :
    if statefile.endswith('.h5'):
        #print(statefile)
        # Flip PartData
        GetInfoFromDataset(statefile,'PartData',NbrOfFiles,species_info_read)
    else:
        print(statefile," skipping")

print(132*"-")
#print(data)
print("The charge data has been written to CurrentOverTime_PartStateBoundary.csv")
np.savetxt('CurrentOverTime_PartStateBoundary.csv', data, delimiter=',', header='"time", "Current [A]"')
print(132*"-")

