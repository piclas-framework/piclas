import numpy as np
import math
from timeit import default_timer as timer
import argparse
import re
import shutil
import os.path
import configparser
import types

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


def GetInfoFromDataset(statefile,data_set,NbrOfFiles,species_info_read,CutOff) :
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

        # Create empy data array
        #data = np.empty([NbrOfFiles, 2], dtype=float)
        data = np.zeros([NbrOfFiles, 2+len(Charge)], dtype=float)
        #print("data = %s" % (data))
        #print("Charge = %s" % (Charge))
        #exit(0)

    # 2   Sum the charges in PartData and divide by the time difference to the current
    TotalCharge = 0.
    SumCharge = np.zeros([len(Charge)], dtype=float)
    #print("SumCharge = %s" % (SumCharge))
    for i in range(len(b1)):
        SpecID = int(b1[i][SpeciesIndex])

        if CutOff.Check:
            # Skip coordinate
            if CutOff.Type == 'greater':
                if b1[i][CutOff.Coord] > CutOff.Val:
                    #print(b1[i][0],b1[i][1],b1[i][2]) 
                    continue
            elif CutOff.Type == 'lesser':
                if b1[i][CutOff.Coord] < CutOff.Val:
                    #print(b1[i][0],b1[i][1],b1[i][2]) 
                    continue
            # Skip radius
            if CutOff.Rad is not None:
                if CutOff.Type == 'greater':
                    if np.linalg.norm(b1[i][0:1]) > CutOff.Rad:
                        continue
                elif CutOff.Type == 'lesser':
                    if np.linalg.norm(b1[i][0:1]) < CutOff.Rad:
                        continue

        if SpecID not in Charge :
            txt = input(yellow("    Enter the charge for species %s in multiples of the elementary charge e=1.602176634e-19 : " % SpecID))
            Charge[SpecID] = float(txt)*1.602176634e-19
            print(yellow("    The charges of the species are: %s" % Charge))
            # Add another column to the data array for the new species
            AddData = np.zeros([NbrOfFiles, 1], dtype=float)
            data = np.concatenate((data,AddData),axis=1)
            #print(data)
            # Add another entry for the sum of charge
            AddSumCharge = np.zeros([1], dtype=float)
            #print("AddSumCharge = %s" % (AddSumCharge))
            SumCharge = np.append(SumCharge,[0.])
            #print("SumCharge = %s" % (SumCharge))
            #exit(0)

        myindex = list(Charge.keys()).index(SpecID)
        #print("myindex = %s" % (myindex))
        TotalCharge        += Charge[SpecID] * b1[i][MPFIndex]
        #print("TotalCharge = %s" % (TotalCharge))
        SumCharge[myindex] += Charge[SpecID] * b1[i][MPFIndex]
        #print("SumCharge = %s" % (SumCharge))

    dt                = Time-PreviousTime
    Current           = TotalCharge / dt
    PreviousTime      = Time
    data[FileCount,0] = Time
    data[FileCount,1] = Current
    for SpecID, v in Charge.items():
        myindex = list(Charge.keys()).index(SpecID)
        Current = SumCharge[myindex] / dt
        data[FileCount,myindex+2] = Current

        #print("myindex = %s" % (myindex))

    FileCount        += 1

    # 3   Close .h5 data file
    f1.close()

    print("%s  t=%.2E  dt=%.2E  charge=%.2E  current=%.2E"% (statefile, Time, dt, TotalCharge, Current))

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

# Start the timer
start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for adding up the charges of the particles in the PartState container of multiple .h5 files.\nSupply a single state file or a group of state files by using the wildcard "*", e.g. MySimulation_000.000000* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Files (.h5) that are to be merged together.', nargs='+')
parser.add_argument('-d', '--debug', action='store_true', help='Print additional information regarding the files onto screen.')
parser.add_argument('-p', '--parameter', default='parameter.ini', help='Name of the parameter file that contains the species information corresponding to the PartState containers. Default=parameter.ini')
parser.add_argument('-c', '--cutoff', action='store_true', help='Exclude particles by given a cut-off coordinate and removing particles that are either > or < as compared with this value.')
parser.add_argument('-o', '--coordinate', type=int, help='Cut-off coordinate 0: x-dir, 1: y-dir, 2: z-dir.')
parser.add_argument('-v', '--value', type=float, help='Cut-off value for removing particles.')
parser.add_argument('-r', '--radius', type=float, help='Cut-off value for removing particles with radius grater/lesser than the supplied value (uses -t, --type).')
parser.add_argument('-t', '--type', type=str, help='Cut-off type. Choose "greater" or "lesser" to remove particles with a coordinate > or < as compared with the cut-off value.')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
print('='*132)
print("Running with the following command line options")
for arg in list(args.__dict__) :
    print(arg.ljust(15)+" = [ "+str(getattr(args,arg))+" ]")
print('='*132)

# Check cut-off parameters
CutOff= types.SimpleNamespace()
CutOff.Check = args.cutoff
if args.cutoff:
    CutOff.Coord = args.coordinate
    CutOff.Val   = args.value
    CutOff.Rad   = args.radius
    CutOff.Type  = args.type
    if not CutOff.Coord in [0,1,2]:
        print(red("Cut-Off: Coordinate [-o] must either be 0, 1 or 2."))
        exit(0)
    if not CutOff.Type in ['greater','lesser']:
        print(red('Cut-Off: Type [-t] must either be "greater" or "lesser".'))
        exit(0)
    if CutOff.Val is None and CutOff.Val is None:
        print(red('Cut-Off: Value [-v] and/or Radius [-r] must be chosen.'))
        exit(0)
    mystr='xyz'
    print(red("WARNING: Coordinate cut-off is activated. Removing particles with %s-coordinate %s than %s" % (mystr[CutOff.Coord], CutOff.Type, CutOff.Val)))
    if CutOff.Rad is not None:
        print(red("WARNING: Radius cut-off is activated. Removing particles with radius %s than %s" % (CutOff.Type, CutOff.Rad)))

# Get maximum number of characters in h5 file names
max_length=0
for statefile in args.files :
    max_length = max(max_length,len(statefile))

#parameterFile="parameter.ini"
config = ReadConfig(args.parameter)

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
    print(yellow("    Automatically fetched species ID's and charges from %s" % args.parameter))
    print(yellow("    The charges of the species are: %s" % Charge))
else:
    species_info_read=False

if not species_info_read:
    print(yellow("Parameter file [%s] not found. Please enter the required species information by hand or create a symbolic link to the parameter.ini file which was used to create the data" % args.parameter))

InitialDataRead = True
NbrOfFiles = len(args.files)
for statefile in args.files :
    if statefile.endswith('.h5') and not statefile.endswith('_merged.h5'):
        # check if file exists
        if not os.path.exists(statefile):
            print(red('File does not exist: [%s]' % statefile))
            exit(0)
    else:
        NbrOfFiles -= 1


print("Reading %s files" % NbrOfFiles)
for statefile in args.files :
    if statefile.endswith('.h5') and not statefile.endswith('_merged.h5'):

        # Read data from .h5 file
        GetInfoFromDataset(statefile,'PartData',NbrOfFiles,species_info_read,CutOff)
    else:
        print(statefile," skipping")

#print( )
#print("Charge = %s" % (Charge))
#print("len(Charge) = %s" % (len(Charge)))
#x = ',"'.join(myTuple)
#print("x = %s" % (x))

myheader = '"time","Current [A]"'
for k, v in Charge.items():
    #print(k,v)
    #print(list(Charge.keys()).index(k))

    #print( )
    myheader = myheader+str(',"Spec-(%s)"' % k)

print( )
print("myheader = %s" % (myheader))
print(132*"-")
#print(data)
if 'data' in locals():
    print("The charge data has been written to CurrentOverTime_PartStateBoundary.csv")
    np.savetxt('CurrentOverTime_PartStateBoundary.csv', data, delimiter=',', header=myheader, comments='') # comments='' is required to prevent the "#" from being written in the header line
else:
    print("No output created.")
print(132*"-")

