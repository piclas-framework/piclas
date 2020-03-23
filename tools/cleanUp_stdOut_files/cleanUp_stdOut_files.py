import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil

# Bind raw_input to input in Python 2
try:
    input = raw_input
except NameError:
    pass

def getFirst(line):
    first = line.strip("\n").split()
    if first:
        return first[0]
    else:
        return None

def getSecond(line):
    line = line.strip("\n").split()
    try:
        second = line[1]
    except Exception as e:
        print(e)
        second = None

    return second

def getPartInfo(line):
    #print(line.strip("\n"))
    #print(line.strip("\n").split())
    line = line.strip("\n").split()
    PartID = 0
    Element = 0
    SpecID = 0
    #print(line)

    NbrOfLost_Index = line.index('lost.')
    Element_Index = line.index('Element:')
    SpecID_Index = line.index('(species:')

    PartID = int(line[NbrOfLost_Index-1])
    Element = int(line[Element_Index+1])
    SpecID = int(line[SpecID_Index+1])

            #['32', 'Error', 'in', 'Particle', 'TriaTracking!', 'Particle', 'Number', '442065', 'lost.', 'Element:', '2881', '(species:', '1', ')']
    return PartID, Element, SpecID 


def CleanFile(stdfile) :

    killList = {}

    n=0
    nLostParts=0
    meshFound=False
    with open(stdfile, "r") as input:
        with open(stdfile+".lost", "w") as output_lost: 
            with open(stdfile+".new", "w") as output_new: 
                for line in input:
                    n+=1
                    first = getFirst(line)

                    if 'Error in Particle TriaTracking! Particle Number' in line.strip("\n"):
                        nLostParts+=1
                        if first in killList:
                            print("Error: the rank is already in the list. The lines might overlap. Resolve this in \nLine : %s" % n)
                            exit(1)
                        killList.update( {first : 1} )
                        output_lost.write(line)
                    elif first in killList:
                        output_lost.write(line)
                        killList[first] += 1
                        if killList[first] == 5:
                            del killList[first]
                    else:
                        output_new.write(line)
                        if not meshFound:
                            if 'MeshFile' in line.strip("\n"):
                                meshFound=True
                                myline=line.replace(" ", "")
                                tmp=myline.split('|')
                                for x in tmp:
                                    if x.endswith('.h5'):
                                        MeshFile=x

                #print(red("Lost %s particles" % nLostParts))

    killList = {}

    # Read std.out.lost file with all lost particles and write the data to a .h5 file
    if nLostParts > 0:
        # 1. Open .h5 file
        f1 = h5py.File(stdfile+"-lost-particles.h5",'w')

        data = np.zeros(( nLostParts,12))

        # 2. Read lost particles file
        with open(stdfile+".lost", "r") as f:
            lines = f.readlines()

        # 3. Read lost particles and sort particle data into array
        n=0
        for line in lines:
            first = getFirst(line)

            if 'Error in Particle TriaTracking! Particle Number' in line.strip("\n"):
                if first in killList:
                    print("Error: the rank is already in the list. The lines might overlap. Resolve this in \nLine : %s" % n)
                    exit(1)

                killList.update( {first : 1} )
                PartID, Element, SpecID = getPartInfo(line)
                #print(PartID, Element, SpecID)

            elif first in killList:
                second = getSecond(line)
                linesplit = line.strip("\n").split()
                #print('[%s]' % second)
                if second == 'LastPos:':
                    xLastPos = float(linesplit[linesplit.index('LastPos:')+1])
                    yLastPos = float(linesplit[linesplit.index('LastPos:')+2])
                    zLastPos = float(linesplit[linesplit.index('LastPos:')+3])
                elif second == 'Pos:':
                    x = float(linesplit[linesplit.index('Pos:')+1])
                    y = float(linesplit[linesplit.index('Pos:')+2])
                    z = float(linesplit[linesplit.index('Pos:')+3])
                elif second == 'Velo:':
                    vx = float(linesplit[linesplit.index('Velo:')+1])
                    vy = float(linesplit[linesplit.index('Velo:')+2])
                    vz = float(linesplit[linesplit.index('Velo:')+3])

                killList[first] += 1
                if killList[first] == 5:
                    #print(n,x,y,z,vx,vy,vz,Element,SpecID,PartID,xLastPos,yLastPos,zLastPos)
                    data[n,0] = x
                    data[n,1] = y
                    data[n,2] = z
                    data[n,3] = vx
                    data[n,4] = vy
                    data[n,5] = vz
                    data[n,6] = SpecID
                    data[n,7] = Element
                    data[n,8] = PartID 
                    data[n,9] = xLastPos
                    data[n,10] = yLastPos
                    data[n,11] = zLastPos

                    n+=1

                    del killList[first]

            elif 'Warning: Particle located inside of face and moves parallel to side. Undefined position.' in line.strip("\n"):
                print('Warning: Particle located inside of face and moves parallel to side. Undefined position.')
                print('Found lost particle within tracing method. This has not yet been implemented. Please contact the developer!.')
                exit(1)

            elif 'Tolerance issue during tracing! Unable to locate particle inside computational domain' in line.strip("\n"):
                print('Tolerance issue during tracing! Unable to locate particle inside computational domain')
                print('Found lost particle within tracing method. This has not yet been implemented. Please contact the developer!.')
                exit(1)


        # 4. Write dummy DG_Solution container
        data1 = np.zeros(( 0,0))
        dset1 = f1.create_dataset('DG_Solution', shape=data1.shape, dtype=np.float64)

        # 5. Create new dataset 'dset'
        dset = f1.create_dataset('PartData', shape=data.shape, dtype=np.float64)
        #dset = f1.create_dataset(data, shape=data.shape, dtype=getattr(np, str(dataType)))

        # 6. Write as C-continuous array via np.ascontiguousarray()
        if not data.any() :
            print(" %s has dimension %s. Skipping" % (data_set,data.shape))
            pass
        else :
            dset.write_direct(np.ascontiguousarray(data))

        # 7. Write attributes
        f1.attrs.modify('File_Version',[1.5])   # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('NodeType',[b'GAUSS'])  # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('File_Type',[b'PartStateBoundary']) # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('Program',[b'PICLas'])  # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('VarNames',[b'empty'])  # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('Time',[0.])            # these brackets [] are required for ParaView plugin !
        f1.attrs.create('MeshFile', [MeshFile], None, dtype='<S255')
        f1.attrs.create('N', [1], None, dtype='i4')
        f1.attrs.create('Project_Name',[b'Lost-particles'], None, dtype='<S255')  # these brackets [] are required for ParaView plugin ! 
        f1.attrs.create('VarNamesParticles', [b'ParticlePositionX'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'ParticlePositionY'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'ParticlePositionZ'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'VelocityX'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'VelocityY'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'VelocityZ'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'Species'.ljust(255),    # .ljust(255) is required for ParaView plugin !
                                              b'ElementID'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'PartID'.ljust(255),     # .ljust(255) is required for ParaView plugin !
                                              b'LastPos-X'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'LastPos-Y'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'LastPos-Z'.ljust(255)], None, dtype='<S255')   # .ljust(255) is required for ParaView plugin !

        # 8. Close .h5 data file
        f1.close()


        return nLostParts


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
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for cleaning std*.out files.\nSupply a single file or a group of files by using the wildcard "*", e.g. std* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Files (std*.out) that are to be cleaned.', nargs='+')
parser.add_argument('-d', '--debug', action='store_true', help='Print additional imformation regarding the files onto screen.')

# Get command line arguments
args = parser.parse_args()

# Display all command line arguments
print('='*132)
print("Running with the following command line options")
for arg in list(args.__dict__) :
    print(arg.ljust(15)+" = [ "+str(getattr(args,arg))+" ]")
print('='*132)


# Get maximum number of characters in h5 file names
max_length=0
for stdfile in args.files :
    max_length = max(max_length,len(stdfile))

InitialDataRead = True
NbrOfFiles = len(args.files)

ext = ['.new', '.lost', '.h5']
for stdfile in args.files :
    if stdfile.endswith(tuple(ext)):
        print("%s " % stdfile + yellow("(skipping)"))
        continue

    nLostParts = CleanFile(stdfile)

    if nLostParts > 0:
        print("%s " % stdfile + red("Lost %s particles" % nLostParts) + " Written particles to %s-lost-particles.h5" % stdfile)
    else:
        print("%s" % stdfile)

print(132*"-")
