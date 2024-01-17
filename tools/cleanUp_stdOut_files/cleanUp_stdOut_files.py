import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil
import os
import grp

import pwd
from pwd import getpwuid

# Bind raw_input to input in Python 2
try:
    input = raw_input
except NameError:
    pass

def hasNumbers(inputString):
    return bool(re.search(r'\d', inputString))

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

def get_owner_and_group(stdfile):
    try:
        # get stat of file/folder or symbolic link
        if os.path.islink(stdfile):
            status = os.lstat(stdfile)
        else:
            status = os.stat(stdfile)

        uid = status.st_uid
        gid = status.st_gid

        user = pwd.getpwuid(uid)[0]
        group = grp.getgrgid(gid)[0]
    except Exception as e:
        return None, None

    return user, uid, group, gid


def RenameFiles(differences, stdfile, stdfile_backup, stdfile_new, args):

    # Check user/group name vs. original file
    userOrig,uidOrig,groupOrig,gidOrig = get_owner_and_group(stdfile)

    # Check if differences exist (nLostParts or changedLines)
    if differences > 0:

        # Only create backup file if it does not exist (i.e. prevent over-writing backup files from other clean-up functions in the tool)
        if not os.path.exists(stdfile_backup):
            os.rename(stdfile, stdfile_backup) # backup original file (only once)
        os.rename(stdfile_new, stdfile)        # replace original file with cleaned file

        # Check group name vs. original file
        user,uid,group,gid = get_owner_and_group(stdfile)
        if groupOrig is not None and group is not None:
            if groupOrig == group:
                pass
            else:
                # Change the owner and group to [original owner] and [original group]
                if gid != gidOrig:
                    os.chown(stdfile, uidOrig, gidOrig, follow_symlinks=False)

    else :
        os.remove(stdfile_new) # remove new file (it is empty when no particles were lost)

    # Check if user has supplied the flag for backup file removal and the file actually exists
    if not args.save and os.path.exists(stdfile_backup):
        os.remove(stdfile_backup) # remove backup file


def CleanSingleLines(stdfile,args):
    '''Remove lines such as
     1. [ REACTION       85084           5           1]
     2. [ CALCULATING THE ADAPTIVE SURFACE FLUX VALUES... Number of sampled particles:   0.0000000000000000     ]
     3. [   3382.0000000000000        1694.4798966228368]
     OPTIONAL 4. [iter:                  702 time:   3.2151600000000616E-008]
    '''

    stdfile_new    = stdfile+".new"
    stdfile_backup = stdfile+".bak"

    # 1. Read the original file content
    with open(stdfile) as input:
        lines = input.readlines()

    # REACTION       85084           5           1
    n=0
    changedLines=0
    with open(stdfile_new, "w") as output_new:
        for line in lines:
            n+=1
            #if all(c in line.strip("\n") for c in arr):
            line_stripped = line.strip("\n")
            line_split = line_stripped.split()
            if line_stripped.startswith(' REACTION') and hasNumbers(line_stripped):
                # 1. [ REACTION       85084           5           1]
                # Ignore this line and increase the counter by 1
                changedLines+=1
            elif line_stripped.startswith(' CALCULATING THE ADAPTIVE SURFACE FLUX VALUES') and hasNumbers(line_stripped):
                # 2. [ CALCULATING THE ADAPTIVE SURFACE FLUX VALUES... Number of sampled particles:   0.0000000000000000     ]
                # Ignore this line and increase the counter by 1
                changedLines+=1
            elif line_stripped.count(' ') > 9 and len(line_split) == 2 and all(x.replace('.', '', 1).isdigit() for x in line_split):
                # 3. [   3382.0000000000000        1694.4798966228368]
                # Ignore this line and increase the counter by 1
                changedLines+=1
            elif args.iter and line_stripped.count(' ') > 16 and 'iter:' in line and 'time:' in line and len(line_split) > 3 and hasNumbers(line_stripped):
                # OPTIONAL 4. [iter:                  702 time:   3.2151600000000616E-008]
                changedLines+=1
                #print(line_split)
                #exit(0)
            elif line_stripped.count(' ') > 0 and any(substring in line for substring in (' Reason:',' Iterations:',' Norm:')) and hasNumbers(line_stripped):
                #[ Reason:            4]
                #[ Iterations:            1]
                #[ Norm:   0.0000000000000000]
                changedLines+=1
            elif any(substring in line for substring in ('to mpool ucp_requests','UCX  WARN','mpool.c:','ucp_requests')):
                # remove UCX warnings (e.g. on hawk)
                #[[1669126882.059241] [r34c2t5n4:1727877:0]           mpool.c:54   UCX  WARN  object 0x1dce980 {flags:0x20040 recv length 16 host memory} was not returned to mpool ucp_requests]
                changedLines+=1
            elif line_stripped.startswith('myrank') and hasNumbers(line_stripped) and ('MPI_WIN_UNLOCK_ALL' in line_stripped or 'MPI_WIN_FREE' in line_stripped or 'with WIN_SIZE =' in line_stripped):
                # remove [myrank=      0                         Unlocking SampWallState_Shared_Win with MPI_WIN_UNLOCK_ALL()]
                #     or [myrank=      0                    Freeing window SampWallState_Shared_Win with       MPI_WIN_FREE()]
                #     or [myrank=      0     Allocated CNTotalSide2GlobalSide_Shared_Win with WIN_SIZE =              8769432]
                changedLines+=1
            elif 'WARNING: Memory reaching maximum, RAM is at' in line_stripped:
                # Remove [ WARNING: Memory reaching maximum, RAM is at *****%]
                changedLines+=1
            else:
                # Write the line to the new (clean) file
                output_new.write(line)

    # 3. Rename files
    RenameFiles(changedLines, stdfile, stdfile_backup, stdfile_new, args)

    return changedLines

def CleanDoPrintStatusLine(stdfile,args):

    stdfile_new    = stdfile+".new"
    stdfile_backup = stdfile+".bak"

    # 1. Search for carriage-return characters

    changedLines=0
    with open(stdfile, "r") as input:
            with open(stdfile_new, "w") as output_new:
                for line in input:
                    # Remove carriage-return
                    output_new.write(line.rstrip()+"\n")


    with open(stdfile_new) as output_new:
        lines = output_new.readlines()

    #Time = 0.8601E-07    dt = 0.1000E-10   eta =      0:44:32     |=========================================>        | [ 82.51%]
    arr = ['Time', 'dt', 'eta', '%', '|']
    n=0
    with open(stdfile_new, "w") as output_new:
        for line in lines:
            n+=1
            if all(c in line.strip("\n") for c in arr):
                # Ignore this line and increase the counter by 1
                changedLines+=1
            else:
                # Write the line to the new (clean) file
                output_new.write(line)

    # 3. Rename files
    RenameFiles(changedLines, stdfile, stdfile_backup, stdfile_new, args)

    return changedLines

def CleanLostParticles(stdfile,args):

    if not os.path.exists(stdfile):
        print("Error: the file does not exist : %s" % stdfile)
        exit(1)

    stdfile_lost   = stdfile+".lost"
    stdfile_new    = stdfile+".new"
    stdfile_backup = stdfile+".bak"
    h5_output_lost = stdfile+"-lost-particles.h5"

    # 1. Locate the lines corresponding to lost particles (write them to a different file), print the remaining lines to a "new" file
    killList = {} # store the myrank for checking blocks that belong together

    n=0
    nLostParts=0
    meshFound=False
    with open(stdfile, "r") as input:
        with open(stdfile_lost, "w") as output_lost:
            with open(stdfile_new, "w") as output_new:
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
                            # Check for MeshFile, as it is required for the .h5 file that is written when lost particles are found
                            if 'MeshFile' in line.strip("\n"):
                                meshFound=True
                                myline=line.replace(" ", "")
                                tmp=myline.split('|')
                                for x in tmp:
                                    if x.endswith('.h5'):
                                        MeshFile=x

                #print(red("Lost %s particles" % nLostParts))

    # 2. If lost particles have been found, they now are written to a separate .h5 file (PartData) for access via ParaView or VisIT
    killList = {} # store the myrank for checking blocks that belong together

    # Read std.out.lost file with all lost particles and write the data to a .h5 file
    if nLostParts > 0:

        # 2.1. Open .h5 file
        f1 = h5py.File(h5_output_lost,'w')

        data = np.zeros(( nLostParts,12))

        # 2.2. Read lost particles file
        with open(stdfile_lost, "r") as f:
            lines = f.readlines()

        # 2.3. Read lost particles and sort particle data into array
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
                    data[n,0] = xLastPos #x -> store LastPartPos instead of PartPos (the position that was still within the domain)
                    data[n,1] = yLastPos #y -> store LastPartPos instead of PartPos (the position that was still within the domain)
                    data[n,2] = zLastPos #z -> store LastPartPos instead of PartPos (the position that was still within the domain)
                    data[n,3] = vx
                    data[n,4] = vy
                    data[n,5] = vz
                    data[n,6] = SpecID
                    data[n,7] = Element
                    data[n,8] = PartID
                    data[n,9] = x #LastPos -> store new position (outside of the domain)
                    data[n,10] = y #LastPos -> store new position (outside of the domain)
                    data[n,11] = z #LastPos -> store new position (outside of the domain)

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


        # 2.4. Write dummy DG_Solution container
        data1 = np.zeros(( 0,0))
        dset1 = f1.create_dataset('DG_Solution', shape=data1.shape, dtype=np.float64)

        # 2.5. Create new dataset 'dset'
        dset = f1.create_dataset('PartData', shape=data.shape, dtype=np.float64)
        #dset = f1.create_dataset(data, shape=data.shape, dtype=getattr(np, str(dataType)))

        # 2.6. Write as C-continuous array via np.ascontiguousarray()
        if not data.any() :
            print(" %s has dimension %s. Skipping" % (data_set,data.shape))
            pass
        else :
            dset.write_direct(np.ascontiguousarray(data))

        # 2.7. Write attributes
        f1.attrs.modify('File_Version',[1.5])   # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('NodeType',[b'GAUSS'])  # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('File_Type',[b'PartStateBoundary']) # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('Program',[b'PICLas'])  # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('VarNames',[b'empty'])  # these brackets [] are required for ParaView plugin !
        f1.attrs.modify('Time',[0.])            # these brackets [] are required for ParaView plugin !
        f1.attrs.create('MeshFile', [MeshFile], None, dtype='<S255')
        f1.attrs.create('N', [1], None, dtype='i4')
        f1.attrs.create('Project_Name',[b'Lost-particles'], None, dtype='<S255')  # these brackets [] are required for ParaView plugin !

        # old: PartPos as x,y, z for particles
        # f1.attrs.create('VarNamesParticles', [b'ParticlePositionX'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'ParticlePositionY'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'ParticlePositionZ'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'VelocityX'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'VelocityY'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'VelocityZ'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'Species'.ljust(255),    # .ljust(255) is required for ParaView plugin !
        #                                       b'ElementID'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'PartID'.ljust(255),     # .ljust(255) is required for ParaView plugin !
        #                                       b'LastPos-X'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'LastPos-Y'.ljust(255),  # .ljust(255) is required for ParaView plugin !
        #                                       b'LastPos-Z'.ljust(255)], None, dtype='<S255')   # .ljust(255) is required for ParaView plugin !

        # new: LastPartPos as x,y, z for particles (because this position is still in the simulation domain)
        f1.attrs.create('VarNamesParticles', [b'LastPartPos-X'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'LastPartPos-Y'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'LastPartPos-Z'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'VelocityX'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'VelocityY'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'VelocityZ'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'Species'.ljust(255),    # .ljust(255) is required for ParaView plugin !
                                              b'ElementID'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'PartID'.ljust(255),     # .ljust(255) is required for ParaView plugin !
                                              b'ParticlePositionX'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'ParticlePositionY'.ljust(255),  # .ljust(255) is required for ParaView plugin !
                                              b'ParticlePositionZ'.ljust(255)], None, dtype='<S255')   # .ljust(255) is required for ParaView plugin !

        # 2.8. Close .h5 data file
        f1.close()

    # 3. Rename files
    RenameFiles(nLostParts, stdfile, stdfile_backup, stdfile_new, args)
    os.remove(stdfile_lost) # remove lost file (is empty when no particles were lost)

    return nLostParts


def CleanFile(stdfile,args) :

    changedLines = 0

    # 1.  Remove lost particles from std.out file and write them to a .h5 file
    nLostParts = CleanLostParticles(stdfile,args)

    # 2.  Remove remnants of DoPrintStatusLine=T
    changedLines += CleanDoPrintStatusLine(stdfile,args)

    # 3.  Remove lines containing: "REACTION       85084           5           1"
    changedLines += CleanSingleLines(stdfile,args)

    return nLostParts, changedLines


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
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for cleaning std*.out files.\nSupply a single file or a group of files by using the wildcard "*", e.g. std* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Files (std*.out) that are to be cleaned.', nargs='+')
parser.add_argument('-d', '--debug', action='store_true', help='Print additional information regarding the files onto screen.')
parser.add_argument('-s', '--save', action='store_true', help='Save *.bak backup files to see what was actually removed from the std-x.out files.')
parser.add_argument('-i', '--iter', action='store_false', help='Do not remove lines matching "iter:   702 time:   3.2151600000000616E-008" (default=False).')

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

# Ignore:  - cleaned (".new")
#          - lost particles text files (".lost")
#          - state files (".h5")
#          - backup files (".bak")
ext = ['.new', '.lost', '.h5', '.bak']

for stdfile in args.files :
    # Check if file can be skipped, see variable "ext" with all file extensions that are to be ignored
    if stdfile.endswith(tuple(ext)):
        print("%s " % stdfile + yellow("(skipping)"))
        continue

    # Clean-up the current file
    nLostParts, changedLines = CleanFile(stdfile,args)

    # Display results
    if nLostParts > 0:
        if changedLines > 0:
            print("%s " % stdfile + red("Lost %s particles" % nLostParts) + " Written particles to %s-lost-particles.h5" % stdfile + red(" and removed %s lines" % changedLines))
        else:
            print("%s " % stdfile + red("Lost %s particles" % nLostParts) + " Written particles to %s-lost-particles.h5" % stdfile)
    else:
        if changedLines > 0:
            print("%s " % stdfile + red("Removed %s lines" % changedLines) )
        else:
            print("%s" % stdfile)

print(132*"-")
