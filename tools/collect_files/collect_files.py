import numpy as np
from timeit import default_timer as timer
import argparse
import re
import shutil
import os
from collections import OrderedDict

# Bind raw_input to input in Python 2
try:
    input = raw_input
except NameError:
    pass

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
























# Start the timer
start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\n\
Collects files (e.g. ElemTimeStatistics.csv) from a user-supplied folder recursively by sweeping through all sub-directories.\n\
All files found will be copied into the current directory and are renamed using the host directory name as prefix.\n\n\
Example:\n\
  -  files: ElemTimeStatistics.csv\n\
  - source: path/to/TestCases/\n\n\
where TestCases/ contains the sub-directories folder1, folder2, folder3, each containing a separate version of ElemTimeStatistics.csv would result in\n\
  - folder1_ElemTimeStatistics.csv\n\
  - folder2_ElemTimeStatistics.csv\n\
  - folder3_ElemTimeStatistics.csv\n\n\
in the directory where the script is executed.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Name of the files that are to be collected.')
parser.add_argument('source', type=str, help='Directory path where the recursive search starts.')
#parser.add_argument('-d', '--debug', action='store_true', help='Print additional information regarding the files onto screen.')
#parser.add_argument('-c', '--clean', action='store_true', help='Clean-up afterwards by removing any *.bak backup files.')

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
if args.files:
    for stdfile in args.files :
        max_length = max(max_length,len(stdfile))

    InitialDataRead = True
    NbrOfFiles = len(args.files)

#print(132*"-")



paths   = []
home = os.getcwd()
#print(args.source)
print(yellow(os.path.abspath(args.source)))
print( )

for (dirpath, dirnames, filenames) in os.walk(args.source, topdown=True):
    for file in filenames:
        path = os.path.join(dirpath, file)
        path = os.path.abspath(path)
        oldfile = path
        head,tail = os.path.split(path) 
        relpath = os.path.relpath(path,args.source)
        #if tail == 'ElemTimeStatistics.csv':
        if tail == args.files:
            with open(oldfile) as input:
                #print("load %s" % oldfile)
                lines = input.readlines()

                # The file should have at least two data points
                if len(lines) > 2:
                    basename = os.path.basename(head)
                    newfile = relpath.replace(os.path.sep, '_')
                    newfile = "%s/%s" % (home,newfile)
                    #newfile = "%s/%s_%s" % (home,basename, tail)
                    print(blue(basename)+'/'+green(tail))
                    #print(tail)
                    shutil.copyfile(oldfile, newfile)

#print(132*"-")

print('='*132)









