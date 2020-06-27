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
parser = argparse.ArgumentParser(description='DESCRIPTION:\nTool for cleaning std*.out files.\nSupply a single file or a group of files by using the wildcard "*", e.g. std* for a list of file names.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('files', type=str, help='Files (std*.out) that are to be cleaned.')
parser.add_argument('source', type=str, help='Files (std*.out) that are to be cleaned.')
parser.add_argument('-d', '--debug', action='store_true', help='Print additional information regarding the files onto screen.')
parser.add_argument('-c', '--clean', action='store_true', help='Clean-up afterwards by removing any *.bak backup files.')

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
        #if tail == 'ElemTimeStatistics.csv':
        if tail == args.files:

            with open(oldfile) as input:
                #print("load %s" % oldfile)
                lines = input.readlines()
                #print(len(lines))
                if len(lines) > 1:
                    basename = os.path.basename(head)
                    newfile = "%s/%s_%s" % (home,basename, tail)
                    #print(newfile)
                    print(blue(basename)+'/'+green(tail))
                    #print(tail)
                    shutil.copyfile(oldfile, newfile)

#print(132*"-")

print('='*132)









