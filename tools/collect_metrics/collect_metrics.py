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

def assignOrder(order):
  #@decorator
  def do_assignment(to_func):
    to_func.order = order
    return to_func
  return do_assignment

class CollectedParameters():

  #don't decorate functions you don't want called
  def __init__(self):
    #don't call this one either!
    self.egg = 2
    #self.data = type('data',(object,),{})()
    #self.data = {}


  @assignOrder(1)
  def GridCells(self,line):
    if '#GridCells' in line :
        content = line.split(':')
        self.data['GridCells'] = np.int_(np.double(content[1]))
        return True
    else :
        return False

  @assignOrder(2)
  def Procs(self,line):
    if '#Procs' in line :
        content = line.split(':')
        self.data['Procs'] = np.int_(np.double(content[1]))
        return True
    else :
        return False

  @assignOrder(3)
  def DOFs(self,line):
    if '#DOFs' in line :
        content = line.split(':')
        self.data['DOFs'] = np.int_(np.double(content[1]))
        return True
    else :
        return False

  @assignOrder(4)
  def Timestep(self,line):
    if 'Initial Timestep' in line :
        content = line.split(':')
        self.data['timestep'] = np.double(content[1])
        return True
    else :
        return False

  @assignOrder(5)
  def DOFs(self,line):
    if 'PICLAS FINISHED!' in line :
        #print(line)
        content = re.search(r"\[([A-Za-z0-9_\.\s]+)\]", line)
        content = content.group(1)
        if 'sec' in content :
            content = content.replace('sec', '')
        self.data['PICLasTime'] = np.double(content)
        #self.data['PICLasTime'] = np.double(re.sub('[^0-9]','', content.group(1)))
        return True
    else :
        return False

  @assignOrder(6)
  def Efficiency(self,line):
    if 'EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]:' in line :
        #print(yellow(line))
        line = line.replace('EFFICIENCY: SIMULATION TIME PER CALCULATION in [s]/[Core-h]:', '')
        if 'sec/h' in line :
            line = line.replace('sec/h', '')
        content = re.search(r"\[([A-Za-z0-9_\.\s\-]+)\]", line)
        content = content.group(1)
        self.data['Efficiency'] = np.double(content)
        return True
    else :
        return False

  @assignOrder(7)
  def Timesteps(self,line):
    if '#Timesteps' in line :
        content = line.split(':')
        self.data['Timesteps'] = np.int_(np.double(content[1]))
        return True
    else :
        return False

  @assignOrder(8)
  def PID(self,line):
    if 'PID: CALCULATION TIME PER TSTEP/DOF:' in line :
        #print(line)
        line = line.replace('PID: CALCULATION TIME PER TSTEP/DOF:', '')
        content = re.search(r"\[([A-Za-z0-9_\.\s\-]+)\]", line)
        content = content.group(1)
        if 'sec' in content :
            content = content.replace('sec', '')
        self.data['PID'] = np.double(content)
        #self.data['PICLasTime'] = np.double(re.sub('[^0-9]','', content.group(1)))
        return True
    else :
        return False



def AnalyzeLine(functions,line):
    for func in functions:
      #line='#GridCells :    3.6320000E+03'
      found = func(line)
      if found:
          break










# Start the timer
start = timer()

"""get command line arguments"""
parser = argparse.ArgumentParser(description='DESCRIPTION:\n\
Collects metrics from std*.out files.\n\
Supply a single file or a group of files by using the wildcard "*", e.g. std* for a list of file names or supply nothing to analyze all std.out files in the current directory.\n\n\
Collected information includes\n\
  - Number of grid cells (#GridCells)\n\
  - Number of processors (#Procs)\n\
  - Number of field degrees of freedom (#DOFs)\n\
  - Time step (Initial Timestep)\n\
  - Total PICLas wall time (PICLAS FINISHED!)\n\
  - Simulation efficiency (SIMULATION TIME PER CALCULATION in [s]/[Core-h])\n\
  - Number of time steps (#Timesteps)\n\
  - Performance index (PID: CALCULATION TIME PER TSTEP/DOF)\n\n\
The output is stored in .stdfile.collected.csv', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-f', '--files', type=str, help='Files (std*.out) that are to be cleaned.', nargs='+')
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

# Ignore:  - cleaned (".new")
#          - lost particles text files (".lost")
#          - state files (".h5")
#          - backup files (".bak")
ext = ['.new', '.lost', '.h5', '.bak']
files = []
if args.files is not None:
    for stdfile in args.files :
        # Check if file can be skipped, see variable "ext" with all file extensions that are to be ignored
        if stdfile.endswith(tuple(ext)):
            print("%s " % stdfile + yellow("(skipping)"))
            continue

        files.append(stdfile)
    # Over-write files list with cleaned list
    args.files = files

print(132*"-")





MyCollParams = CollectedParameters()
functions = sorted(
             #get a list of fields that have the order set
             [
               getattr(MyCollParams, field) for field in dir(MyCollParams)
               if hasattr(getattr(MyCollParams, field), "order")
             ],
             #sort them by their order
             key = (lambda field: field.order)
            )


total = []
if args.files :
    for stdfile in args.files :
        MyCollParams.data = OrderedDict()
        MyCollParams.data['stdfile'] = stdfile
        with open(stdfile) as input:
            print("load %s" % stdfile)
            lines = input.readlines()

        for line in lines:
            AnalyzeLine(functions,line.strip())

        total.append(MyCollParams.data)
else :
    directory = os.getcwd()
    for stdfile in os.listdir(directory):
        if stdfile.startswith("std") and stdfile.endswith(".out"): 
            MyCollParams.data = OrderedDict()
            MyCollParams.data['stdfile'] = stdfile
            with open(stdfile) as input:
                print("load %s" % stdfile)
                lines = input.readlines()

            for line in lines:
                AnalyzeLine(functions,line.strip())

            total.append(MyCollParams.data)
        else:
            continue


n=0
OutputFile = ".stdfile.collected.csv"
with open(OutputFile, "w") as output:
    for item in total:

        if n==0:
            header = ''.join(',"%s"' % key for key, value in item.items() )
            header = header[1:]
            print(header)
            # Write the line to the file
            output.write(header+'\n')
            

        n+=1
        print("")
        for key, value in item.items() :
            print(blue(key),green(str(value)))
        line = ''.join(',%s' % value for key, value in item.items() )
        line = line[1:]
        print(line)
        output.write(line+'\n')
