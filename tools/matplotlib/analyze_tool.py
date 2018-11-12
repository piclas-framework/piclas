#!/usr/bin/python
import os, sys
import numpy as np
import re
import csv

#----------------------------------------------------------------------------------------------------------------------
# classes
#----------------------------------------------------------------------------------------------------------------------
foldername='CFL' #raw_input(' Folder-name for range: ')
lowerLimit=1 #int(raw_input(' Lower-Limit: '))
upperLimit=401 #int(raw_input(' Upper-Limit: '))+1
inc=1 #int(raw_input(' Increment: '))
outputfile='results_explicitO4-p03.csv'

# pattern in std.out for extraction
pattern=[]
pattern.append('CFLScale')                        # number of iterations
pattern.append('iter:')                            # number of iterations
pattern.append('PICLAS FINISHED')              # runtime
pattern.append('Total iteration Linear Solver')    # iteration lin-solver
pattern.append('Total iteration outer-Newton')     # outer Newton
scaleFloat=[1.,1.,0.00379665135350620752,1.,1.,]

# information from csv file
doAnalyzeCSV=1
# definition for sinus comparison
Amplitude = 2.3178e-5
frequency = 5.3591e7
tDelay = 2e-7
column=[1,5]

# --------------------------------------------------------------------------
# get float from string
# --------------------------------------------------------------------------
numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

#---------------------------------------------------------------------------------------
# define analytic error function (here sinus)
# input: data
# output: local quadratic error
#---------------------------------------------------------------------------------------
def errorfunction( data):
    value=Amplitude*np.cos(2.*np.pi*frequency*data[0])**2
    # compute quatratic error
    error=(data[1]-value)**2
    return error

#---------------------------------------------------------------------------------------
# read in csv and compute error norm
# input:  filename
# output: error-norm L2
#---------------------------------------------------------------------------------------
def readingcsv( filename ):
    with open(filename,'r') as csvfile:
        reader = csv.reader(csvfile)# , delimiter=',')
        row_count = sum(1 for row in reader)
        #n=int(np.log(row_count)/np.log(2))
        data=np.zeros(2) # default is numpy.float64
        ii=-1
        jj=0
        csvfile.seek(0)
        time=0
        error=0
        n=0
        for row in reader:
            if(ii==row_count-1):
                break 
            if(ii>=0):
                for iCol in column:
                    ind=column.index(iCol)
                    data[ind] = float(row[iCol-1])
                if data[0]>=tDelay:
                    # compute L2 error
                    data[0]=data[0]-tDelay
                    n=n+1
                    error=error+errorfunction(data)
            ii+=1
        errornorm=np.sqrt(error/n)
    return errornorm

# --------------------------------------------------------------------------
# get size of output array
# --------------------------------------------------------------------------
# number of processed dirs
ndirs=upperLimit-lowerLimit
# initialize data array
lenpat=len(pattern)
if(doAnalyzeCSV>0):
    data=np.zeros((ndirs,lenpat+1))
    pattern.append('L2')
else:
    data=np.zeros((ndirs,lenpat))

iFile=0
# --------------------------------------------------------------------------
# loop folder structure 
# --------------------------------------------------------------------------
for ii in range(lowerLimit,upperLimit,inc):
    dirname= ''.join([foldername, str(ii)]) # inserting a list
    # get output-file
    nout=0
    for file in os.listdir(dirname):
        if file.startswith('std'):
            outfile=file
            nout=nout+1
    if(nout==0):
        sys.exit(' Exit: No std-XXX.out found.')
    elif(nout>1):
        sys.exit(' Exit: Too many std-XXXX.out files')
    #else:
    #     print ' Found std-XXX.out: ', outfile
    # --------------------------------------------------------------------------
    # and get date from output file
    # --------------------------------------------------------------------------
    filename=''.join([dirname,'/',outfile])
    input_file = open(filename,'r')
    for line in input_file:
        for ii in range(0,lenpat):
            if pattern[ii] in line.strip():
                value =rx.findall(line)
                data[iFile,ii] = scaleFloat[ii]*float(value[0])
    # --------------------------------------------------------------------------
    # and get date from output file
    # --------------------------------------------------------------------------
    if(doAnalyzeCSV>0):
        filename=''.join([dirname,'/Database.csv'])
        data[iFile,ii+1]= readingcsv(filename)
    # count files
    iFile=iFile+1

#----------------------------------------------------------------------------------------------------------------------
# output csv
#----------------------------------------------------------------------------------------------------------------------

with open(outputfile,'wt') as f:
    writer = csv.writer(f)
    writer.writerow(pattern)
    writer.writerows(data)
