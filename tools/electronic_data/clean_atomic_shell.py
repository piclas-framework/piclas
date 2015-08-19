#---------------------------------------------------------------------------------------
# This script is written by Philip Ortwein Anno Domini 2015-08-19
#
#!/usr/bin/env python
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# import modules
import csv, math, sys
import commands
import numpy as np
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
input_filename  = 'ArIon2.csv'
output_filename = 'Database_reduced.csv'
#---------------------------------------------------------------------------------------

def BubbleSort( gin, tin):
    # modified from rosetta code
    changed = True
    while changed:
        changed = False
        for ii in xrange(len(tin) -1):
            if tin[ii] > tin[ii+1]:
                tin[ii], tin[ii+1] = tin[ii+1], tin[ii]
                gin[ii], gin[ii+1] = gin[ii+1], gin[ii]
                changed = True
    return None

def DataInfo( data):
    Data=np.array(data)
    deltaMax=float(0.0)
    deltaMax2=float(0.0)
    deltaMin=float(1e99)
    deltaMean=float(0.0)
    for ii in xrange(len(Data)-2):
        delta= float(Data[ii+1]-Data[ii])
        if ii>1:
            deltaMax2=max(delta,deltaMax2)
        deltaMax=max(delta,deltaMax)
        deltaMin=min(delta,deltaMin)
        deltaMean=deltaMean+delta
    deltaMean2=(deltaMean-deltaMax)/(len(Data)-1)
    deltaMean=deltaMean/(len(Data)-1)
    print ' delta_min:    ', deltaMin#, '\n'
    print ' delta_max:    ', deltaMax#, '\n'
    print ' delta_max/1:  ', deltaMax2#, '\n'
    print ' delta_mean:   ', deltaMean#, '\n'
    print ' delta_mean/1: ', deltaMean2#, '\n'

#def ReduceData( delta,gin, Tin, gout, Tout):
#    gIn=np.array(gin)
#    TIn=np.array(Tin)
#    changed= True
#    iMin=0
#    iMax=0
#    while changed:
#        for ii in range(iMin,len(Tin)-1):
#            deltaLoc= float(TIn[ii+1])-float(TIn[ii])
#            if deltaLoc<delta:
#                iMax=ii
#            else:
#                break
#        iMax=iMax+1
#        gnew=float(0.0)
#        Tnew=float(0.0)
#        for ii in range(iMin,iMax):
#            gnew=gnew+gIn[ii]
#            Tnew=Tnew+TIn[ii]
#            #print ii, gnew, Tnew
#        if iMax-iMin>0:
#            Tnew=Tnew/(iMax-iMin)
#        gout.append( [ float(gnew)] )
#        Tout.append( [ float(Tnew)] )
#        iMin=iMax
#        if iMax==len(TIn)-1:
#            changed= False

def ReduceData( delta,gin, Tin, gout, Tout):
    delta0=delta
    gIn=np.array(gin)
    TIn=np.array(Tin)
    changed= True
    iMin=0
    iMax=0
    while changed:
        for ii in range(iMin,len(Tin)-1):
            deltaLoc= float(TIn[ii+1])-float(TIn[ii])
            if deltaLoc<delta:
                iMax=ii
            else:
                break
        iMax=iMax+1
        gnew=float(0.0)
        Tnew=float(0.0)
        for ii in range(iMin,iMax):
            gnew=gnew+gIn[ii]
            Tnew=Tnew+TIn[ii]
            #print ii, gnew, Tnew
        if iMax-iMin>0:
            Tnew=Tnew/(iMax-iMin)
        gout.append( [ float(gnew)] )
        Tout.append( [ float(Tnew)] )
        iMin=iMax
        delta=delta0*float(1.0)/(float(iMin)**0.4)
        if iMax==len(TIn)-1:
            changed= False




print '\n Reading from file: ', input_filename, '\n'
input_file   = open(input_filename, 'r')
input_reader = csv.reader(input_file, delimiter=',')

gIn = []
TIn = []
for row in input_reader:
  gIn.append( [ float(row[0])] )
  TIn.append( [ float(row[1])] )
input_file.close()

print ' Found number of entries: ', len(TIn), '\n'

print ' ... Readin done.'

print ' Sorting entries...'

BubbleSort( gIn,TIn)

print ' ... Done.','\n'

print ' DataSet Information:'
DataInfo(TIn)
delta = float(raw_input(" Please enter delta for mergin: "))
print " You entered", delta

print ' Reducing Dataset...'
gOut=[]
TOut=[]
ReduceData( delta,gIn, TIn, gOut, TOut)
print ' ... Done.','\n'

print ' Writing Dataset...'

output_file = open(output_filename, 'w')
print ' Reduced Datasetlength:', len(TOut)
for ii in range(0,len(TOut)):
    # debug
    #output_file.write( " ".join(map(str,[ii]))+" ".join(" ")+" ".join(map(str,gOut[ii]))+" ".join(" ")+" ".join(map(str,TOut[ii]))  + "\n")
    # use
    output_file.write( " ".join(map(str,gOut[ii]))+" ".join(" ")+" ".join(map(str,TOut[ii]))  + "\n")
output_file.close()

print ' ... Done.'
