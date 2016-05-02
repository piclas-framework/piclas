#---------------------------------------------------------------------------------------
# import modules
import csv, math, sys
import numpy as np
import matplotlib.pyplot as plt
import os
import h5py
sys.path.append('./')
#from format_settings import *
#---------------------------------------------------------------------------------------

print (' -------------------------------------------------------------------------------')
print ('     Analyzing  Load ImBalance                                                  ')
print (' -------------------------------------------------------------------------------')

#---------------------------------------------------------------------------------------
# filenames
#filename = 'PlasmaPlume_State_000.000000000007500.h5'
filename = 'PlasmaPlume_State_000.000000000000510.h5'
outputfile='PlasmaPlume-Elemload.csv'
outputtocsv = 0
doloaddistri= 1
PartMPIWeight=0.02
nProcs      = 9216
#---------------------------------------------------------------------------------------

print (' Reading ....')

file=h5py.File(filename,'r')

PartInt = file['PartInt']
nElems  = PartInt.shape[1]

print (' nElems:        ', nElems )

PartsInElem=np.zeros(nElems)

PartsInElem = PartInt[1,:] - PartInt[0,:]
print(' ElemWeight ....')
node = 'ElemWeight'
if node in file.keys():
    print(' Reading ElemWeight...')
    ElemWeight = file['ElemWeight']
else:
    print(' Build ElemWeight out of PartInt...')
    for iElem in range(0,nElems):
        ElemWeight[iElem]=1.0+float(PartInt[1,iElem]-PartInt[0,iElem])*PartMPIWeight
print(' ...DONE')

ElemID=np.zeros(nElems)
for row in range(0,nElems):
    ElemID[row]=row+1
#----------------------------------------------------------------------------------------------------------------------
# output of data
#----------------------------------------------------------------------------------------------------------------------
if(outputtocsv==1):
    print (' Writing ElemWeight to csv ....')
    csv_file = open(outputfile,'w')
    csv_writer= csv.writer(csv_file, delimiter=',')
    csv_writer.writerow([str('Element,ElementLoad,PartsInElem')])
    
    for row in range(0,nElems):
        csv_writer.writerow([row+1,ElemWeight[row],PartsInElem])
    csv_file.close()
    print (' Done.')


#----------------------------------------------------------------------------------------------------------------------
# compute load distribution
#----------------------------------------------------------------------------------------------------------------------
if(doloaddistri==1):
    print(' Performing load distribution step: ')
    SumWeight=np.sum(ElemWeight)
    MeanWeight=SumWeight/float(nProcs)
    
    print(' Total   weight: ', SumWeight )
    print(' Average weight: ', MeanWeight )

    finished = "false"
    targetWeight=MeanWeight#*(iProc+1)
    LastProcDiff=0.
    nTrys=0

    while finished != "true":
        targetWeight=targetWeight+(LastProcDiff)/float(nProcs)
        print(' new target weight ', targetWeight)
        curiElem=1
        offSetElem=np.zeros(nProcs+1)
        LoadDistri=np.zeros(nProcs)
        LoadDiff=np.zeros(nProcs)
        oldCurWeight=0.
        for iProc in range(0,nProcs):
            offSetElem[iProc] = curiElem-1
            if(iProc>0 and offSetElem[iProc]==offSetElem[iProc-1]):
                print('blubb')
            #print('iproc',iProc,offSetElem[iProc])
            CurWeight=0.
            getElem=0
            for iElem in range(curiElem-1, nElems-nProcs + 1+ iProc):
                getElem=getElem+1
                CurWeight=CurWeight+ElemWeight[iElem]
                #if(iProc==nProcs-2):
                #    print('iElem',iElem)
                if(CurWeight>=targetWeight) or (iElem == nElems-nProcs+iProc):
                    diffLower=(CurWeight-ElemWeight[iElem]-targetWeight)
                    diffCurr=(CurWeight-targetWeight)
                    #print (' lower, curr', diffLower,diffCurr)
                    if(getElem>1):
                        if(iProc==nProcs-1):
                                LoadDiff[iProc]=diffCurr
                                curiElem=iElem+2
                                LoadDistri[iProc]=CurWeight
                                break
                        else:
                            if(abs(diffLower)<abs(diffCurr) and iElem <nElems-nProcs+iProc):
                                LoadDiff[iProc]=diffLower
                                # caution for last elem
                                curiElem=iElem+1
                                LoadDistri[iProc]=CurWeight-ElemWeight[iElem]
                                break
                            else:
                                LoadDiff[iProc]=diffCurr
                                curiElem=iElem+2
                                LoadDistri[iProc]=CurWeight
                                break
                    else:
                        LoadDiff[iProc]=diffCurr
                        curiElem=iElem+2
                        LoadDistri[iProc]=CurWeight
                        break
        
        ElemDistri=np.zeros(nProcs)
        offSetElem[nProcs]=nElems
        for iProc in range(0,nProcs):
            ElemDistri[iProc]=offSetElem[iProc+1]-offSetElem[iProc]
        ElemDistri[nProcs-1]=nElems-offSetElem[nProcs-1]
        LoadDistri[nProcs-1]=np.sum(ElemWeight[int(offSetElem[nProcs-1]):nElems])
        LoadDiff[nProcs-1]=(LoadDistri[nProcs-1]-targetWeight)
        
        #print(' first distri')
        #print('      Elem Distri: ', ElemDistri)
        #print('      OffsetElem : ', offSetElem)
        #print('      Load Distri: ', LoadDistri)
        #print('      Diff Load  : ', LoadDiff)
        print('      Sum  Load : ', np.sum(LoadDistri))
    
        MaxLoadDiff=np.max(LoadDiff[0:nProcs-1])
        MinLoadDiff=np.min(LoadDiff[0:nProcs-1])
        LastLoadDiff=LoadDiff[nProcs-1]
        print(' MaxLoadDiff: ', MaxLoadDiff)
        print(' MinLoadDiff: ', MinLoadDiff)
        print(' LastLoadDiff:', LastLoadDiff)
        LastProcDiff=LastLoadDiff-MaxLoadDiff
        print(' LastProcDiff', LastProcDiff)
        if(LastProcDiff>0):
            finished="false"
            nTrys=nTrys+1
        else:
            finished="true"
        if(nTrys==10):
            print('no good load distribution found!')
            break
    #-----------------
    # create plot
    #-----------------
    RMS=0
    for iProc in range(0,nProcs):
        RMS=RMS+(LoadDistri[iProc]-targetWeight)**2
    RMS = RMS/float(nProcs-1)
    Mean=np.sum(LoadDistri)/float(nProcs)
    print( ' Mean   : ', Mean)
    print( ' sigma  : ', RMS)
        
    ProcID=np.zeros(nProcs)
    for row in range(0,nProcs):
        ProcID[row]=row
    fig, ax1= plt.subplots()
    #x1.set_xlim(1, nProcs)
    ax1.set_ylim(0, 1.25)
    ax2 = ax1.twinx()
    maxDiff=np.max(np.abs(LoadDiff))
    iDiff=int(maxDiff/RMS)+1
    #iDiff=int(maxDiff/Mean)+1
    if iDiff<5:
        iDiff=5
    ax2.set_ylim(-iDiff, iDiff)
    ax2.set_xlim(1, nProcs)
    #p1  = ax1.plot(ProcID,LoadDistri/Mean ,linestyle='none', marker='o', color='g', label = 'distribution')
    p1  = ax1.plot(ProcID,LoadDistri/Mean ,linestyle='-',color='g', label = 'distribution')
    #p1  = ax1.plot(ProcID,LoadDistri, 'g-', label = 'distribution')
    #p2  = ax2.plot(ProcID,LoadDiff/RMS ,linestyle='none', marker='o', color='b', label = 'imbalance')
    p2  = ax2.plot(ProcID,LoadDiff/RMS ,linestyle='-', color='b', label = 'imbalance')
    ptot =p1+p2
    legs=[l.get_label() for l in ptot]
    ax1.legend(ptot, legs, loc=0)
    ax1.set_xlabel('$n_{Processor}$ / -')
    ax1.set_ylabel('load/load$_{mean}$')
    ax2.set_ylabel('load diff / $\sigma$ / -')
    #ax2.set_ylabel('load diff / load$_{mean}$ / -')
        


#----------------------------------------------------------------------------------------------------------------------
# create plots
#----------------------------------------------------------------------------------------------------------------------

print(" Create plots")

fig, ax1 = plt.subplots()
ax1.set_xlim(1, nElems)

#ax2 = ax1.twinx()
#ax2.set_xlim(1, nElems)
#p1  = ax1.plot(ElemID,ElemWeight,linestyle='none', marker='o', color='g', label = 'time')
p1  = ax1.plot(ElemID,ElemWeight,linestyle='-', color='g', label = 'time')
#p2  = ax2.plot(ElemID,PartsInElem, 'b-', label = 'nParts')
 
ptot =p1#+p2
legs=[l.get_label() for l in ptot]
ax1.legend(ptot, legs, loc=0)

ax1.set_xlabel('elements / -')
ax1.set_ylabel('time / s')
#ax2.set_ylabel('parts in element / -')

fig, ax1 = plt.subplots()
ax1.set_xlim(1, nElems)

#ax2 = ax1.twinx()
#ax2.set_xlim(1, nElems)
#p1  = ax1.plot(ElemID,PartsInElem,linestyle='none', marker='o', color='b', label = 'nParts')
p1  = ax1.plot(ElemID,PartsInElem, color='b', label = 'nParts')
 
ptot =p1#+p2
legs=[l.get_label() for l in ptot]
ax1.legend(ptot, legs, loc=0)

ax1.set_xlabel('elements / -')
ax1.set_ylabel('parts in element / -')


plt.show()


