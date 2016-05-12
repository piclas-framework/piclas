#---------------------------------------------------------------------------------------
# import modules
import csv, math, sys
import numpy as np
from matplotlib.pyplot import *
import os
sys.path.append('./')
from format_settings import *
#---------------------------------------------------------------------------------------

print (' -------------------------------------------------------------------------------')
print ('     Computing FFTs                                                             ')
print (' -------------------------------------------------------------------------------')

#---------------------------------------------------------------------------------------
# filenames
filename1 = './../imex2/Database.csv'
filename2 = './../implicit/Database.csv'
column=[1,5,7]
scale =[1e9,1e6,1e6]
imagename = 'comparison_plasmawave-ekin.pdf'
Variables=['$t$','$E$']
Units=['s','J']
#---------------------------------------------------------------------------------------

print (' Reading ....')

print (' File1: ', filename1)
with open(filename1,'r') as csvfile:
    reader = csv.reader(csvfile)# , delimiter=',')
    row_count = sum(1 for row in reader)
    print( ' Number of rows: ', row_count)
    n=int(np.log(row_count)/np.log(2))
    data1=np.zeros((row_count-1,3)) # default is numpy.float64
    ii=-1
    jj=0
    csvfile.seek(0)
    for row in reader:
        if(ii>=0):
            for iCol in column:
                ind=column.index(iCol)
                data1[ii,ind] = float(row[iCol-1])*float(scale[ind])
        ii+=1


print (' File2: ', filename2)
with open(filename2,'r') as csvfile:
    reader = csv.reader(csvfile)# , delimiter=',')
    row_count = sum(1 for row in reader)
    print( ' Number of rows: ', row_count)
    n=int(np.log(row_count)/np.log(2))
    data2=np.zeros((row_count-1,3)) # default is numpy.float64
    ii=-1
    csvfile.seek(0)
    for row in reader:
        if(ii>=0):
            for iCol in column:
                ind=column.index(iCol)
                data2[ii,ind] = float(row[iCol-1])*float(scale[ind])
        ii+=1


# figure with energy
figure(1,(fig_width,fig_height))
xName= ''.join([Variables[0]," $\cdot$ ", str('{:2.4}'.format(scale[0])), " \ " ,Units[0]])
yName= ''.join([Variables[1]," $\cdot$ ", str('{:2.4}'.format(scale[1])), " \ " ,Units[1]])
xlabel(xName)
ylabel(yName)

grid(True)
xgridlines = getp(gca(), 'xgridlines')
ygridlines = getp(gca(), 'ygridlines')
setp(xgridlines, 'linestyle', '-')
setp(ygridlines, 'linestyle', '-')

ax = gca()
#ax.ticklabel_format(style='sci', scilimits=(1,0), axis='y')
#ax.ticklabel_format(style='sci', scilimits=(1,0), axis='x')

matplotlib.rc('lines', linewidth=0.5, markersize=5)
#plot(data1[:,0],data1[:,1],label='E_{el,1}', color='r', marker='s', markevery=360, markerfacecolor='w')
#plot(data2[:,0],data2[:,1],label='E_{el,2}', color='g', marker='^', markevery=440, markerfacecolor='w')
plot(data1[:,0],data1[:,1],label='$E_{el,imex}$', color='b')
plot(data2[:,0],data2[:,1],label='$E_{el,implicit}$', color='r', marker='^', markevery=10, markerfacecolor='w')
plot(data1[:,0],data1[:,2],label='$E_{kin,imex}$', color='g')
plot(data2[:,0],data2[:,2],label='$E_{kin,implicit}$', color='k', marker='^', markevery=10, markerfacecolor='w')

#majoryLocator   = MultipleLocator(0.02)
#ax.yaxis.set_major_locator(majoryLocator)
#ax.set_xlabel("$t$ [ns]")
#ax.xaxis.set_label_coords(0.5, -0.05)
#ax.yaxis.set_label_coords(-0.08, 0.5)
##ax.set_yscale('log')
#ax.set_xlim([10,60])
#ax.set_ylim([0,20])

ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=3,prop = font)
savefig(imagename) 
show()
#close(1)
