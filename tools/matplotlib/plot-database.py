#---------------------------------------------------------------------------------------
# import modules
import csv, math, sys
import numpy as np
from matplotlib.pyplot import *
import os
sys.path.append('./')
from format_settings import *
# export to tikz
from matplotlib2tikz import save as tikz_save
#---------------------------------------------------------------------------------------

print (' -------------------------------------------------------------------------------')
print ('     Plot      FFTs                                                             ')
print (' -------------------------------------------------------------------------------')

#---------------------------------------------------------------------------------------
# filenames
filenames = ['./../implicit-test-2/Database.csv']
nfiles    =len(filenames)
legendentries = ['He$^+$','e$^-$','E$_{el}$']
column=[1,2,3,5]
scale =[1e6,1e0,1e0,1e0]
plotstart=[1,3]
nPlots=len(plotstart)
doscaleX= 0
doscaleY= 0
ncolumn=len(column)
dotikz=1
imagename =['sheath-particle-number','sheath-electric-field']
#Variables=['$t$','$E_{pot}$','$E_{kin}$']
Variables=['$t$','$N$','$E$']
Units=['$\mu$s','-','J']
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# define function for reading data from cvs file
#---------------------------------------------------------------------------------------
def readingcsv( filename, ncolumn):
    print (' File: ', filename)
    with open(filename,'r') as csvfile:
        reader = csv.reader(csvfile)# , delimiter=',')
        row_count = sum(1 for row in reader)
        print( ' Number of rows:    ', row_count)
        n=int(np.log(row_count)/np.log(2))
        data=np.zeros((row_count-1,ncolumn)) # default is numpy.float64
        ii=-1
        jj=0
        csvfile.seek(0)
        for row in reader:
            if(ii>=0):
                for iCol in column:
                    ind=column.index(iCol)
                    data[ii,ind] = float(row[iCol-1])*float(scale[ind])
            ii+=1
    return data

#---------------------------------------------------------------------------------------
# generate graph
#---------------------------------------------------------------------------------------
def plotdata( iPlot,plotstart,data):
    # figure with energy
    figure(1,(fig_width,fig_height))
    if(iPlot<len(plotstart)-1):
        nGraphs=plotstart[iPlot+1]-plotstart[iPlot]
    else:
        nGraphs=1
    if(doscaleX==1):
        xName= ''.join([Variables[0]," $\cdot$ ", str('{:2.4}'.format(scale[0])), " / " ,Units[0]])
    else:
        xName= ''.join([Variables[0], " / " ,Units[0]])
    if(doscaleY==1):
        yName= ''.join([Variables[iPlot+1]," $\cdot$ ", str('{:2.4}'.format(scale[1])), " / " ,Units[iPlot+1]])
    else:
        yName= ''.join([Variables[iPlot+1], " / " ,Units[iPlot+1]])
    xlabel(xName)
    ylabel(yName)
    # grid and gridlines
    grid(True)
    xgridlines = getp(gca(), 'xgridlines')
    ygridlines = getp(gca(), 'ygridlines')
    setp(xgridlines, 'linestyle', '-')
    setp(ygridlines, 'linestyle', '-')
    # position of orrigin 
    ax = gca()
    #ax.ticklabel_format(style='sci', scilimits=(1,0), axis='y')
    #ax.ticklabel_format(style='sci', scilimits=(1,0), axis='x')
    matplotlib.rc('lines', linewidth=0.5, markersize=5)
    # plot data via loop
    PlotGraph0 = plotstart[iPlot]
    for dataset in data:
        for iGraph in range(0,nGraphs):
            PlotGraph = PlotGraph0+iGraph
            plot(dataset[:,0],dataset[:,PlotGraph],label=legendentries[PlotGraph-1])

    # grid locations 
    #majoryLocator   = MultipleLocator(0.02)
    #ax.yaxis.set_major_locator(majoryLocator)
    #ax.set_xlabel("$t$ [ns]")
    #ax.xaxis.set_label_coords(0.5, -0.05)
    #ax.yaxis.set_label_coords(-0.08, 0.5)
    ##ax.set_yscale('log')
    #ax.set_xlim([10,60])
    #ax.set_ylim([0,20])
    # legend 
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=int(2),prop = font)
    savefig(''.join([imagename[iPlot],'.pdf'])) 
    if(dotikz==1):
        tikz_save(''.join([imagename[iPlot],'.tikz']),
                  figureheight = '\\figureheight',
                  figurewidth = '\\figurewidth'
                  )
    show()
    #close(1)
    return 

#---------------------------------------------------------------------------------------
# HERE: script starts
#---------------------------------------------------------------------------------------
print (' Reading ....')
print (' number of files: ', nfiles)

data=[]
for filename in filenames:
    data_in= readingcsv(filename,ncolumn)
    data.append(data_in)

for iPlot in range(nPlots):
    # plot data
    plotdata( iPlot,plotstart,data)
