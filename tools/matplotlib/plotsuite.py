#---------------------------------------------------------------------------------------
# Matplotlib plotting script
# 
# Tools to plot PICLas csv files with matplotlib into tikz format files, pdfs, etc.
# * plot several lines or multiple files
# * output to pdf format
# * output to tikz format
# * simple reduction of data size for tikz format
# * trapezoidal rule for time-average for frequencies -> periods
# * second-order first time-derivative go compute temporal derivatives from given data
#
# Maintainer:
# Philip Ortwein
# ortwein@iag.uni-stuttgart.de
# 
# format_settings:
# Maintainer:
# J. Neudorfer
# A. Stock
#---------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------
# import modules
#---------------------------------------------------------------------------------------
import csv, math, sys
import numpy as np
from matplotlib.pyplot import *
import os
sys.path.append('./')
from format_settings import *
# export to tikz
from matplotlib2tikz import save as tikz_save
from itertools import cycle
#---------------------------------------------------------------------------------------

print (' -------------------------------------------------------------------------------')
print ('     PICLas dataplot-tool                                                       ')
print (' -------------------------------------------------------------------------------')

#---------------------------------------------------------------------------------------
# filenames
#---------------------------------------------------------------------------------------
filenames = ['./results_explicitO4.csv','./results_cn.csv','./results_implicitO4.csv','./results_implicitO4-eisenstat2.csv']
nfiles    =len(filenames)
legendentries = ['LSERKO4','ERKO4','IRKO4','IRKO4-Eisen']
#column=[1,9]
column=[1,3,4,5]
scale =[1.,1.,1.,1.]
nYAxis=len(column)-1
defineLimits=0
axisLimit=[1,72,0,800,5000,100000,0,5000]
plotstart=[1]
nPlots=len(plotstart)
doscaleX= 0
doscaleY= 0
ncolumn=len(column)
plotMarker=0
dotikz=0
Frequency=-1      # 35e9 ! zero nothing to do
TimeDerivative=0 # 0/1 no/yes
TimeMode=0 # mode=0: input value=abs value, e.g. Ekin, mode=1: input=diff value
StartTime=2e-9
imagename =['runtime']
Variables=['CFL','t', 'iterLS', 'iter Newton' ]
Units=[' -',' s', '-', '-']
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
        #n=int(np.log(row_count)/np.log(2))
        data=np.zeros((row_count-1,ncolumn)) # default is numpy.float64
        ii=-1
        jj=0
        csvfile.seek(0)
        for row in reader:
            if(ii==row_count-1):
                break 
            if(ii>=0):
                for iCol in column:
                    ind=column.index(iCol)
                    data[ii,ind] = float(row[iCol-1])*float(scale[ind])
            ii+=1
    return data


#---------------------------------------------------------------------------------------
# generate integral average
# 
# time average is computed by trapezoidal rule, modified to PICLas
# the corresponding time point is set in between the interval for the trapezoidal rule
# integral is not perfectly used, hence a  small error is introduced + small
# frequency shift
#---------------------------------------------------------------------------------------
def averagedata( Frequency,tScale,time,locdata):
    if(Frequency>0):
        print(' average: yes ' )
        tPeriode=1./Frequency
        # rescale period to print time
        time=time/tScale
        # average data
        # get number of periods
        n=size(time)
        iSample=0
        tSample=0.
        for ii in range(1,n):
            dt=time[ii]-time[ii-1]
            # ignore double values of dt
            if(dt>0):
                tSample=tSample+dt
                if(tSample>=tPeriode):
                    iSample=iSample+1
                    tSample=0.
        nPeriods=iSample 
        print(' tPeriode, last-dt:  ',tPeriode,dt)
        print(' nPeriods,nSamples:  ',nPeriods,n)
        dataPlot=np.zeros((nPeriods,2))
        # simple, integration
        tSample=0.
        tStep=time[0]
        averageValue=0.
        dtOld=0.
        dtAvg=0.
        iSample=0
        SamplingFrequenc=0.
        for ii in range(1,n):
            dt=time[ii]-time[ii-1]
            # ignore double values of dt
            if(dt>0):
                tSample=tSample+dt
                dtStep=(dtOld+dt)*0.5
                if(tSample>=tPeriode):
                    dtStep=0.5*dt
                dtOld=dt
                dtAvg=dtAvg+dtStep 
                averageValue=averageValue+locdata[ii]*dtStep
                if(tSample>=tPeriode):
                    dataPlot[iSample,0]=(tStep+0.5*dtAvg)*tScale
                    dataPlot[iSample,1]=averageValue/dtAvg
                    SamplingFrequenc=SamplingFrequenc+dtAvg
                    tSample=0.
                    tStep=time[ii]
                    averageValue=0. 
                    dtOld=0.
                    dtAvg=0.
                    iSample=iSample+1
        print(' check ', iSample)
        print(' mean sampling period: ', SamplingFrequenc/iSample)
    else:
        # copy data
        print(' average: no ' )
        if(Frequency==-1):
            print(' mean value: ', np.mean(locdata[:]))
        n=size(time)
        dataPlot=np.zeros((n,2))
        for ii in range(n):
            dataPlot[ii,0]=time[ii]
            dataPlot[ii,1]=locdata[ii]
    return dataPlot


#---------------------------------------------------------------------------------------
# time derivative 
# 
# the time-derivative is computed by a central differnces between t^n and t^n+1 
# for the new time level t^n+1/2
# TimeMode=0: input abs value like E_kin
# TimeMode=1: input is difference like E_kin_in, which is the kinetic energy inflow over
#             a certain time-interval
#---------------------------------------------------------------------------------------
def timederivative( TimeDerivative,TimeMode,StartTime,tScale,time,locdata):
    if(TimeDerivative>0):
        print(' time-derivative: yes ' )
        # rescale period to print time
        time=time/tScale
        # average data
        # get number of periods
        n=size(time)
        iSample=0
        for ii in range(0,n-1):
            if(time[ii]>=StartTime):
                dt=time[ii+1]-time[ii]
                # ignore double values of dt
                # second condition removes strong neg. peaks
                if(dt>0):
                    iSample=iSample+1
        nPeriods=iSample 
        print(' nSamples of n', nPeriods/n)
        dataPlot=np.zeros((nPeriods,2))
        # second order time derivative of first order term
        iSample=0
        mean=0.
        for ii in range(0,n-1):
            if(time[ii]>=StartTime):
                dt=time[ii+1]-time[ii]
                # ignore double values of dt
                # second condition removes strong neg. peaks
                if(dt>0):
                    tStep=0.5*(time[ii+1]+time[ii])
                    dataPlot[iSample,0] = tStep*tScale
                    if(TimeMode==0):
                        tDeri = (locdata[ii+1]-locdata[ii])/dt
                    else:
                        tDeri = locdata[ii+1]/dt
                    mean=mean+tDeri
                    dataPlot[iSample,1] = tDeri
                    iSample=iSample+1
        print(' check ', iSample/n)
    else:
        # copy data
        print(' time-derivative: no ' )
        n=size(time)
        dataPlot=np.zeros((n,2))
        for ii in range(n):
            dataPlot[ii,0]=time[ii]
            dataPlot[ii,1]=locdata[ii]
    return dataPlot


#---------------------------------------------------------------------------------------
# generate graph
#---------------------------------------------------------------------------------------
def plotdata( iPlot,plotstart,data):
    # figure with energy
    if(nYAxis==0):
        figure(1,(fig_width,fig_height))
        if(iPlot<len(plotstart)-1):
            if(plotstart[iPlot+1]-plotstart[iPlot]>1):
                condense=plotstart[iPlot+1]-plotstart[iPlot]-1
                if(condense>1):
                    sys.exit(" Script cannot combine more than two columns!")
            else:
                condense=0
        else:
            condense=0
        #print(' condense', condense)
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
        if(plotMarker>0):
            matplotlib.rc('lines', linewidth=0.5, markersize=5)
            # plot data via loop
            if(condense>0):
                ii=0
                for dataset in data:
                    plot(dataset[:,0],[np.sum(i) for i in zip(dataset[:,plotstart[iPlot]],dataset[:,plotstart[iPlot+condense]])],label=legendentries[ii])
                    ii+=1
            else:
                ii=0
                for dataset in data:
                    dataPlot0 = timederivative( TimeDerivative,TimeMode,StartTime,scale[0],dataset[:,0],dataset[:,plotstart[iPlot]])
                    dataPlot  = averagedata(Frequency,scale[0],dataPlot0[:,0],dataPlot0[:,1])
                    plot(dataPlot[:,0],dataPlot[:,1],label=legendentries[ii],marker='o')
                    ii+=1
        else:
            matplotlib.rc('lines', linewidth=1.0)
            # plot data via loop
            if(condense>0):
                ii=0
                for dataset in data:
                    plot(dataset[:,0],[np.sum(i) for i in zip(dataset[:,plotstart[iPlot]],dataset[:,plotstart[iPlot+condense]])],label=legendentries[ii])
                    ii+=1
            else:
                ii=0
                for dataset in data:
                    dataPlot0 = timederivative( TimeDerivative,TimeMode,StartTime,scale[0],dataset[:,0],dataset[:,plotstart[iPlot]])
                    dataPlot  = averagedata(Frequency,scale[0],dataPlot0[:,0],dataPlot0[:,1])
                    plot(dataPlot[:,0],dataPlot[:,1],label=legendentries[ii])
                    ii+=1
        savefig(''.join([imagename[iPlot],'.pdf'])) 
        if(dotikz==1):
            tikz_save(''.join([imagename[iPlot],'.tikz']),
                      figureheight = '\\figureheight',
                      figurewidth = '\\figurewidth'
                      )
        show()
        #close(1)
        return 
    else:
        if(plotMarker>0):
            matplotlib.rc('lines', linewidth=0.5, markersize=5)
            # plot data via loop
            if(condense>0):
                ii=0
                for dataset in data:
                    plot(dataset[:,0],[np.sum(i) for i in zip(dataset[:,plotstart[iPlot]],dataset[:,plotstart[iPlot+condense]])],label=legendentries[ii])
                    ii+=1
            else:
                ii=0
                for dataset in data:
                    dataPlot0 = timederivative( TimeDerivative,TimeMode,StartTime,scale[0],dataset[:,0],dataset[:,plotstart[iPlot]])
                    dataPlot  = averagedata(Frequency,scale[0],dataPlot0[:,0],dataPlot0[:,1])
                    plot(dataPlot[:,0],dataPlot[:,1],label=legendentries[ii],marker='o')
                    ii+=1
        else:
            fig, ax = subplots()
            axes = [ax]
            for y in range(nYAxis-1):
                # twin the x-axis twice to make independent y-axis
                axes.append(ax.twinx())
            extra_ys = len(axes[2:])
            # Make some space on the right side for the extra y-axes.
            if extra_ys>0:
                temp = 0.85
                if extra_ys<=2:
                    temp = 0.75
                elif extra_ys<=4:
                    temp = 0.6
                if extra_ys>5:
                    print('you are being ridiculous')
                fig.subplots_adjust(right=temp)
                right_additive = (0.98-temp)/float(extra_ys)
            # Move the last y-axis spine over to the right by x% of the width of the axes
            i = 1.
            for ax in axes[2:]:
                ax.spines['right'].set_position(('axes', 1.+right_additive*i))
                ax.set_frame_on(True)
                ax.patch.set_visible(False)
                ax.yaxis.set_major_formatter(matplotlib.ticker.OldScalarFormatter())
                i +=1.
            cols = []
            lines = []
            line_styles = ['-','--', '-.', ':', '.', ',', 'o', 'v', '^', '<', '>',
                       '1', '2', '3', '4', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|', '_']
            colors = ['blue','green','red','orange','yellow']
            # and finaly, plot
            #for axis y in range(nYAxis-1):
            iLine=0
            iPlot=0
            for ax,y in zip(axes,data):
                ls=line_styles[iLine]
                iPlot=iPlot+1
                print('iPlot', iPlot,iLine)
                if len(y)==1:
                    color = colors[iLine]
                    ii=0
                    for dataset in data:
                        dataPlot0 = timederivative( TimeDerivative,TimeMode,StartTime,scale[0],dataset[:,0],dataset[:,iPlot])
                        dataPlot  = averagedata(Frequency,scale[0],dataPlot0[:,0],dataPlot0[:,1])
                        lines.append(ax.plot(dataPlot[:,0],dataPlot[:,iPlot],label=legendentries[ii],linestyle =ls,color=color))
                        ax.set_ylabel('test',color=color)
                        #ax.tick_params(axis='y', colors=color)
                        ax.spines['right'].set_color(color)
                        ii +=1
                else:
                    ii=0
                    for dataset in data:
                        color = colors[ii]
                        dataPlot0 = timederivative( TimeDerivative,TimeMode,StartTime,scale[0],dataset[:,0],dataset[:,iPlot])
                        dataPlot  = averagedata(Frequency,scale[0],dataPlot0[:,0],dataPlot0[:,1])
                        if(iPlot==1):
                            lines.append(ax.plot(dataPlot[:,0],dataPlot[:,1],label=legendentries[ii],linestyle =ls,color=color))
                        else:
                            lines.append(ax.plot(dataPlot[:,0],dataPlot[:,1],linestyle =ls,color=color))
                        #ax.set_ylabel('test',color=color)
                        yName= ''.join([Variables[iLine+1], " / " ,Units[iLine+1]])
                        ax.set_ylabel(yName)
                        #ax.tick_params(axis='y', colors=color)
                        #ax.spines['right'].set_color(color)
                        ii +=1
                    #for col in y:
                    #    color = colors.next()
                    #    lines.append(ax.plot(d[col],linestyle =ls,label = col,color=color))
                    #    cols.append(col)
                    # ax.set_ylabel(', '.join(y))
                    #ax.tick_params(axis='y')
                iLine +=1
            xName= ''.join([Variables[0], " / " ,Units[0]])
            axes[0].set_xlabel(xName)
            axes[0].set_xlim([axisLimit[0],axisLimit[1]])
            iRange=2
            for ii in range(nYAxis):
                axes[ii].set_ylim([axisLimit[iRange],axisLimit[iRange+1]])
                iRange=iRange+2
            axes[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=int(nfiles/2),prop = font)

        savefig(''.join([imagename[0],'.pdf'])) 
        if(dotikz==1):
            tikz_save(''.join([imagename[iPlot],'.tikz']),
                      figureheight = '\\figureheight',
                      figurewidth = '\\figurewidth'
                      )


        show()
        #close(1)
        return 


#            # plot data via loop
#            if(condense>0):
#                ii=0
#                for dataset in data:
#                    plot(dataset[:,0],[np.sum(i) for i in zip(dataset[:,plotstart[iPlot]],dataset[:,plotstart[iPlot+condense]])],label=legendentries[ii])
#                    ii+=1
#            else:
#                ii=0
#                for dataset in data:
#                    dataPlot0 = timederivative( TimeDerivative,TimeMode,StartTime,scale[0],dataset[:,0],dataset[:,plotstart[iPlot]])
#                    dataPlot  = averagedata(Frequency,scale[0],dataPlot0[:,0],dataPlot0[:,1])
#                    plot(dataPlot[:,0],dataPlot[:,1],label=legendentries[ii])
#                    ii+=1

    #plot(data1[:,0],data1[:,1]+data1[:,2], label='$\zeta=0$')
    #plot(data1[:,0],data1[:,1],label='E_{el,1}', color='r', marker='s', markevery=360, markerfacecolor='w')
    # grid locations 
    #majoryLocator   = MultipleLocator(0.02)
    #ax.yaxis.set_major_locator(majoryLocator)
    #ax.set_xlabel("$t$ [ns]")
    #ax.xaxis.set_label_coords(0.5, -0.05)
    #ax.yaxis.set_label_coords(-0.08, 0.5)
    ##ax.set_yscale('log')
    #ax.set_xlim([0.65,0.68])
    #ax.set_ylim([0,200])
    #ax.set_xlim([0.65,0.68])
    #ax.set_ylim([13,17.5])
    # legend 
    #ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=int(nfiles/2),prop = font)
    #ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=2,prop = font)


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
