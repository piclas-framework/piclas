RenderView1 = GetRenderView()
print RenderView1.CameraParallelScale
#RenderView1.CameraParallelScale = 0.5
#RenderView1.CameraParallelScale = 0.0001 # small means zoom in !
#print RenderView1.CameraParallelScale
Render()

from paraview.simple import *
import os
import glob
import sys
from paraview import simple
view = GetActiveView()
[width,height]= view.ViewSize

print view.ViewSize
if height%2==1:
    height=height+1
if width%2==1:
    width=width+1
view.ViewSize = [width, height]
print view.ViewSize
# [1031, 687] Standard Settings for Samsung E172 Laptop
view.CameraParallelScale=RenderView1.CameraParallelScale
print view.CameraParallelScale

# get active pipeline object (i.e. the statefile which has to be selected by hand)
#S = GetActiveSource()
#print S
#P = S.FileName
#Spath,Sfile = os.path.split(P)
#os.chdir(Spath)

pwd=os.getcwd() # get current path
print pwd

PSources=GetSources() # search through existing pipeline objects
#sys.exit() # stop and close paraview
#print PSources

# get all state files in current folder
List = glob.glob(pwd+"/*State*.h5")
List = sorted(List)
#print List
K=0
for file in List:
    K=K+1 # count all state files
N=0
for file in List:
    Fpath,Ffile = os.path.split(file)    
    N=N+1
    print N," of ",K,"  (",Ffile,")"

    # old
    #TE23 = GetActiveSource()
    #TE23.FileName = file # das State file muss ausgewählt sein
    # new
    for object in PSources:
        # check all objects in pipiline and find out which are State files
        if hasattr(PSources[object], 'FileName'):
          if (PSources[object].FileName.find("State")>0.0 and 
              PSources[object].FileName.find(".h5") > 0.0):
              PSources[object].FileName=file
    
    SetActiveView(view)
    Render()
    simple.WriteImage("new"+"%0.4d" % N+".bmp")
    #N=N-1

print "Creating Video ..."
os.system('avconv -r 10 -y -i new%4d.bmp -b:v 4000k output.mp4')
print "Done."
