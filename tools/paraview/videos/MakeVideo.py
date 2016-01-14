RenderView1 = GetRenderView()
RenderView1.CameraParallelScale = 0.5
Render()

from paraview.simple import *
import os
import glob
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


S = GetActiveSource()
print S
P = S.FileName
Spath,Sfile = os.path.split(P)
os.chdir(Spath)

pwd=os.getcwd() 
print pwd

List = glob.glob(pwd+"/*State*.h5")
List = sorted(List)
#print List
K=0
for file in List:
    K=K+1
N=0
for file in List:
    Fpath,Ffile = os.path.split(file)    
    N=N+1
    print N," of ",K,"  (",Ffile,")"
    TE23 = GetActiveSource()
    TE23.FileName = file # das State file muss ausgewählt sein
    SetActiveView(view)
    Render()
    simple.WriteImage("new"+"%0.3d" % N+".bmp")
    #N=N-1

print "Creating Video ..."
os.system('avconv -r 10 -y -i new%3d.bmp -b:v 4000k output.mp4')
print "Done."
