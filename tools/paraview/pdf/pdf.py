# Programmable Filter for computing the particle distribution

Name = 'PDF'
Label = 'Propability Distribution Function'
Help = 'Computes the propability distribution function.'

NumberOfInputs = 1
OutputDataType = 'vtkImageData'
ExtraXml = ''

Properties = dict(
  NumberOfVeloBins = 250,
  NumberOfSpaceBins = 250,
  iDirect = 0,
  xMin=0.,
  xMax=12.57,
  iVelocity=0,
  maxVelo =3E8
  )

def RequestData():
 import math
 import numpy
 from paraview.vtk.dataset_adapter import numpyTovtkDataArray
 # This script computes the particle distribution function. Missing: Selection of 
 # spacial coordinate and velocity
 deltaX=(xMax-xMin)/float(NumberOfSpaceBins)
 deltaV=(maxVelo-minVelo)/float(NumberOfVeloBins)
 # input
 pdi = self.GetInput()
 nParts = pdi.GetNumberOfPoints()
 # output
 pdo = self.GetOutputDataObject(0)
 pdo = self.GetOutput()
 # generate 2d grid
 pdo.SetDimensions(NumberOfSpaceBins+1,NumberOfVeloBins+1,0)
 deltaXplot=1./float(NumberOfSpaceBins)
 deltaVplot=1./float(NumberOfVeloBins)
 if(maxVelo==-minVelo):
    pdo.SetOrigin(0,-0.5,0.)
 else:
    pdo.SetOrigin(0,0.0,0.)
 pdo.SetSpacing(deltaXplot,deltaVplot,1.)
 
 # On ParaView 3.98, 4.0 and 4.1
 pdo.SetExtent(0,NumberOfSpaceBins,0,NumberOfVeloBins,0,1)
 
 PDF = numpy.zeros((NumberOfSpaceBins, NumberOfVeloBins), dtype='float')
 # generate array
 # loop over all particles 
 for i in range(0, nParts):
   coord = pdi.GetPoint(i)
   pos   = coord[iDirect]
   velo  = pdi.GetPointData().GetArray("Velocity").GetValue(3*i + iVelocity)
   # check max value
   if(xMin>pos):
     print " Position < xmin. iPart ", i
   if(xMax<pos):
     print " Position > xmax. iPart  ", i
   if(minVelo>velo):
     print " Velocity < minVelo. iPart ", i
   if(maxVelo<velo):
     print " Velocity > maxVelo. iPart ", i
   # compute position in 2d-space pos-velo array
   ipos = int((pos-xMin)/deltaX)#+1
   ivelo= int((velo-minVelo)/deltaV)#+1
   PDF[ipos,ivelo]=PDF[ipos,ivelo]+1.#/float(nParts)
 array=vtk.vtkFloatArray()
 array.SetName("PDF")
 array.SetNumberOfComponents(1)
 ncells  = NumberOfSpaceBins*NumberOfVeloBins
 array.SetNumberOfTuples(ncells)
 pdo.GetCellData().AddArray(array)
 ipos=0
 for j in range(0,NumberOfVeloBins):
    for i in range(0,NumberOfSpaceBins):
     # caution: transpoesed index because of storage
     array.SetValue(ipos,PDF[i,j]/float(nParts))
     ipos=ipos+1

def RequestInformation():
 from paraview import util
 pdi = self.GetInput()
 util.SetOutputWholeExtent(self, [0,NumberOfSpaceBins,0,NumberOfVeloBins,0,0])
