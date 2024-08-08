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
  xMax=150.,
  iVelocity=0,
  minVelo =-3E8,
  maxVelo =3E8
  )

def RequestData():
  import math
  import numpy
  import paraview
  import vtk.numpy_interface.dataset_adapter
  import vtk.numpy_interface.algorithms
  # -- this will import vtkMultiProcessController and vtkMPI4PyCommunicator

  # This script computes the particle distribution function. Missing: Selection of 
  # spacial coordinate and velocity
  deltaX=float(xMax-xMin)/float(NumberOfSpaceBins)
  deltaV=float(maxVelo-minVelo)/float(NumberOfVeloBins)
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
  nPartsMin=0
  nPartsMax=0
  nPartsIn=0
  nXmin=0
  nXmax=0
  for i in range(0, nParts):
    coord = pdi.GetPoint(i)
    pos   = coord[iDirect]
    if (iVelocity!=3):
      velo  = pdi.GetPointData().GetArray("Velocity").GetValue(3*i + iVelocity)
    else:
      vx =pdi.GetPointData().GetArray("Velocity").GetValue(3*i    )
      vy =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 1)
      vz =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 2)
      velo=math.sqrt(vx**2+vy**2+vz**2)
    # check x range
    if(xMin>pos):
      nXmin=nXmin+1
    elif(xMax<pos):
      nMax=nXmax+1
    else:
      # particle in x-range
      if(minVelo>velo):
        nPartsMin=nPartsMin+1
      elif(maxVelo<velo):
        nPartsMax=nPartsMax+1
      else:
        # compute position in 2d-space pos-velo array
        ipos = int((pos-xMin)/deltaX)#+1
        ivelo= int((velo-minVelo)/deltaV)#+1
        nPartsIn=nPartsIn+1
        PDF[ipos,ivelo]=PDF[ipos,ivelo]+1.0#/float(nParts)
  array=vtk.vtkFloatArray()
  nPartsOut=nParts-nPartsIn
  if(nXmin>0) or (nXmax>0):
    print " Particles out of coordinate range."
    print " nMinOut: ", nXmin
    print " nMaxOut: ", nXmax
    print " Percent coord out: ", float(nXmin+nXmax)/float(nParts)*100.0
  if(nPartsMin>0) or (nPartsMax>0):
    print " Particles out of velocity range. Velocity truncated!!!"
    print " nPartsMin: ", nPartsMin
    print " nPartsMax: ", nPartsMax
    print " Percent velo out of nPartsIn:     ", float(nPartsMin+nPartsMax)/float(nPartsIn)*100.0
    print " Percent velo out of nTotalParts: ", float(nPartsMin+nPartsMax)/float(nParts)*100.0
  if(nPartsOut>0):
    print " nPartsIn: ", nPartsIn
    print " total out: ", nPartsOut
    print " Percent nPartIn:  ", float(nPartsIn)/float(nParts)*100.0
    print " Percent nPartOut: ", float(nPartsOut)/float(nParts)*100.0
  array.SetName("PDF")
  array.SetNumberOfComponents(1)
  ncells  = NumberOfSpaceBins*NumberOfVeloBins
  array.SetNumberOfTuples(ncells)
  pdo.GetCellData().AddArray(array)
  ipos=0
  for j in range(0,NumberOfVeloBins):
     for i in range(0,NumberOfSpaceBins):
      # caution: transpoesed index because of storage
      array.SetValue(ipos,PDF[i,j]/float(nPartsIn))
      ipos=ipos+1

def RequestInformation():
  from paraview import util
  pdi = self.GetInput()
  util.SetOutputWholeExtent(self, [0,NumberOfSpaceBins,0,NumberOfVeloBins,0,0])
