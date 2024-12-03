# Programmable Filter for computing the particle distribution

Name = 'PDF'
Label = 'Probability Distribution Function'
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
  minVelo =-3E8,
  maxVelo =3E8
  )

def RequestData():
  import math
  import numpy
  import paraview
  import vtk.numpy_interface.dataset_adapter
  import vtk.numpy_interface.algorithms
  #from mpi4py import MPI
  try:
      from vtk.vtkParallelCore import vtkMultiProcessController
      from vtk.vtkParallelMPI4Py import vtkMPI4PyCommunicator
  except ImportError:
      vtkMultiProcessController = None
      vtkMPI4PyCommunicator = None
  # -- this will import vtkMultiProcessController and vtkMPI4PyCommunicator
  #if controller is None and vtkMultiProcessController is not None:
  #    controller = vtkMultiProcessController.GetGlobalController()
  controller = vtkMultiProcessController.GetGlobalController()
  nProcs =controller.GetNumberOfProcesses()
  print ' nProcs: ', nProcs
  if controller and controller.IsA("vtkMPIController") and controller.GetNumberOfProcesses() > 1:
      from mpi4py import MPI
      comm = vtkMPI4PyCommunicator.ConvertToPython(controller.GetCommunicator())
      rank = comm.Get_rank()
  else:
      rank = 0

  # This script computes the particle distribution function. Missing: Selection of 
  # spacial coordinate and velocity
  deltaX=float(xMax-xMin)/float(NumberOfSpaceBins)
  deltaV=float(maxVelo-minVelo)/float(NumberOfVeloBins)
  # input
  input=self.GetInputDataObject(0,0)
  if input.IsA("vtkMultiBlockDataSet"):
      # here: new format with vtk-multiblock
      print(" vtkMultiBlockDataSet")
      iter = input.NewIterator()
      iter.UnRegister(None)
      iter.InitTraversal()
      pdi=iter.GetCurrentDataObject()
  else:
      # old format without multiblock
      pdi=input.GetInput()  
  nParts= pdi.GetNumberOfPoints()
  if nProcs>1:
      #nTotalParts = numpy.array(0, 'i')
      #nTotalParts=[]
      #comm.Allreduce([nParts, MPI.INT], [nTotalParts, MPI.INT], op=MPI.SUM)
      nTotalParts=comm.allreduce(nParts  , op=MPI.SUM)
  else:
      nTotalParts=nParts
  if rank ==0:
      print ' nTotalParts:', nTotalParts
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
    if (iVelocity==3):
        # vabs
        vx =pdi.GetPointData().GetArray("Velocity").GetValue(3*i    )
        vy =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 1)
        vz =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 2)
        velo=math.sqrt(vx**2+vy**2+vz**2)
    elif(iVelocity==4):
       # v-tang to x
       vy =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 1)
       vz =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 2)
       velo=math.sqrt(vy**2+vz**2)
    elif(iVelocity==5):
       # vabv-tang to y
       vx =pdi.GetPointData().GetArray("Velocity").GetValue(3*i    )
       vz =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 2)
       velo=math.sqrt(vx**2+vz**2)
    elif(iVelocity==6):
       # vabv-tang to z
       vx =pdi.GetPointData().GetArray("Velocity").GetValue(3*i    )
       vy =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 1)
       velo=math.sqrt(vx**2+vy**2)
    else:
       # velocity in x,y or z
       velo  = pdi.GetPointData().GetArray("Velocity").GetValue(3*i + iVelocity)

    #if (iVelocity!=3):
    #  velo  = pdi.GetPointData().GetArray("Velocity").GetValue(3*i + iVelocity)
    #else:
    #  vx =pdi.GetPointData().GetArray("Velocity").GetValue(3*i    )
    #  vy =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 1)
    #  vz =pdi.GetPointData().GetArray("Velocity").GetValue(3*i + 2)
    #  velo=math.sqrt(vx**2+vy**2+vz**2)
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
  totalPDF = numpy.zeros((NumberOfSpaceBins, NumberOfVeloBins), dtype='float')
  if nProcs>1:
      # mpi stuff
      nTotalPartsIn=comm.allreduce(nPartsIn  , op=MPI.SUM)
      totalPDF=comm.allreduce(PDF  , op=MPI.SUM)
      for j in range(0,NumberOfVeloBins):
           for i in range(0,NumberOfSpaceBins):
               totalPDF[i,j]=totalPDF[i,j]/float(nTotalPartsIn)
  else:
      nTotalPartsIn=nPartsIn
      for j in range(0,NumberOfVeloBins):
         for i in range(0,NumberOfSpaceBins):
             totalPDF[i,j]=PDF[i,j]/float(nTotalPartsIn)
  array=vtk.vtkFloatArray()
  if nProcs>1:
      # sum up to total
      nTotalXmin=comm.reduce(nXmin , op=MPI.SUM)
      nTotalXmax=comm.reduce(nXmax , op=MPI.SUM)
      nTotalPartsMin=comm.reduce(nPartsMin , op=MPI.SUM)
      nTotalPartsMax=comm.reduce(nPartsMax , op=MPI.SUM)
  else:
      nTotalXmin=nXmin
      nTotalXmax=nXmax
      nTotalPartsMin=nPartsMin
      nTotalPartsMax=nPartsMax
  nTotalPartsOut=nTotalParts-nTotalPartsIn
  if rank==0:
      # output
      if(nTotalXmin>0) or (nTotalXmax>0):
        print " Particles out of coordinate range."
        print " nMinOut: ", nTotalXmin
        print " nMaxOut: ", nTotalXmax
        print " Percent coord out: ", float(nTotalXmin+nTotalXmax)/float(nTotalParts)*100.0
      if(nTotalPartsMin>0) or (nTotalPartsMax>0):
        print " Particles out of velocity range. Velocity truncated!!!"
        print " nPartsMin: ", nTotalPartsMin
        print " nPartsMax: ", nTotalPartsMax
        print " Percent velo out of nPartsIn:     ", float(nTotalPartsMin+nTotalPartsMax)/float(nTotalPartsIn)*100.0
        print " Percent velo out of nParts:       ", float(nTotalPartsMin+nTotalPartsMax)/float(nTotalParts)*100.0
      if(nTotalPartsOut>0):
        print " nPartsIn:  ", nTotalPartsIn
        print " total out: ", nTotalPartsOut
        print " Percent nPartIn:  ", float(nTotalPartsIn) /float(nTotalParts)*100.0
        print " Percent nPartOut: ", float(nTotalPartsOut)/float(nTotalParts)*100.0
  array.SetName("PDF")
  array.SetNumberOfComponents(1)
  ncells  = NumberOfSpaceBins*NumberOfVeloBins
  array.SetNumberOfTuples(ncells)
  pdo.GetCellData().AddArray(array)
  ipos=0
  for j in range(0,NumberOfVeloBins):
     for i in range(0,NumberOfSpaceBins):
      # caution: transpoesed index because of storage
      #array.SetValue(ipos,totalPDF[i,j]/float(nPartsIn))
      array.SetValue(ipos,totalPDF[i,j])
      ipos=ipos+1

def RequestInformation():
  from paraview import util
  pdi = self.GetInput()
  util.SetOutputWholeExtent(self, [0,NumberOfSpaceBins,0,NumberOfVeloBins,0,0])
