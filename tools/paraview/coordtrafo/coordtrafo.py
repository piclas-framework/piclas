# Programmable Filter for computing the particle distribution

Name = 'Coordinate Trafo'
Label = 'Coordinate Trafo'
Help = 'Computes vectors in cylindrical or spherical coordinates'

NumberOfInputs = 1
#OutputDataType = 'vtkPointsData'
ExtraXml = ''

Properties = dict(
  WhichField=0,
  Trafo=0
  )

def RequestData():
  import math
  import numpy
  import paraview
  import vtk.numpy_interface.dataset_adapter
  import vtk.numpy_interface.algorithms
  # -- this will import vtkMultiProcessController and vtkMPI4PyCommunicator
  input=self.GetInputDataObject(0,0)
  output= self.GetOutputDataObject(0)
  if input.IsA("vtkMultiBlockDataSet"):
      # here: new format with vtk-multiblock
      print(" plugin is new-plugin with vtkMultiBlockDataSet")
      output.CopyStructure(input)
      iter = input.NewIterator()
      iter.UnRegister(None)
      iter.InitTraversal()
      pdi=iter.GetCurrentDataObject()
      pdo=pdi.NewInstance()
      pdo.UnRegister(None)
      output.SetDataSet(iter, pdo)
  else:
      # old format without multiblock
      print(" using old paraview plugin.")
      pdi=input.GetInput() 
      pdo = self.GetOutput() 
  #pdo = self.GetOutput()
  pdo.ShallowCopy(pdi)
  numPoints = pdi.GetNumberOfPoints()
  newField=vtk.vtkDoubleArray()
  if WhichField == 0:
      VarName    = "ElectricField"
  elif WhichField == 1:
      VarName    = "MagneticField"
  elif WhichField == 2:
      VarName    = "Velocity"
  newField.SetNumberOfComponents(3)
  if Trafo == 0:
      # cylindrical coord
      VarNameOut=''.join([VarName,'_CC'])
      newField.SetComponentName(0,'r')
      newField.SetComponentName(1,'theta')
      newField.SetComponentName(2,'z')
  else:
      # spherical coord
      VarNameOut=''.join([VarName,'_SC'])
      newField.SetComponentName(0,'r')
      newField.SetComponentName(1,'theta')
      newField.SetComponentName(2,'phi')
  newField.SetName(VarNameOut)
  newField.SetNumberOfTuples(numPoints)
  if Trafo == 0:
      # computes trafo in cylindrical coordinates
      for i in range(0, numPoints):
          coord = pdi.GetPoint(i)
          x, y, z  = coord[:3]
          FieldX = pdi.GetPointData().GetArray(VarName).GetValue(3*i  )
          FieldY = pdi.GetPointData().GetArray(VarName).GetValue(3*i+1)
          FieldZ = pdi.GetPointData().GetArray(VarName).GetValue(3*i+2)
          theta=math.atan2(y,x)
          FieldTheta=-FieldX*math.sin(theta)+FieldY*math.cos(theta)
          FieldR    = FieldX*math.cos(theta)+FieldY*math.sin(theta)
          #newData.InsertNextValue(theta)
          newField.SetValue(3*i  ,FieldR)
          newField.SetValue(3*i+1,FieldTheta)
          newField.SetValue(3*i+2,FieldZ)
  else:
      # spherical coord 
      for i in range(0, numPoints):
          coord = pdi.GetPoint(i)
          x, y, z  = coord[:3]
          FieldX = pdi.GetPointData().GetArray(VarName).GetValue(3*i  )
          FieldY = pdi.GetPointData().GetArray(VarName).GetValue(3*i+1)
          FieldZ = pdi.GetPointData().GetArray(VarName).GetValue(3*i+2)
          rho=math.sqrt(x*x+y*y+z*z)
          theta=math.acos(z/rho)
          phi=math.atan2(y,x)
          FieldRho=math.sin(theta)*math.cos(phi)*FieldX+math.sin(theta)*math.sin(phi)*FieldY+math.cos(theta)*FieldZ
          FieldTheta=math.cos(theta)*math.cos(phi)*FieldX+math.cos(theta)*math.sin(phi)*FieldY-math.sin(theta)*FieldZ
          FieldPhi=-math.sin(phi)*FieldX+math.cos(phi)*FieldY
          #newData.InsertNextValue(theta)
          newField.SetValue(3*i  ,FieldRho)
          newField.SetValue(3*i+1,FieldTheta)
          newField.SetValue(3*i+2,FieldPhi)
  pdo.GetPointData().AddArray(newField)

def RequestInformation():
  from paraview import util
  pdi = self.GetInput()
