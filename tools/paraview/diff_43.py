# Programmable Filter for computing the particle distribution

Name = 'Diff'
Label = 'Difference between to DataFiles'
Help = ''

NumberOfInputs = 2
OutputDataType = 'vtkUnstructuredGrid'
ExtraXml = ''

Properties = dict()

def RequestData():
  import paraview

  for key in inputs[0].PointData.keys() :
     in1 = inputs[0].PointData[key]
     in2 = inputs[1].PointData[key]
     output.PointData.append(in1 - in2, "diff(%s)" % key)

def RequestInformation():
  from paraview import util
  pdi = self.GetInput()
