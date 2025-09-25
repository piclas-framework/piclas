#==================================================================================================================================
# Copyright (c) 2010 - 2025 Raphael Tietz
#
# This file is part of PICLas (piclas.boltzplatz.eu/piclas/piclas). PICLas is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.
#
# PICLas is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
#
# You should have received a copy of the GNU General Public License along with PICLas. If not, see <http://www.gnu.org/licenses/>.
#==================================================================================================================================

import numpy as np
import h5py
import glob
import argparse
import random
import os
import re
import subprocess as sp
import shutil
import sys
import textwrap
import scipy.optimize

# Constants
pi=3.14159
kb=1.38E-23
Na=6.022E23

# global variables
ElemIDsorted = None
Barycenters = None

# Class Definition

class CAutoadjust:
  tend = None
  dt = None
  Mesh = None
  MPF = None
  maxPartNum = None
Autoadjust=CAutoadjust()
class CMesh:
  dxElem = None
  Length = None
  HalfWidth = None
  def nElems(self):
    return int(self.Length/self.dxElem)
Mesh = CMesh()
class CSimulation:
  tend = None
  dt = None
  MPF = None
  Median = None
  nSpecies = None
  SpeciesMass = None
  ProjectName = None
  Redosimulation = False
  PartperElemmin = None
  MCSoverMFPmax = None
  VeloTempEps = None
  MaxCollProbMax = None
  UsedValues = None
  UseMPI = None
  nMPICores = 1
  InflowSpecies = None
  InflowPressure = None
  ElemsPerWrite = None
  MaxPartNum = None
  ElectronSpecies = None
Simulation = CSimulation()
class CIteration:
  TargetVelo = None
  TargetDev = None
  CurrentVelo = None
  ResponseRatio = None
  iter = -1
  LastVelo = None
  LastLagShockSpeed = None
  CurrentLagShockSpeed = None
  StandingShock = False
  VIn=False
  VOut=False
Iteration = CIteration()
class CPath:
  piclas = None
  hopr = None
  DSMCSpecies = None
  KeepSimulationData = None
  MasterDir = None
  CurrentSimDir = None
  LastSimDir = None
  AdditionalCopyFiles = None
Path = CPath()
class CSolution:
  Shockspeed = None
  times = None
  xShock = None
  xShockValid = None
  MCSoverMFPmax = None
  TotalSimPartNummin = None
  MaxCollProb = None
  MaxVelo = None
  TooFewValues = False
Solution = CSolution
RestartStanding = False
DoRestart = False

class ReadinValue:
  """
  Readin Values
  """
  name = None
  datatype = None
  default = None
  value = None
  help = None
  def __init__(self, name, datatype, default, help):
    self.name = name.lower()
    self.datatype = datatype
    self.default = default
    self.isdefault = None
    self.help = help

class CReadinValues:
  l = []
  def SetValues(self):
    global ShockSimulationParameter
    for Entry in self.l:
      if Entry.name.lower() in ShockSimulationParameter.keys():
        Entry.value = convert(ShockSimulationParameter[Entry.name],Entry.datatype)
        Entry.isdefault = False
      else:
        if Entry.default == None:
          raise Exception('%s has to be set' % Entry.name)
        else:
          Entry.value = convert(Entry.default,Entry.datatype)
          Entry.isdefault = True
    ShockSimulationParameter = None
  def GetValue(self,name):
    for Entry in self.l:
      if Entry.name.lower() == name.lower():
        if Entry.isdefault:
          print('  '+Entry.name + ' '*(100-len(Entry.name)-len(str(Entry.value)))+ str(Entry.value) + ' DEFAULT')
          # print(str(type(Entry.value)) + '  ' + str(Entry.value) + '  ' + Entry.name + '  Default')
        else:
          print('  '+Entry.name + ' '*(100-len(Entry.name)-len(str(Entry.value)))+ str(Entry.value) + '  CUSTOM')
          # print(str(type(Entry.value)) + '  ' + str(Entry.value) + '  ' + Entry.name + '  Custrom')
        output = Entry.value
        self.l.remove(Entry)
        return output
    raise Exception('parameter %s is not in ReadinValues list' % name)
  def ListUnusedValues(self):
    if len(self.l) == 0:
      return
    print('Variables set but unused:')
    for Entry in self.l:
      print(' '*5+Entry.name)
    self.l=[]
  def ShocktubeiniHelp(self):
    helpstr = textwrap.dedent('''\
      Help for Shocktube.ini:

      Set all neccesary piclas simulation parameter in Shocktube.ini like nSpecies or useOctree.
      Set Particles-Species[$]-nInit=1 for all incident Species and set temperature(s), numberdensity and velovec for those inits.
      All Shocktube related parameter are listed below with their default value. If default is None there is no default value and has to be set.

      ''')
    for Entry in self.l:
      helpentry = Entry.help.split('\n')
      for i in range(len(helpentry)):
        if i==0:
          helpstr += '%s = %s '%(Entry.name,str(Entry.default))+' '*(50-len(Entry.name) - 4 - len(str(Entry.default))) + helpentry[0] + '\n'
        else:
          helpstr += ' '*50 + helpentry[i] + '\n'
    return helpstr
ReadinValues = CReadinValues()

# Function Definition

def CreateTitle():
  """
  Print Title
  """
  CharMaxFreestream = 15
  CharMaxShocked = 5
  TitleLength = 128
  TitleHeight = 6

  Title=[]

  for i in range(TitleHeight):
    Title.append('')
    for ii in range(int(TitleLength/2)):
      Title[i]+=TitleCharsFreestream(random.randint(0, CharMaxFreestream-1))
    for ii in range(int(TitleLength/2)+1,TitleLength+1):
      Title[i]+=TitleCharsShocked(random.randint(0, CharMaxShocked-1))


  print('=' * TitleLength)
  print('')
  print(' '*int((TitleLength-117)/2) + ' _______           _______  _______  _        _______  ______   _______ _________ _______  _        _______  _______ ')
  print(' '*int((TitleLength-117)/2) + '(  ____ \|\     /|(  ___  )(  ____ \| \    /\(  ____ \(  __  \ (  ____ )\__   __/(  ____ \( \      (  ___  )(  ____ \\')
  print(' '*int((TitleLength-117)/2) + '| (    \/| )   ( || (   ) || (    \/|  \  / /| (    \/| (  \  )| (    )|   ) (   | (    \/| (      | (   ) || (    \/')
  print(' '*int((TitleLength-117)/2) + '| (_____ | (___) || |   | || |      |  (_/ / | (__    | |   ) || (____)|   | |   | |      | |      | (___) || (_____ ')
  print(' '*int((TitleLength-117)/2) + '(_____  )|  ___  || |   | || |      |   _ (  |  __)   | |   | ||  _____)   | |   | |      | |      |  ___  |(_____  )')
  print(' '*int((TitleLength-117)/2) + '      ) || (   ) || |   | || |      |  ( \ \ | (      | |   ) || (         | |   | |      | |      | (   ) |      ) |')
  print(' '*int((TitleLength-117)/2) + '/\____) || )   ( || (___) || (____/\|  /  \ \| (____/\| (__/  )| )      ___) (___| (____/\| (____/\| )   ( |/\____) |')
  print(' '*int((TitleLength-117)/2) + '\_______)|/     \|(_______)(_______/|_/    \/(_______/(______/ |/       \_______/(_______/(_______/|/     \|\_______)')


  print('')
  print('')
  print('=' * TitleLength)
  for i in range(TitleHeight):
    print(Title[i])
  print('=' * TitleLength)
  print('')

def TitleCharsShocked(i):
  switcher={
          0:'%',
          1:':',
          2:'"',
          3:';',
          4:'*',
        }
  return switcher.get(i)

def TitleCharsFreestream(i):
  switcher={
          0:"'",
          1:'-',
          2:'`',
          3:'.',
          4:',',
          5:' ',
          6:' ',
          7:' ',
          8:' ',
          9:' ',
          10:' ',
          11:' ',
          12:' ',
          13:' ',
          14:' ',
        }
  return switcher.get(i)

def ArgumentsParser():
  global ReadinValues
  parser = argparse.ArgumentParser(description='Python programm to run shock tube simulations automatically.\n',epilog=ReadinValues.ShocktubeiniHelp(),formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('InputFile', help='File containing input parameters')
  parser.add_argument('--Equicon',action='store_true', help='Calculate only equilibrium condition')
  parser.add_argument('--restart','-r',action='store', help='Restart from the given simulation (directory)')
  args = parser.parse_args()
  return args

def ReadShocktubeIni(InputFile):

  global SimulationParameter, ShockSimulationParameter, ReadinValues
  if not os.path.exists(InputFile) :
    raise Exception("input file '%s' not found." % InputFile)
  SimulationParameter={}
  ShockSimulationParameter={}
  with open(InputFile) as f :
    for line in f.readlines() :   # iterate over all lines of the file
        line = re.sub(r"\s+", "", line)        # remove all whitespaces ("\s" is the whitespac symbol)
        line = re.sub(r"\\s", " ", line)       # add new whitespaces for all occurrances of "\s" in the string ("\s" is NOT the whitespace symbol here)
        if line.startswith('!') : continue     # skip lines starting with a comment
        line = line.split('!')[0]              # remove comments
        if '=' in line :                                # reading of option finished -> go on with next line

          if line.lower().startswith('shock-'):
            (key,value) = line.split('=',1)          # split line at '='
            ShockSimulationParameter[key.lower()]=value
            continue

          # if line.lower().startswith('Particles') :
          if line.lower().startswith('particles-symmetry-order') :
            print('WARNING: Particles-Symmetry-Order will be set as 1 and is looped in the input file %s' % InputFile)
            continue
          (key,value) = line.split('=',1)          # split line at '='
          SimulationParameter[key.lower()]=value
          # option = Option(key,value)               # generate new Option with a value
          # SimulationParameters.append(option)      # append option to options list, where
          continue
  ReadinValues.SetValues()

def convert(Value,datatype):
  if datatype == str:
    return str(Value)
  elif datatype == int:
    return int(Value)
  elif datatype == float:
    return float(Value)
  elif datatype == bool:
    if Value.lower() in ['true','t','1']:
      return True
    else:
      return False

def GenerateReadinValues():

  global ReadinValue
  global ReadinValues
  ReadinValues = CReadinValues()
  ReadinValues.l.append(ReadinValue("Shock-Autoadjust-Mesh".lower(),bool,'False','If true: mesh length will be adjusted by last simulation results, might be unstable'))
  ReadinValues.l.append(ReadinValue("Shock-Autoadjust-dt".lower(),bool,'False','If true: dt will be adjusted by last simulation results.'))
  ReadinValues.l.append(ReadinValue("Shock-Autoadjust-tend".lower(),bool,'False','If true: tend will be adjusted by last simulation results.'))
  ReadinValues.l.append(ReadinValue("Shock-Autoadjust-maxPartNum".lower(),bool,'False','If true: maxpartnum will be adjusted by last simulation results.'))
  ReadinValues.l.append(ReadinValue("Shock-Autoadjust-MPF".lower(),bool,'False','If true: MPF will be adjusted by last simulation results.'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-dt".lower(),float,None,'dt of the simulation, if autoadjust then dt of the first simulation'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-tend".lower(),float,None,'tend of the simulation, if autoadjust then tend of the first simulation'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-MPF".lower(),float,None,'MPF of the simulation, if autoadjust then MPF of the first simulation'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-MaxPartNum".lower(),int,-1,'MaxPartNum of the simulation, if autoadjust then MaxPartNum of the first simulation'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-Median".lower(),float,0.01,'allowed relative median in Integrated File (smoothness condition of overall solution)'))
  ReadinValues.l.append(ReadinValue("Shock-Mesh-dxElem".lower(),float,None,'x extend of a Element/Cell'))
  ReadinValues.l.append(ReadinValue("Shock-Mesh-Length".lower(),float,None,'length of the simulation, if autoadjust then length of the first simulation'))
  ReadinValues.l.append(ReadinValue("Shock-Mesh-HalfWidth".lower(),float,'0.05','half y and z extend of the mesh'))
  ReadinValues.l.append(ReadinValue("Shock-PiclasPath".lower(),str,None,'path to the piclas executeable or cl command'))
  ReadinValues.l.append(ReadinValue("Shock-HoprPath".lower(),str,None,'path to the hopr executeable or cl command'))
  ReadinValues.l.append(ReadinValue("Shock-DSMCSpecies".lower(),str,'None','Name of the DSMCSpecies file. has not to be set if  DSMCSpecies file is not needed'))
  ReadinValues.l.append(ReadinValue("Shock-Target-Velo".lower(),float,None,'Intended Velocity of the Shock'))
  ReadinValues.l.append(ReadinValue("Shock-Target-dev".lower(),float,None,'Allowed deviation of the target velocity'))
  ReadinValues.l.append(ReadinValue("Shock-Start-Velo".lower(),float,'0.0','Initial velocity of the flow field. If 0 then "Shock-Target-Velo" is used'))
  ReadinValues.l.append(ReadinValue("Shock-ResponseRatio".lower(),float,'1.0','= delta Inflowvelocity/delta Shockspeed'))
  ReadinValues.l.append(ReadinValue("Shock-KeepAllSimulationData".lower(),bool,'False','If true: for every simulation a new folder is created.\nIf false: last simulation data will be deleted'))
  ReadinValues.l.append(ReadinValue("Shock-PartperElemmin".lower(),float,10,'Only if "Shock-Autoadjust-MPF": lowest Simpartnum allowed in Elems'))
  ReadinValues.l.append(ReadinValue("Shock-MCSoverMFPmax".lower(),float,0.2,'Only if "Shock-Autoadjust-MPF": highest MCSoverMFP allowed in Elems'))
  ReadinValues.l.append(ReadinValue("Shock-VeloTempEps".lower(),float,5,'Only if "Shock-Autoadjust-dt": maximum considered multiple thermal velocity of the mean thermal velocity'))
  ReadinValues.l.append(ReadinValue("Shock-MaxCollProbMax".lower(),float,0.9,'Only if "Shock-Autoadjust-dt": highest MaxCollProbMax allowed in Elems'))
  ReadinValues.l.append(ReadinValue("Shock-UseMPI".lower(),bool,'False','If piclas should run with mpi'))
  ReadinValues.l.append(ReadinValue("Shock-nMPICores".lower(),int,1,'Only if "Shock-UseMPI": number of used cores'))
  ReadinValues.l.append(ReadinValue("Shock-AdditionalCopyFiles".lower(),str,'','Only if "KeepAllSimulationData": additional neccessary files like electronic data base'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-ElemsPerWrite".lower(),float,1,'Average travel distance of the shock in Elems between piclas macrowrites'))
  ReadinValues.l.append(ReadinValue("Shock-Simulation-ElectronSpecies".lower(),int,-1,'Species number of electrons =-1 if there is no electron species'))

def Init(KeepOldData,RestartFolder):
  """
  First Initailization of Parameter
  """
  global nSpecies, SimulationParameter, ReadinValues
  global Autoadjust, Mesh, Simulation, Iteration, Path, Solution


  # Autoadjust
  Autoadjust.tend = ReadinValues.GetValue("Shock-Autoadjust-tend")
  Autoadjust.dt = ReadinValues.GetValue("Shock-Autoadjust-dt")
  Autoadjust.Mesh = ReadinValues.GetValue("Shock-Autoadjust-Mesh")
  Autoadjust.MPF = ReadinValues.GetValue("Shock-Autoadjust-MPF")
  Autoadjust.maxPartNum = ReadinValues.GetValue("Shock-Autoadjust-maxPartNum")

  # Mesh
  Mesh.dxElem = ReadinValues.GetValue("Shock-Mesh-dxElem")
  Mesh.Length = ReadinValues.GetValue("Shock-Mesh-Length")
  Mesh.HalfWidth = ReadinValues.GetValue("Shock-Mesh-HalfWidth")

  # Simulation
  Simulation.tend = ReadinValues.GetValue("Shock-Simulation-tend")
  Simulation.dt = ReadinValues.GetValue("Shock-Simulation-dt")
  Simulation.MPF = ReadinValues.GetValue("Shock-Simulation-MPF")
  Simulation.Median = ReadinValues.GetValue("Shock-Simulation-Median")
  Simulation.nSpecies = int(SimulationParameter['Part-nSpecies'.lower()])
  Simulation.ElemsPerWrite = ReadinValues.GetValue("Shock-Simulation-ElemsPerWrite")
  Simulation.MaxPartNum = ReadinValues.GetValue("Shock-Simulation-MaxPartNum")
  Simulation.ElectronSpecies = ReadinValues.GetValue("Shock-Simulation-ElectronSpecies")
  if ((not Simulation.ElectronSpecies==-1) and (Simulation.ElectronSpecies<1 or Simulation.ElectronSpecies>Simulation.nSpecies)):
    raise Exception("'Shock-Simulation-ElectronSpecies' has to be -1 for no electron species or species ID of the electrons")
  if not "ProjectName".lower() in SimulationParameter.keys():
    SimulationParameter["ProjectName".lower()] = "Shocktube"
  Simulation.ProjectName = SimulationParameter["ProjectName".lower()]
  Simulation.SpeciesMass = np.zeros(Simulation.nSpecies+1)
  for iSpec in range(1,Simulation.nSpecies+1):
    Simulation.SpeciesMass[iSpec] = SimulationParameter[("Part-Species"+str(iSpec)+"-MassIC").lower()]
  if Autoadjust.MPF:
    Simulation.PartperElemmin = ReadinValues.GetValue("Shock-PartperElemmin")
    Simulation.MCSoverMFPmax = ReadinValues.GetValue("Shock-MCSoverMFPmax")
  if Autoadjust.dt:
    Simulation.VeloTempEps = ReadinValues.GetValue("Shock-VeloTempEps")
    Simulation.MaxCollProbMax = ReadinValues.GetValue("Shock-MaxCollProbMax")
  Simulation.UseMPI = ReadinValues.GetValue("Shock-UseMPI")
  if Simulation.UseMPI:
    Simulation.nMPICores = ReadinValues.GetValue("Shock-nMPICores")
    if Simulation.nMPICores<2:
      raise Exception("'Shock-nMPICores' has to be at least 2 or higher")
  else:
    Simulation.nMPICores = 1

  # Iteration
  Iteration.TargetVelo = ReadinValues.GetValue("Shock-Target-Velo")
  Iteration.TargetDev = ReadinValues.GetValue("Shock-Target-dev")
  Iteration.CurrentVelo = ReadinValues.GetValue("Shock-Start-Velo")
  if Iteration.CurrentVelo == 0.0:
    print("use 'Shock-Target-Velo' as 'Shock-Start-Velo'")
    Iteration.CurrentVelo = Iteration.TargetVelo
  Iteration.ResponseRatio = ReadinValues.GetValue("Shock-ResponseRatio")

  # Paths
  Path.piclas = os.path.expanduser(ReadinValues.GetValue("Shock-PiclasPath"))
  Path.hopr = os.path.expanduser(ReadinValues.GetValue("Shock-HoprPath"))
  Path.DSMCSpecies = os.path.expanduser(ReadinValues.GetValue("Shock-DSMCSpecies"))
  Path.KeepSimulationData = ReadinValues.GetValue("Shock-KeepAllSimulationData")
  if Path.KeepSimulationData:
    Path.AdditionalCopyFiles = ReadinValues.GetValue("Shock-AdditionalCopyFiles").split(',')
    while '' in Path.AdditionalCopyFiles:
      Path.AdditionalCopyFiles.remove('')

  Path.MasterDir = os.getcwd()

  ReadinValues.ListUnusedValues()

  # Remove old data
  if(not KeepOldData):
    if(RestartFolder==None):
      Files = glob.glob("*State*")
      for File in Files:
        os.remove(File)
      Directories = glob.glob("iter_[0-9]*")
      for Directory in Directories:
        shutil.rmtree('%s/%s'%(Path.MasterDir,Directory))
      Directories = glob.glob("standing_iter_[0-9]*")
      for Directory in Directories:
        shutil.rmtree('%s/%s'%(Path.MasterDir,Directory))
    else:
      tmp = RestartFolder.split("_")[-1]
      if(tmp[-1]=='/'):
        tmp = tmp[:-1]
      lastIteration = int(tmp)
      if("standing" in RestartFolder):
        Directories = glob.glob("standing_iter_[0-9]*")
        for Directory in Directories:
          tmp = Directory.split("_")[-1]
          if(int(tmp)>lastIteration):
            shutil.rmtree('%s/%s'%(Path.MasterDir,Directory))
      else:
        Directories = glob.glob("standing_iter_[0-9]*")
        for Directory in Directories:
          shutil.rmtree('%s/%s'%(Path.MasterDir,Directory))
        Directories = glob.glob("iter_[0-9]*")
        for Directory in Directories:
          tmp = Directory.split("_")[-1]
          if(int(tmp)>lastIteration):
            shutil.rmtree('%s/%s'%(Path.MasterDir,Directory))


  # Set boundarys & mesh
  SimulationParameter['Part-nBounds'.lower()] = '3'

  SimulationParameter['Part-Boundary1-SourceName'.lower()] = 'BC_Sym'
  SimulationParameter['Part-Boundary1-Condition'.lower()] = 'symmetric_dim'
  SimulationParameter['Part-Boundary2-SourceName'.lower()] = 'BC_Inflow'
  SimulationParameter['Part-Boundary2-Condition'.lower()] = 'open'
  SimulationParameter['Part-Boundary3-SourceName'.lower()] = 'BC_Outflow'
  SimulationParameter['Part-Boundary3-Condition'.lower()] = 'reflective'


  SimulationParameter['MeshFile'.lower()] = SimulationParameter["ProjectName".lower()] + '_mesh.h5'
  SimulationParameter['useCurveds'.lower()] = 'F'


  # Set constants

  #SimulationParameter['c0'.lower()] = '299792458.'
  #SimulationParameter['eps'.lower()] = '8.8541878176E-12'
  #SimulationParameter['mu'.lower()] = '12.566370614e-7'

  SimulationParameter['N'.lower()] = '1'
  SimulationParameter['IniExactFunc'.lower()] = '0'

  SimulationParameter['Particles-Symmetry-Order'.lower()] = '1'
  SimulationParameter['Part-WriteMacroVolumeValues'.lower()] = 'true'

  if Autoadjust.MPF:
    SimulationParameter['Particles-DSMC-CalcQualityFactors'.lower()] = 'true'

  # Set constant particle properties
  Simulation.InflowSpecies = []
  for iSpec in range(1,Simulation.nSpecies+1):
    if ('Part-Species'+str(iSpec)+'-nInits').lower() in SimulationParameter.keys():
      if int(SimulationParameter[('Part-Species'+str(iSpec)+'-nInits').lower()])>=2:
        raise Exception("Part-Species%s-nInits has to be 1 or 0, if its present in freestream or not" % iSpec)
      if int(SimulationParameter[('Part-Species'+str(iSpec)+'-nInits').lower()])==1:
        Simulation.InflowSpecies.append(iSpec)
        print("InflowSpecies %i added" % iSpec)
        # Set Init
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-SpaceIC').lower()] = 'cuboid'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BaseVector1IC').lower()] = '(/0.,1.,0./)'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BaseVector2IC').lower()] = '(/0.,0.,1./)'
        # Set Surfaceflux
        SimulationParameter[('Part-Species'+str(iSpec)+'-nSurfaceFluxBCs').lower()] = '1'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-BC').lower()] = '2'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-velocityDistribution').lower()] = 'maxwell_lpn'
        if ('Part-Species'+str(iSpec)+'-Init1-VeloVecIC').lower() in SimulationParameter.keys():
          SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloVecIC').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloVecIC').lower()]
        else:
          raise Exception('Part-Species%s-Init1-VeloVecIC has to be set' % iSpec)
        if ('Part-Species'+str(iSpec)+'-Init1-MWTemperatureIC').lower() in SimulationParameter.keys():
          SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-MWTemperatureIC').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-MWTemperatureIC').lower()]
        else:
          raise Exception('Part-Species%s-Init1-MWTemperatureIC has to be set' % iSpec)
        if ('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower() in SimulationParameter.keys():
          SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-PartDensity').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()]
        else:
          raise Exception('Part-Species%s-Init1-PartDensity has to be set' % iSpec)
        if ('Part-Species'+str(iSpec)+'-Init1-TempRot').lower() in SimulationParameter.keys():
          SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempRot').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempRot').lower()]
        if ('Part-Species'+str(iSpec)+'-Init1-TempVib').lower() in SimulationParameter.keys():
          SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempVib').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempVib').lower()]
        if ('Part-Species'+str(iSpec)+'-Init1-TempElec').lower() in SimulationParameter.keys():
          SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempElec').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempElec').lower()]

  Simulation.InflowPressure = 0.0
  for iSpec in Simulation.InflowSpecies:
    Simulation.InflowPressure += float(SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()]) \
      * kb * float(SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-MWTemperatureIC').lower()])
  print('Simulation.InflowPressure:  ' + str(Simulation.InflowPressure))
  print(Simulation.InflowSpecies)


  Solution.MaxCollProb=0
  Solution.MaxVelo=0
  Solution.MCSoverMFPmax=0
  # Solution.TotalSimPartNummin=sys.maxint-1
  Solution.TotalSimPartNummin=sys.maxsize
  Solution.TooFewValues=False


def InitRestart(RestartFolder):
  """
  Initialize Restart Properties
  """
  global RestartStanding, DoRestart, Iteration, Mesh

  DoRestart = True
  print("")
  print("Read specified iteration %s for restart" % RestartFolder)

  if not os.path.exists(RestartFolder) :
    raise Exception("Restart directory '%s' not found." % RestartFolder)

  if("standing" in RestartFolder):
    RestartStanding=True
    Iteration.StandingShock=True
    print("Perform standing shock restart")
  else:
    RestartStanding=False

  Path.CurrentSimDir = '%s/%s' % (Path.MasterDir,RestartFolder)
  # os.chdir(Path.CurrentSimDir)

  tmp = RestartFolder.split("_")[-1]
  if(tmp[-1]=='/'):
    tmp = tmp[:-1]
  Iteration.iter = int(tmp)


  RestartParameter={}
  with open("%s/parameter.ini" % Path.CurrentSimDir) as f :
    for line in f.readlines() :   # iterate over all lines of the file
        line = re.sub(r"\s+", "", line)        # remove all whitespaces ("\s" is the whitespac symbol)
        line = re.sub(r"\\s", " ", line)       # add new whitespaces for all occurrances of "\s" in the string ("\s" is NOT the whitespace symbol here)
        if line.startswith('!') : continue     # skip lines starting with a comment
        line = line.split('!')[0]              # remove comments
        if '=' in line :                                # reading of option finished -> go on with next line
          (key,value) = line.split('=',1)          # split line at '='
          RestartParameter[key.lower()]=value
          continue


  if(Autoadjust.tend):
    Simulation.tend = float(RestartParameter['TEnd'.lower()])
    print("updated tend to %e" % Simulation.tend)
  else:
    print("No autoadjust for tend defined, use it from Shocktube.ini: %e" % Simulation.tend)
  if(Autoadjust.dt):
    Simulation.dt = float(RestartParameter['Particles-ManualTimeStep'.lower()])
    print("updated dt to %e" % Simulation.dt)
  else:
    print("No autoadjust for dt defined, use it from Shocktube.ini: %e" % Simulation.dt)
  if(Autoadjust.Mesh):
    if(Iteration.StandingShock):
      # SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-CuboidHeightIC').lower()] = '%f' % (Mesh.Length*1./2.)
      Mesh.Length = 2*float(RestartParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init1-CuboidHeightIC').lower()])
    else:
      # SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-CuboidHeightIC').lower()] = '%f' % (Mesh.Length)
      Mesh.Length = float(RestartParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init1-CuboidHeightIC').lower()])
    print("updated Mesh-Length to %e" % Mesh.Length)
  else:
    print("No autoadjust for Mesh-Length defined, use it from Shocktube.ini: %e" % Mesh.Length)
  if(Autoadjust.MPF):
    Simulation.MPF = float(RestartParameter['Part-Species1-MacroParticleFactor'.lower()])
    print("updated MPF to %e" % Simulation.MPF)
  else:
    print("No autoadjust for MPF defined, use it from Shocktube.ini: %e" % Simulation.MPF)
  if(Autoadjust.maxPartNum):
    Simulation.MaxPartNum = int(RestartParameter['Part-MaxParticleNumber'.lower()])
    print("updated MaxPartNum to %i" % Simulation.MaxPartNum)
  else:
    print("No autoadjust for MaxPartNum defined, use it from Shocktube.ini: %i" % Simulation.MaxPartNum)
  if(not RestartStanding):
    Iteration.CurrentVelo = float(RestartParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init1-VeloIC').lower()])
    print("set Iteration.CurrentVelo to %f" % Iteration.CurrentVelo)
  else:
    Iteration.VIn=float(RestartParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init1-VeloIC').lower()])
    print("set Iteration.VIn to %f" % Iteration.VIn)
    Iteration.VOut=float(RestartParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init2-VeloIC').lower()])
    print("set Iteration.VOut to %f" % Iteration.VOut)


  if(Iteration.StandingShock):
    SimulationParameter['Part-Boundary3-Condition'.lower()] = 'open'
    for iSpec in range(1,Simulation.nSpecies+1):
      if iSpec in Simulation.InflowSpecies:
        #Init
        SimulationParameter[('Part-Species'+str(iSpec)+'-nInits').lower()] = '2'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-SpaceIC').lower()] = 'cuboid'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-BaseVector1IC').lower()] = '(/0.,1.,0./)'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-BaseVector2IC').lower()] = '(/0.,0.,1./)'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-VeloVecIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-VeloVecIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-velocityDistribution').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-velocityDistribution').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-PartDensity').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-PartDensity').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-MWTemperatureIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-MWTemperatureIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-TempRot').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-TempRot').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-TempVib').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-TempVib').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-TempElec').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init2-TempElec').lower()]
        #Surface Flux
        SimulationParameter[('Part-Species'+str(iSpec)+'-nSurfaceFluxBCs').lower()] = '2'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-BC').lower()] = '3'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-velocityDistribution').lower()] = 'maxwell_lpn'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-VeloVecIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-VeloVecIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-PartDensity').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-PartDensity').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-MWTemperatureIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-MWTemperatureIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempRot').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempRot').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempVib').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempVib').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempElec').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempElec').lower()]
      else:
        #Init
        SimulationParameter[('Part-Species'+str(iSpec)+'-nInits').lower()] = '1'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-SpaceIC').lower()] = 'cuboid'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BaseVector1IC').lower()] = '(/0.,1.,0./)'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BaseVector2IC').lower()] = '(/0.,0.,1./)'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloVecIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-VeloVecIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-velocityDistribution').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-velocityDistribution').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-MWTemperatureIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-MWTemperatureIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempRot').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-TempRot').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempVib').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-TempVib').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempElec').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Init1-TempElec').lower()]
        #Surface Flux
        SimulationParameter[('Part-Species'+str(iSpec)+'-nSurfaceFluxBCs').lower()] = '1'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-BC').lower()] = '3'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-velocityDistribution').lower()] = 'maxwell_lpn'
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloVecIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloVecIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-PartDensity').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-PartDensity').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-MWTemperatureIC').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-MWTemperatureIC').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempRot').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempRot').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempVib').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempVib').lower()]
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempElec').lower()] = RestartParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempElec').lower()]


    print("")

def InitPiclasParameters():
  """
  Calculate the parameters for the next piclas simulation
  """
  global SimulationParameter, Simulation, Mesh, Iteration

  # Calculate Init area
  if not Iteration.StandingShock:

    SimulationParameter['TEnd'.lower()] = '%e' % Simulation.tend

    SimulationParameter['Part-FIBGMdeltas'.lower()] = '(/%e,1,1/)' % Mesh.dxElem
    if Simulation.MaxPartNum > 0:
      SimulationParameter['Part-MaxParticleNumber'.lower()] = '%i' % Simulation.MaxPartNum


    if Iteration.StandingShock==False:
      SimulationParameter['Part-IterationForMacroVal'.lower()] = '%i' % max(int(Simulation.tend/Simulation.dt/Mesh.nElems()*Simulation.ElemsPerWrite),1)
    SimulationParameter['Particles-ManualTimeStep'.lower()] = '%e' % Simulation.dt
    for iSpec in range(1,Simulation.nSpecies+1):
      SimulationParameter[('Part-Species'+str(iSpec)+'-MacroParticleFactor').lower()] = '%e' % (Simulation.MPF)
      if iSpec in Simulation.InflowSpecies:
        # Set Init
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BasePointIC').lower()] = '(/-%f,-0.5,-0.5/)' % (Mesh.Length)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-CuboidHeightIC').lower()] = '%f' % (Mesh.Length)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloIC').lower()] = '%f' % Iteration.CurrentVelo
        # Set Surfaceflux
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloIC').lower()] = '%f' % Iteration.CurrentVelo


  else:
    SimulationParameter['TEnd'.lower()] = '%e' % Simulation.tend

    SimulationParameter['Part-FIBGMdeltas'.lower()] = '(/%e,1,1/)' % Mesh.dxElem
    if Simulation.MaxPartNum > 0:
      SimulationParameter['Part-MaxParticleNumber'.lower()] = '%i' % Simulation.MaxPartNum

    SimulationParameter['Part-IterationForMacroVal'.lower()] = '%i' % max(int(Simulation.tend/Simulation.dt/Mesh.nElems()*Simulation.ElemsPerWrite),1)
    SimulationParameter['Particles-ManualTimeStep'.lower()] = '%e' % Simulation.dt
    for iSpec in range(1,Simulation.nSpecies+1):
      SimulationParameter[('Part-Species'+str(iSpec)+'-MacroParticleFactor').lower()] = '%e' % (Simulation.MPF)
      if iSpec in Simulation.InflowSpecies:
        # Set Init
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BasePointIC').lower()] = '(/-%f,-0.5,-0.5/)' % Mesh.Length
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-CuboidHeightIC').lower()] = '%f' % (Mesh.Length*1./2.)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloIC').lower()] = '%f' % Iteration.VIn
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-BasePointIC').lower()] = '(/-%f,-0.5,-0.5/)' % (Mesh.Length*1./2.)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-CuboidHeightIC').lower()] = '%f' % (Mesh.Length*1./2.)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-VeloIC').lower()] = '%f' % Iteration.VOut
        # Set Surfaceflux
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloIC').lower()] = '%f' % Iteration.VIn
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-VeloIC').lower()] = '%f' % Iteration.VOut
      else:
        # Set Init
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BasePointIC').lower()] = '(/-%f,-0.5,-0.5/)' % (Mesh.Length*1./2.)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-CuboidHeightIC').lower()] = '%f' % (Mesh.Length*1./2.)
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloIC').lower()] = '%f' % Iteration.VOut
        # Set Surfaceflux
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloIC').lower()] = '%f' % Iteration.VOut

def WriteParameterIni():
  """
  Write the parameter.ini for piclas in the order of the parameters of "piclas --help"
  """
  global SimulationParameter, Path
  # if not os.path.exists("piclas.help") :

  #   f=open("piclas.help","w")
    # sp.Popen([os.path.expanduser("~/masterdev_ref_piclas/build/bin/piclas","--help")])
  p=sp.Popen([Path.piclas,"--help"],stdout=sp.PIPE)
  i=-1
  f=open("parameter.ini","w")
  output=[]
  # Read stdout of "piclas --help"
  for line in p.stdout:
    if sys.version_info[0] == 3:
      line = line.decode()
    i+=1
    if i<22:
      # Remove header
      continue
    if line.startswith(' ') :
      # Line with no parameter
      continue
    if line.startswith('!') :
      # New section in piclas --help
      if len(output)<3:
        # Add line to output
        output.append(line)
      elif len(output)==3:
        # Last section without setted parameters --> remove section header
        output=[]
        output.append(line)
      else:
        # Last section with parameter --> write section to file
        for ii in range(len(output)):
          f.write(output[ii])
        output=[]
        output.append(line)

      # f.write(line)
    else:
      # Line with parameter
      line = re.sub(r"\s+", "", line)
      (line,dumb)=line.split('=',1)
      if line.count('$') == 0:
        # No mutiple parameter
        if line.lower() in SimulationParameter.keys():
          # If in SimulationParameter dictionary add to output
          output.append(line + '=' + SimulationParameter[line.lower()] + '\n')
      else:
        # Multiple parameter (with [$])
          indexes=np.ones(line.count('$'),dtype=int)
          maxindex = max(Simulation.nSpecies,3)
          lastchange=0
          while indexes[0]<=maxindex:
            lineSpec=line
            for ii in range(len(indexes)):
              lineSpec = lineSpec[0:lineSpec.find('$')-1] + str(indexes[ii]) + lineSpec[lineSpec.find('$')+2:]
            # print(lineSpec)
            if lineSpec.lower() in SimulationParameter.keys():
              # f.write(lineSpec + '=' + SimulationParameter[lineSpec.lower()] + '\n')
              output.append(lineSpec + '=' + SimulationParameter[lineSpec.lower()] + '\n')
            indexes[-1]+=1
            lastchange=len(indexes)-1
            if indexes[-1]>maxindex:
              ii=1
              while indexes[-ii]>maxindex:
                if ii == len(indexes):
                  break
                else:
                  indexes[-ii] = 1
                  indexes[-(ii+1)]+=1
                  ii+=1
          # elif lastchange == 0:
          #   # output.append(line+'\n')
          #   break
          # else:
          #   indexes[lastchange-1]+=1
          #   for ii in range(lastchange,len(indexes)):
          #     indexes[ii]=1
          #   lastchange+=-1

  if len(output)>3:
    # Write last output if it is not only a header
    for ii in range(len(output)):
      f.write(output[ii])
  f.close()

def RunPiclas():
  """
  Run Piclas with the ini and DSMCSpecies file
  """
  global Path, Simulation
  RunCmd = []
  if Simulation.UseMPI:
    RunCmd.append('mpirun')
    RunCmd.append('-np')
    RunCmd.append(str(Simulation.nMPICores))
  RunCmd.append(Path.piclas)
  RunCmd.append('parameter.ini')
  if Path.DSMCSpecies != 'None':
    RunCmd.append(Path.DSMCSpecies)
  if Simulation.Redosimulation:
    StateFiles = sorted(glob.glob("*_State_*"))
    RunCmd.append(StateFiles[-1])
  else:
    StateFiles = sorted(glob.glob("*State_*"))
    for File in StateFiles:
      os.remove(File)
  RunString = ''
  for cmd in RunCmd:
    RunString += ' ' + cmd
  print('Run piclas with:%s' % RunString)
  p=sp.Popen(RunCmd,stdout=sp.PIPE)
  output=[]
  with open('piclas.out','w') as f:
    for line in p.stdout:
      if sys.version_info[0] == 3:
        line = line.decode()
      line = line[:-1]
      output.append(line)
      f.write(line+'\n')



  # output,_ = p.communicate()
  # with open('piclas.out','w') as f:
  #   f.write(output)
  # output=output.split('\n')
  if not output[-2].startswith(' PICLAS FINISHED!'):
    print('PICLAS TERMINATED WITH AN ERROR!!!')
    i=1
    for line in reversed(output):
      i+=1
      if line.startswith(' _____________________________________________________________________________'):
        break
    for ii in reversed(range(1,i)):
      print(output[-ii])
    raise Exception("An Error occurs during piclas execution. See above.")

def GenerateMeshFile():
  """
  Generate the Hopr.ini File and executes hopr
  """
  global SimulationParameter, Mesh, Path
  f=open("hopr.ini","w")
  f.write("! " + "="* 127 + " !\n")
  f.write("! OUTPUT\n")
  f.write("! " + "="* 127 + " !\n")
  f.write("  ProjectName  = %s\n" % SimulationParameter["ProjectName".lower()])
  f.write("  Debugvisu    = F\n")

  f.write("! " + "="* 127 + " !\n")
  f.write("! MESH\n")
  f.write("! " + "="* 127 + " !\n")
  f.write("  Mode         =1\n")
  f.write("  nZones       =1\n")
  f.write("  jacobianTolerance=1e-30\n")
  f.write("  Corner=(/")
  f.write("-%f,-%f,-%f,,"   % (Mesh.Length,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write(" %f,-%f,-%f,,"   % (0,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write(" %f, %f,-%f,,"   % (0,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write("-%f, %f,-%f,,"   % (Mesh.Length,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write("-%f,-%f, %f,,"   % (Mesh.Length,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write(" %f,-%f, %f,,"   % (0,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write(" %f, %f, %f,,"   % (0,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write("-%f, %f, %f/)\n" % (Mesh.Length,Mesh.HalfWidth,Mesh.HalfWidth))
  f.write("  nElems       =(/%i,1,1/)\n" % Mesh.nElems())
  f.write("  BCIndex      =(/1,1,3,1,2,1/)\n")
  f.write("  elemtype     =108\n")

  f.write("! " + "="* 127 + " !\n")
  f.write("! BOUNDARY CONDITIONS\n")
  f.write("! " + "="* 127 + " !\n")
  f.write("    BoundaryName=BC_Sym\n")
  f.write("    BoundaryType=(/4,0,0,0/)\n")
  f.write("    BoundaryName=BC_Inflow\n")
  f.write("    BoundaryType=(/3,0,0,0/)\n")
  f.write("    BoundaryName=BC_Outflow\n")
  f.write("    BoundaryType=(/3,0,0,0/)\n")

  f.close()

  DEVNULL = open(os.devnull, 'w')
  print('Run hopr with: %s hopr.ini' % Path.hopr)
  sp.call([Path.hopr,"hopr.ini"],stdout=DEVNULL)

def ReadMeshFile(MeshfileName):
  """
  Create the by x position ascending sorted list of ElemIDs out of the Mesh File
  """
  global ElemIDsorted
  global Barycenters, Mesh
  try:
    MeshFile = h5py.File(MeshfileName, 'r')
  except:
    raise Exception("generated mesh file '%s' not found." % MeshfileName)
  Barycenters=MeshFile['ElemBarycenters'][:,0].copy()
  ElemIDsorted=np.arange(Mesh.nElems())

  n = len(ElemIDsorted)
  swapped=True
  while swapped:
    swapped = False
    for i in range(n-1):
      if Barycenters[i] > Barycenters[i+1]:
        ElemIDTemp = ElemIDsorted[i]
        ElemIDsorted[i] = ElemIDsorted[i+1]
        ElemIDsorted[i+1] = ElemIDTemp
        BarycentersTemp = Barycenters[i]
        Barycenters[i] = Barycenters[i+1]
        Barycenters[i+1] = BarycentersTemp
        swapped = True
    n = n-1

def SortArrayToPosition(Array):
  """
  Sort the incoming array to ascending positions
  """
  global Mesh
  if len(Array)!=Mesh.nElems():
    sys.exit('The to be sorted array contains more or less entrys than nElems')
  ArraySorted=np.zeros(Mesh.nElems())
  for i in range(Mesh.nElems()):
    ArraySorted[i] = Array[ElemIDsorted[i]].copy()
  return ArraySorted.copy()

def DetermineShockSpeed():
  """
  Calculate Shock Speed based on the time and xShock data
  """
  global Solution
  # Filter valid values
  nValues = np.sum(Solution.xShockValid)
  time = np.zeros(nValues,dtype=float)
  xShock = np.zeros(nValues,dtype=float)
  Curvatures = np.zeros(nValues,dtype=float)
  AveragesList = np.zeros([7,nValues],dtype=float)
  ii = 0
  with open('rawdata.dat','w') as f:
    f.write('time[s]\txShock[m]\tValid\n')
    for i in range(len(Solution.times)):
      f.write('%e\t%e\t%r\n' % (Solution.times[i],Solution.xShock[i],Solution.xShockValid[i]))
      if Solution.xShockValid[i]:
        time[ii]=Solution.times[i]
        xShock[ii]=Solution.xShock[i]
        ii+=1
  Solution.times = time
  Solution.xShock = xShock

  print('Number of Valid Values: %i'%nValues)
  # Regression Loop
  # 1. Determine Curvature of the given Data
  # 2. Remove first Data Pair and calculate Curvature again
  # 3. Continue untill curvature changes its sign
  # 4. Linear regression of the remaining Data
  CurvatureFirst=0.0
  CurvatureNew=0.0
  i=-1
  RegressionCoeff=np.zeros([7,nValues],dtype=float)
  Averages=np.zeros(7,dtype=float)
  while CurvatureNew*CurvatureFirst>0 or i<2:
    i+=1
    if i==0:
      # Create Coefficinet Vecotors and their Averages
      for ii in range(nValues):
        RegressionCoeff[0,ii]=time[ii]
        RegressionCoeff[1,ii]=xShock[ii]
        RegressionCoeff[2,ii]=time[ii]*xShock[ii]
        RegressionCoeff[3,ii]=time[ii]**2
        RegressionCoeff[4,ii]=time[ii]**3
        RegressionCoeff[5,ii]=time[ii]**4
        RegressionCoeff[6,ii]=time[ii]**2*xShock[ii]
      for ii in range(7):
        if nValues==0:
          Solution.Shockspeed = 0
          return
        else:
          Averages[ii]=np.average(RegressionCoeff[ii,:])
    else:
      # Adjust Averages by removing first value
      for ii in range(7):
        # Averages[ii]=((k+1)*Averages[ii]-RegressionCoeff[ii,i-1])/k
        # Averages[ii]=np.sum(RegressionCoeff[ii,i:-1])/(nValues-i)
        Averages[ii]=np.average(RegressionCoeff[ii,i:])

    CurvatureNew=((Averages[6]-Averages[1]*Averages[3])*(Averages[3]-Averages[0]**2) \
                -(Averages[2]-Averages[0]*Averages[1])*(Averages[4]-Averages[0]*Averages[3])) \
                /((Averages[5]-Averages[3]**2)*(Averages[3]-Averages[0]**2) \
                -(Averages[4]-Averages[0]*Averages[3])**2)
    Curvatures[i] = CurvatureNew
    AveragesList[:,i] = Averages[:]
    if i==0:
      CurvatureFirst = CurvatureNew
  Simulation.UsedValues = nValues-i+1
  print('Number of Used Values: %i'%Simulation.UsedValues)
  if Simulation.UsedValues<0.1*nValues:
    Simulation.UsedValues=int(0.1*nValues)
    for ii in range(7):
      Averages[ii]=np.average(RegressionCoeff[ii,-Simulation.UsedValues:])
      # Averages[ii]=np.sum(RegressionCoeff[ii,-Simulation.UsedValues:])/nValues
    print('Too low number of valid values in shock speed determination. Use the last 10%% of all values: %i'%Simulation.UsedValues)
    Solution.TooFewValues=True

  # Print to file
  with open('ShockPositions.txt','w') as f:
    f.write('nTotalValues:,%i,nUsedValues:,%i\n' % (nValues,Simulation.UsedValues))
    f.write('t,xShock,Used,t,x,tx,t^2,t^3,t^4,xt^2,a_t,a_x,a_tx,a_t^2,a_t^3,a_t^4,a_xt^2,Curvature\n')
    for i in range(nValues):
      if nValues-Simulation.UsedValues>=i:
        f.write('%e,%e,%r,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n' % (time[i],xShock[i],False,RegressionCoeff[0,i],RegressionCoeff[1,i],RegressionCoeff[2,i],RegressionCoeff[3,i],RegressionCoeff[4,i],RegressionCoeff[5,i],RegressionCoeff[6,i],AveragesList[0,i],AveragesList[1,i],AveragesList[2,i],AveragesList[3,i],AveragesList[4,i],AveragesList[5,i],AveragesList[6,i],Curvatures[i]))
      else:
        f.write('%e,%e,%r,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n' % (time[i],xShock[i],True,RegressionCoeff[0,i],RegressionCoeff[1,i],RegressionCoeff[2,i],RegressionCoeff[3,i],RegressionCoeff[4,i],RegressionCoeff[5,i],RegressionCoeff[6,i],AveragesList[0,i],AveragesList[1,i],AveragesList[2,i],AveragesList[3,i],AveragesList[4,i],AveragesList[5,i],AveragesList[6,i],Curvatures[i]))

  Solution.Shockspeed = (Averages[2]-Averages[0]*Averages[1])/(Averages[3]-Averages[0]**2)
  # if Solution.Shockspeed<0.0:
  #   Solution.Shockspeed = 0.0

def ReadDSMCHOData():
  global Autoadjust, Solution, Simulation
  DSMCHOFiles = sorted(glob.glob("*SMCState*"))
  nFiles = len(DSMCHOFiles)
  times = np.zeros(nFiles,dtype=float)
  xShock = np.zeros(nFiles,dtype=float)
  xShockvalid = np.zeros(nFiles,dtype=bool)
  iTotalTemp = -1
  iTotalNumDens = -1
  iTotalSimPartNum = -1
  iMCSoverMFP = -1
  iSpecTemp = -1*np.ones(Simulation.nSpecies+1, dtype=int)
  iSpecVelo = -1*np.ones(Simulation.nSpecies+1, dtype=int)
  iMaxCollProb = -1
  p2=0.0
  for iFile in range(-1,nFiles):
    sys.stdout.write("\rread and processed %i output files of a total of %i files" % (iFile+1,nFiles))
    sys.stdout.flush()
    if iFile==-1:
      iFile=int(nFiles/2)
    elif iFile==int(nFiles/2):
      continue
    # Read NumDens and Temperature Data and sort them by position
    DSMCHOFile = h5py.File(DSMCHOFiles[iFile], 'r')
    times[iFile] = DSMCHOFile.attrs['Time']
    if iFile==int(nFiles/2):
      for iName in range(len(DSMCHOFile.attrs['VarNamesAdd'])):
        VarName = DSMCHOFile.attrs['VarNamesAdd'][iName]
        if sys.version_info[0] == 3:
          VarName = VarName.decode()
        # print(iName,VarName)
        if VarName=='Total_NumberDensity':
          iTotalNumDens=iName
        if VarName=='Total_TempTransMean':
          iTotalTemp=iName
        if Autoadjust.MPF:
          if VarName=='DSMC_MCS_over_MFP':
            iMCSoverMFP=iName
          if VarName=='Total_SimPartNum':
            iTotalSimPartNum=iName
        if Autoadjust.dt:
          if VarName=='DSMC_MaxCollProb':
            iMaxCollProb = iName
          if Simulation.nSpecies ==1:
            if VarName=='Total_TempTransX':
              iSpecTemp[0] = iName
            if VarName=='Total_VeloX':
              iSpecVelo[0] = iName
          else:
            for iSpec in range(1,Simulation.nSpecies+1):
              if VarName=='Spec%03i_TempTransX' % (iSpec):
                iSpecTemp[iSpec] = iName
              if VarName=='Spec%03i_VeloX' % (iSpec):
                iSpecVelo[iSpec] = iName

      if iTotalTemp==-1 and iTotalNumDens==-1:
        sys.exit("Position of Total_NumberDensity and Total_TempTransMean in DSMCState-File cannot be determined")
      if iTotalNumDens==-1:
        sys.exit("Position of Total_NumberDensity in DSMCState-File cannot be determined")
      if iTotalTemp==-1:
        sys.exit("Position of Total_TempTransMean in DSMCState-File cannot be determined")
      if Autoadjust.MPF:
        if iTotalSimPartNum==-1 or iMCSoverMFP==-1:
          sys.exit("Position of 'Total_SimPartNum' or 'DSMC_MCS_over_MFP' in DSMCState-File cannot be determined")


    TotalTemp=SortArrayToPosition(DSMCHOFile['ElemData'][:,iTotalTemp].copy())
    TotalNumDens=SortArrayToPosition(DSMCHOFile['ElemData'][:,iTotalNumDens].copy())

    # Calculate Position of the shock front
    TotalPress = TotalTemp*TotalNumDens*kb
    if p2==0:
      p2=TotalPress[-1]
    RelaPress = (TotalPress-Simulation.InflowPressure)/(p2-Simulation.InflowPressure)
    # Filter RealPress
    RelaPressFiltered = np.zeros(len(RelaPress),dtype=float)
    intervall=int(0.05*Mesh.nElems())
    for ii in range(Mesh.nElems()):
      iStart=max(0,ii-intervall)
      iEnd  =min(Mesh.nElems()-1,ii+intervall)+1
      for iii in range(iStart,iEnd):
        RelaPressFiltered[ii]+=RelaPress[iii]
      RelaPressFiltered[ii]/=iEnd-iStart
    # Find Elements at the transition from RelaPressFiltered<0.5 to >0.5
    iShock=0
    Found=True
    while RelaPressFiltered[iShock+1]<0.5:
      if iShock<Mesh.nElems()-2:
        iShock+=1
      else:
        Found=False
        break
    if Found and iShock>10:
      xShock[iFile]=Barycenters[iShock] + (Barycenters[iShock+1]-Barycenters[iShock]) \
        / (RelaPressFiltered[iShock+1]-RelaPressFiltered[iShock]) * (0.5-RelaPressFiltered[iShock])
      if xShock[iFile]>Barycenters[intervall] and xShock[iFile]<Barycenters[-intervall]:
        xShockvalid[iFile]=True

    # Autoadjust
    # MPF
    if Autoadjust.MPF:
      MCSoverMFP=DSMCHOFile['ElemData'][:,iMCSoverMFP].copy()
      TotalSimPartNum=DSMCHOFile['ElemData'][:,iTotalSimPartNum].copy()
      Solution.MCSoverMFPmax=max(max(MCSoverMFP),Solution.MCSoverMFPmax)
      if min(TotalSimPartNum)!=0:
        Solution.TotalSimPartNummin=min(min(TotalSimPartNum),Solution.TotalSimPartNummin)
    if Autoadjust.dt:
      Solution.MaxCollProb=max(max(DSMCHOFile['ElemData'][:,iMaxCollProb].copy()),Solution.MaxCollProb)
      for iSpec in range(1,Simulation.nSpecies+1):
        # print(iSpecVelo[iSpec])
        # print(max(DSMCHOFile['ElemData'][:,iSpecVelo[iSpec]].copy()))
        # print(Simulation.VeloTempEps)
        # print(max(DSMCHOFile['ElemData'][:,iSpecTemp[iSpec]].copy()))
        # print(pi)
        # print(max(np.sqrt(2*kb*DSMCHOFile['ElemData'][:,iSpecTemp[iSpec]].copy()/(pi*Simulation.SpeciesMass[iSpec]))))
        # print(Solution.MaxVelo)
        # print("")
        if not iSpec == Simulation.ElectronSpecies:
          # Negelecting Electrons
          Solution.MaxVelo=max(max(DSMCHOFile['ElemData'][:,iSpecVelo[iSpec]].copy()+Simulation.VeloTempEps*np.sqrt(2*kb*DSMCHOFile['ElemData'][:,iSpecTemp[iSpec]].copy()/(pi*Simulation.SpeciesMass[iSpec]))),Solution.MaxVelo)

  Solution.times = times.copy()
  Solution.xShock = xShock.copy()
  Solution.xShockValid = xShockvalid.copy()

def UpdateVariables():
  global Iteration,Simulation, Solution, Autoadjust, SimulationParameter
  if Solution.Shockspeed == 0:
    # Shockspeed cannot be determined
    Simulation.Redosimulation = True
    if Autoadjust.tend:
      Simulation.tend *=2
      print('No Shock can be determined. Extend tend to %e.' % Simulation.tend)
      return False
    else:
      raise Exception("Shockspeed cannot be determined")
  else:
    if Iteration.StandingShock == False:
      Simulation.Redosimulation = False
      Iteration.LastLagShockSpeed = Iteration.CurrentLagShockSpeed
      Iteration.CurrentLagShockSpeed = Iteration.CurrentVelo-Solution.Shockspeed
      if Iteration.CurrentLagShockSpeed<0.0:
        Simulation.Redosimulation = True
        if Autoadjust.tend:
          Simulation.tend *=2
          print('No Shock can be determined. Extend tend to %e.' % Simulation.tend)
          return False
        else:
          raise Exception("Shockspeed cannot be determined")

      print('Shockspeed: %f'%Iteration.CurrentLagShockSpeed)

      if abs(Iteration.CurrentLagShockSpeed-Iteration.TargetVelo)<= Iteration.TargetDev:
        # End of iteration reached
        return True
      if Iteration.LastVelo != None and Iteration.LastLagShockSpeed != None:
        # Update ResponseRatio
        Iteration.ResponseRatio = (Iteration.CurrentVelo-Iteration.LastVelo)/(Iteration.CurrentLagShockSpeed-Iteration.LastLagShockSpeed)
        # Iteration.ResponseRatio = Iteration.CurrentVelo * Iteration.LastLagShockSpeed / (Iteration.LastVelo * Iteration.CurrentLagShockSpeed)
        print('New ResponseRatio: %e' % Iteration.ResponseRatio)

      # Adjust CurrentVelo
      Iteration.LastVelo = Iteration.CurrentVelo
      Iteration.CurrentVelo = Iteration.ResponseRatio*(Iteration.TargetVelo-Iteration.CurrentLagShockSpeed)+Iteration.CurrentVelo
      # Iteration.CurrentVelo = Iteration.LastVelo * Iteration.ResponseRatio * Iteration.TargetVelo / Iteration.CurrentLagShockSpeed
      print('New Velocity: %f' % Iteration.CurrentVelo)

      # Adjust tend
      if Autoadjust.tend:
        xend = Barycenters[int(0.25*Mesh.nElems())]
        Simulation.tend = (xend-Solution.xShock[-1]+Solution.times[-1]*Solution.Shockspeed)/Solution.Shockspeed
        Simulation.tend *= Iteration.CurrentVelo/Iteration.LastVelo
        print('Adjust tend to %e' % Simulation.tend)
      if Autoadjust.MPF:
        if(('Particles-DSMC-UseNearestNeighbour'.lower()) in SimulationParameter):
          if(SimulationParameter[('Particles-DSMC-UseNearestNeighbour').lower()] in ['true', 't','1']):
            Simulation.MPF = min(Simulation.MPF * Solution.TotalSimPartNummin/Simulation.PartperElemmin,Simulation.MPF * Simulation.MCSoverMFPmax/Solution.MCSoverMFPmax)
          else:
            Simulation.MPF = Simulation.MPF * Solution.TotalSimPartNummin/Simulation.PartperElemmin
        else:
          Simulation.MPF = Simulation.MPF * Solution.TotalSimPartNummin/Simulation.PartperElemmin
        print('Adjust MPF to %e' % Simulation.MPF)
      if Autoadjust.dt:
        dt1=Simulation.tend/Mesh.nElems()
        dt2=Mesh.dxElem/Solution.MaxVelo*Iteration.CurrentVelo/Iteration.LastVelo
        dt3=Simulation.dt*Simulation.MaxCollProbMax/Solution.MaxCollProb
        Simulation.dt=min(dt1,dt2,dt3)
        # Simulation.dt = min(Simulation.tend/Mesh.nElems(),Mesh.dxElem/Solution.MaxVelo*Iteration.CurrentVelo/Iteration.LastVelo,Simulation.dt*Simulation.MaxCollProbMax/Solution.MaxCollProb)
        print('Adjust dt to %e' % Simulation.dt)
      if Autoadjust.Mesh:
        if Simulation.UsedValues<=0.25*Mesh.nElems():
          Mesh.Length*=2
          Simulation.tend*=2
          print('Adjust Length to %f and therfore doubled tend to %e' % (Mesh.Length,Simulation.tend))
    else:
      IntegrateSolution2()
      Finished = True
      Iteration.VIn = Iteration.VIn - Solution.Shockspeed
      Iteration.VOut = Iteration.VOut - Solution.Shockspeed
      Iteration.LastVelo = Iteration.CurrentVelo
      Iteration.CurrentVelo=Iteration.VIn
      print("Shockspeed loc: %f" % (Solution.Shockspeed))
      print("Shockspeed: %f" % (Iteration.VIn - Solution.Shockspeed))
      print("Mesh.dxElem: %f" % Mesh.dxElem)
      print("Solution.Shockspeed*Simulation.tend: %f" % (Solution.Shockspeed*Simulation.tend))
      if abs(Solution.Shockspeed*Simulation.tend) >= Mesh.dxElem:
        Finished = False
      print("Convection of Shock front/xElem: %e" % (Solution.Shockspeed*Simulation.tend))
      if Autoadjust.MPF:
        std = np.std(Solution.xShock[len(Solution.xShock)-Simulation.UsedValues:])
        print("Standarddeviation/xElem: %e" % (std/Mesh.dxElem))
        if std >= Mesh.dxElem:
          Finished = False
        if(('Particles-DSMC-UseNearestNeighbour'.lower()) in SimulationParameter):
          if(SimulationParameter[('Particles-DSMC-UseNearestNeighbour').lower()] in ['true', 't']):
            MPFratio = min((Mesh.dxElem/(2*std))**2*Simulation.MPF,Simulation.MPF * Solution.TotalSimPartNummin/Simulation.PartperElemmin,Simulation.MPF * Simulation.MCSoverMFPmax/Solution.MCSoverMFPmax)/(Simulation.MPF*2)
          else:
            MPFratio = min((Mesh.dxElem/(2*std))**2*Simulation.MPF,Simulation.MPF * Solution.TotalSimPartNummin/Simulation.PartperElemmin)/(Simulation.MPF*2)
        else:
          MPFratio = min((Mesh.dxElem/(2*std))**2*Simulation.MPF,Simulation.MPF * Solution.TotalSimPartNummin/Simulation.PartperElemmin)/(Simulation.MPF*2)
        Simulation.MPF *= MPFratio
        print('Adjust MPF to %e' % Simulation.MPF)
        if Autoadjust.tend:
          tend1 = Solution.times[0] + (Simulation.tend-Solution.times[0])*MPFratio
          tend2 = 5./4.* Mesh.Length/(Iteration.VOut)
          Simulation.tend = max(tend1,tend2)
        # if Autoadjust.maxPartNum:
        #   Simulation.MaxPartNum /= MPFratio
        print('Adjust MaxPartNum to %i' % Simulation.MaxPartNum)
      if Autoadjust.dt:
        dt2=Mesh.dxElem/Solution.MaxVelo*Iteration.CurrentVelo/Iteration.LastVelo
        print(Mesh.dxElem,Solution.MaxVelo,Iteration.CurrentVelo,Iteration.LastVelo)
        dt3=Simulation.dt*Simulation.MaxCollProbMax/Solution.MaxCollProb
        print(dt2,dt3)
        dtratio = min(dt2,dt3)/(Simulation.dt*2)/2
        Simulation.dt *= dtratio
        if Solution.MaxCollProb>1:
          Finished= False
        if Autoadjust.tend:
          Simulation.tend = Solution.times[0] + (Simulation.tend-Solution.times[0])*dtratio
        print('Adjust dt to %e' % Simulation.dt)
      if Autoadjust.tend:
        Medians=[]
        with h5py.File("Integrated.h5", 'r') as IntegratedFile:
          ElemData = IntegratedFile['ElemData']
          for i in range(np.shape(ElemData)[1]):
            if sys.version_info[0] == 2:
              if(("VeloY" in IntegratedFile.attrs["VarNamesAdd"][i]) or ("VeloZ" in IntegratedFile.attrs["VarNamesAdd"][i]) or\
                ("TempTransY" in IntegratedFile.attrs["VarNamesAdd"][i]) or ("TempTransZ" in IntegratedFile.attrs["VarNamesAdd"][i]) or\
                ("DSMC" in IntegratedFile.attrs["VarNamesAdd"][i]) or ("SimPartNum" in IntegratedFile.attrs["VarNamesAdd"][i])):
                continue
            elif sys.version_info[0] == 3:
              if(("VeloY" in IntegratedFile.attrs["VarNamesAdd"][i].decode()) or ("VeloZ" in IntegratedFile.attrs["VarNamesAdd"][i].decode()) or\
                ("TempTransY" in IntegratedFile.attrs["VarNamesAdd"][i].decode()) or ("TempTransZ" in IntegratedFile.attrs["VarNamesAdd"][i].decode()) or\
                ("DSMC" in IntegratedFile.attrs["VarNamesAdd"][i].decode()) or ("SimPartNum" in IntegratedFile.attrs["VarNamesAdd"][i].decode())):
                continue
            Diff = []
            Data = SortArrayToPosition(ElemData[:,i].copy())
            Data = abs(Data)
            if max(Data) == 0:
              continue
            Data /= max(Data)
            for ii in range(len(Data)-1):
              Diff.append(abs(Data[ii]-Data[ii+1]))
            Medians.append(np.median(Diff))
        print(Medians)
        Median = np.average(Medians)
        print('Median: %f' %Median)
        if Median>Simulation.Median:
          Finished = False
        Tend1 = Solution.times[0] + (Simulation.tend-Solution.times[0]) * (Median*2/(Simulation.Median))**2
        Tend2 = Mesh.Length/(2*Iteration.VOut)
        if Solution.TooFewValues:
          Tend3 = Simulation.tend*2
          Finished = False
        else:
          Tend3 = 0
        Simulation.tend = max(Tend1,Tend2,Tend3)
        print('Adjust tend to %e' % Simulation.tend)
      if Finished:
        # End of iteration reached
        return True










  return False

def PrepareSimulation():
  """
  Prepare next simulation run. Delete old data or create new directory.
  """
  global Path, Iteration, Simulation, Autoadjust, Solution
  if not Path.KeepSimulationData:
    DSMCHOFiles = glob.glob("*SMCState*")
    for File in DSMCHOFiles:
      os.remove(File)
    StateFiles = glob.glob("*_State_*")
    if Simulation.Redosimulation:
      StateFiles = sorted(StateFiles)
      for File in StateFiles[0:-1]:
        os.remove(File)
    else:
      for File in StateFiles:
        os.remove(File)
  else:
    Path.LastSimDir = Path.CurrentSimDir
    if Iteration.StandingShock==False:
      os.mkdir('%s/iter_%02i' % (Path.MasterDir,Iteration.iter))
      Path.CurrentSimDir = '%s/iter_%02i' % (Path.MasterDir,Iteration.iter)
    else:
      os.mkdir('%s/standing_iter_%02i' % (Path.MasterDir,Iteration.iter))
      Path.CurrentSimDir = '%s/standing_iter_%02i' % (Path.MasterDir,Iteration.iter)
    if Iteration.iter == 0:
      Path.LastSimDir = Path.MasterDir
    shutil.copyfile('%s/%s' % (Path.MasterDir,Path.DSMCSpecies),'%s/%s' % (Path.CurrentSimDir,Path.DSMCSpecies))
    if Simulation.Redosimulation:
      StateFiles = sorted(glob.glob("*_State_*"))
      StateFile = StateFiles[-1]
      shutil.copyfile('%s/%s' % (Path.LastSimDir,StateFile),'%s/%s' % (Path.CurrentSimDir,StateFile))
    if not Autoadjust.Mesh:
      MeshFile=SimulationParameter["ProjectName".lower()] + '_mesh.h5'
      shutil.copyfile('%s/%s' % (Path.MasterDir,MeshFile),'%s/%s' % (Path.CurrentSimDir,MeshFile))
    if len(Path.AdditionalCopyFiles)>0:
      for File in Path.AdditionalCopyFiles:
        shutil.copyfile('%s/%s' % (Path.MasterDir,File),'%s/%s' % (Path.CurrentSimDir,File))
    os.chdir(Path.CurrentSimDir)

  Solution.MaxCollProb=0
  Solution.MaxVelo=0
  Solution.MCSoverMFPmax=0
  # Solution.TotalSimPartNummin=sys.maxint-1
  Solution.TotalSimPartNummin=sys.maxsize
  Solution.TooFewValues=False

  # Output Simparameter:

  print('tend:                     %e' % Simulation.tend)
  print('dt:                       %e' % Simulation.dt)
  print('Velocity:                 %f' % Iteration.CurrentVelo)
  print('iterations per DSMCState: %i' % max(int(Simulation.tend/Simulation.dt/Mesh.nElems()*Simulation.ElemsPerWrite),1))
  print('total iterations:         %i' % (Simulation.tend/Simulation.dt))
  print('length of domain:         %e' % Mesh.Length)
  print('MPF:                      %e' % Simulation.MPF)
  if(Simulation.MaxPartNum>0):
    print('MaxPartNum:               %i' % Simulation.MaxPartNum)

def GetEquilibriumCondition2():
  global Solution, Barycenters, Mesh, SimulationParameter
  # print(glob.glob("*SMCHO*"))
  if(('Particles-DSMC-CollisMode'.lower()) in SimulationParameter):
    if(SimulationParameter['Particles-DSMC-CollisMode'.lower()]=='3'):
      Chem=True
    else:
      Chem=False
  else:
    Chem=False

  DSMCHOFiles = sorted( glob.glob("*SMCState*"))
  iFile = len(DSMCHOFiles)-1
  # Skip unvalid files
  while(not Solution.xShockValid[iFile]):
    iFile -= 1
  # Generate VarNames
  VarNames = []
  VarNames.append('Total_NumberDensity')
  VarNames.append('Total_VeloX')
  VarNames.append('Total_TempTransMean')
  if(Simulation.nSpecies!=1):
    for iSpec in range(1,Simulation.nSpecies+1):
      VarNames.append('Spec%03i_NumberDensity'%iSpec)
  # Find Shock Elem
  Data = []
  DSMCHOFile = h5py.File(DSMCHOFiles[iFile], 'r')
  for VarName in VarNames:
    iData=-1
    for iName in range(len(DSMCHOFile.attrs['VarNamesAdd'])):
      if sys.version_info[0] == 2:
        if DSMCHOFile.attrs['VarNamesAdd'][iName]==VarName:
          iData = iName
      elif sys.version_info[0] == 3:
        if DSMCHOFile.attrs['VarNamesAdd'][iName].decode()==VarName:
          iData = iName
    if(iData==-1):
      raise Exception('Data of "%s" cannot be found in file %s'%(VarName,DSMCHOFiles[iFile]))
    Data.append(SortArrayToPosition(DSMCHOFile['ElemData'][:,iData].copy()))

  # Determine iElemShock by max Transtemp
  iElem = 0
  maxTransTemp = np.max(Data[2])
  while((Data[2])[iElem]!=maxTransTemp):
    iElem+=1
  iElemShock = iElem

  iElem = -1
  iElemNumDens = -1
  # while (iElemNumDens == -1):
  #   iElem += 1
  #   if Data[0][iElem]>Data[0][-1]:
  #     iElemNumDens = iElem
  meanNumDens = (np.max(Data[0])+np.min(Data[0]))*0.5
  while (iElemNumDens == -1):
    iElem += 1
    if Data[0][iElem]>meanNumDens:
      iElemNumDens = iElem

  iElem = iElemShock-1
  iElemTransTemp = -1
  while (iElemTransTemp == -1):
    iElem += 1
    if Data[2][iElem]>Data[2][-1]:
      iElemTransTemp = iElem

  iElemShock = max(iElemTransTemp,iElemNumDens)
  print('iElemShock: %i =MAX(%i,%i)' % (iElemShock,iElemTransTemp,iElemNumDens))


  Values = []

  # DetermineVanderMondeMatrix
  PolynomOrder = 4
  Vandermonde = np.ones([Mesh.nElems()-iElemShock,PolynomOrder+1],dtype=float)

  for ix in range(iElemShock,Mesh.nElems()):
    for n in range(PolynomOrder+1):
      Vandermonde[ix-iElemShock,n] = Barycenters[ix]**n

  Vandermonde = np.matmul(  np.linalg.inv( np.matmul( np.transpose(Vandermonde) ,Vandermonde) )  ,np.transpose(Vandermonde))


  Output = []
  Output.append('x,iElem,')
  for iElem in range(0,Mesh.nElems()):
    Output.append('%f,%i,'%(Barycenters[iElem],iElem))

  listcoeffs=[]
  listregtype=[]
  x0real = []
  if(Chem):
    upperbound = len(VarNames)
  else:
    # Omit species specific number densities
    upperbound = len(VarNames)-Simulation.nSpecies
  for iData in range(upperbound):
    VarName = VarNames[iData]
    DataVar = Data[iData]

    print('Calculating post-shock condition of %s' % VarName)
    # Determine regression function type
    if(np.max(DataVar) == np.min(DataVar)):
      listregtype.append('const')
      listcoeffs.append(np.max(DataVar))
    elif (list(Data[iData][iElemShock:]).count(0)<0.05*(Mesh.nElems()-iElemShock)):
      listregtype.append('poly')
      a = np.matmul(Vandermonde, Data[iData][iElemShock:])
      # print("a: ",a)
      listcoeffs.append(a)
      # Build derevative and find roots of it (Extremwerte)
      a_ = np.zeros(PolynomOrder,dtype=float)
      for i in range(1,PolynomOrder+1):
        a_[i-1]=a[i]*i
      # print("a_:",a_)
      x0=np.roots(np.flip(a_,0))
      # print("x0:",x0)
      for i in range(len(x0)):
        x = x0[i]
        # print("x",x)
        # print("x.imag",x.imag)
        if x.imag == 0:
          x = x.real
          # print("Barycenters[iElemShock]",Barycenters[iElemShock])
          if x > Barycenters[iElemShock] and x<0:
            x0real.append(x)
      # print("x0real:",x0real)
    else:
      listregtype.append('exp')
      relafac = max(DataVar[iElemShock:])
      # print('relafac: %f' % relafac)
      # print(DataVar[iElemShock:])
      a,pcov = scipy.optimize.curve_fit(exp_reg_func,Barycenters[iElemShock:],DataVar[iElemShock:]/relafac,method='trf',max_nfev=1E15,bounds=(0, np.inf))
      a[0]*=relafac;a[2]*=relafac;a[4]*=relafac
      listcoeffs.append(a)

    # print("Regression type: %s" % listregtype[-1])
    # print("a-array: ",listcoeffs[-1])

  x0=np.median(np.array(x0real))
  EquiCon = np.zeros(len(VarNames),dtype=float)
  for iData in range(upperbound):
    VarName = VarNames[iData]
    DataVar = Data[iData]
    a = listcoeffs[iData]
    if(listregtype[iData]=='const'):
      yEqui = a
    if(listregtype[iData]=='poly'):
      yEqui = 0
      # print("a",a)
      # print("x0",x0)
      for n in range(PolynomOrder+1):
        yEqui+=a[n]*x0**n
      # print("yEqui",yEqui)
      Output[0]+=',%s,%s_interpol,Equi'%(VarName,VarName)
      for iElem in range(Mesh.nElems()):
        if(iElem>=iElemShock):
          y=0
          for n in range(PolynomOrder+1):
            y+=a[n]*Barycenters[iElem]**n
          Output[iElem+1]+=',%f,%f,%f'%(DataVar[iElem],y,yEqui)
        else:
          Output[iElem+1]+=',%f,,%f'%(DataVar[iElem],yEqui)
    elif (listregtype[iData]=='exp'):
      yEqui = a[4]
      Output[0]+=',%s,%s_interpol,Equi'%(VarName,VarName)
      for iElem in range(Mesh.nElems()):
        if(iElem>=iElemShock):
          y=exp_reg_func(Barycenters[iElem],*a)
          Output[iElem+1]+=',%f,%f,%f'%(DataVar[iElem],y,yEqui)
        else:
          Output[iElem+1]+=',%f,,%f'%(DataVar[iElem],yEqui)
    EquiCon[iData] = yEqui

  if(not Chem):
    # If no chemical reactions: relative composition should not change
    TotalDens=0
    for iSpec in range(1,Simulation.nSpecies+1):
      TotalDens = TotalDens + float(SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()])
    for iSpec in range(1,Simulation.nSpecies+1):
      EquiCon[iSpec+upperbound-1] = EquiCon[0]/TotalDens*float(SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()])

  if (Simulation.ElectronSpecies>0):
    EquiCon[Simulation.ElectronSpecies+2] = 0.0
    for iSpec in range(1,Simulation.nSpecies+1):
      if(iSpec==Simulation.ElectronSpecies):
        continue
      if(('Part-Species'+str(iSpec)+'-ChargeIC'.lower()) in SimulationParameter):
        if(float(SimulationParameter['Part-Species'+str(iSpec)+'-ChargeIC'.lower()])>0):
          EquiCon[Simulation.ElectronSpecies+2] = EquiCon[Simulation.ElectronSpecies+2] + EquiCon[iSpec+2]

  with open('EquiOutputVandermonde.txt','w') as f:
    for line in Output:
      f.write(line+'\n')

  print('Equilibrium Condition:')
  for i in range(len(VarNames)):
    # print('%s: %e',%(VarNames[i],EquiCon[i]))
    print(VarNames[i]+':'+' '*(len(VarNames[i])-50)+'%e'%EquiCon[i])

  return EquiCon

def exp_reg_func(x,a,b,c,d,e):
  return a*np.exp(-b*x)+c*np.exp(d*x)+e

def ChangeSimulationType(EquiCon):
  """
  Change Parameter for Standing Shock SImulation
  """
  global Iteration, Simulation

  Iteration.StandingShock = True

  Iteration.VOut=-Solution.Shockspeed
  Iteration.VIn=Iteration.CurrentVelo-Solution.Shockspeed
  Iteration.iter = -1
  Simulation.tend = max(Simulation.tend,Mesh.Length/(Iteration.VOut))


  SimulationParameter['Part-Boundary3-Condition'.lower()] = 'open'


  # Set constant particle properties
  for iSpec in range(1,Simulation.nSpecies+1):
    if iSpec in Simulation.InflowSpecies:
      #Init
      SimulationParameter[('Part-Species'+str(iSpec)+'-nInits').lower()] = '2'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-SpaceIC').lower()] = 'cuboid'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-BaseVector1IC').lower()] = '(/0.,1.,0./)'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-BaseVector2IC').lower()] = '(/0.,0.,1./)'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-VeloVecIC').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloVecIC').lower()]
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-velocityDistribution').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-velocityDistribution').lower()]
      if Simulation.nSpecies==1:
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-PartDensity').lower()] = str(EquiCon[0])
      else:
        SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-PartDensity').lower()] = str(EquiCon[iSpec+2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-MWTemperatureIC').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-TempRot').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-TempVib').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init2-TempElec').lower()] = str(EquiCon[2])
      #Surface Flux
      SimulationParameter[('Part-Species'+str(iSpec)+'-nSurfaceFluxBCs').lower()] = '2'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-BC').lower()] = '3'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-velocityDistribution').lower()] = 'maxwell_lpn'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-VeloVecIC').lower()] = \
            SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloVecIC').lower()]
      if Simulation.nSpecies==1:
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-PartDensity').lower()] = str(EquiCon[0])
      else:
        SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-PartDensity').lower()] = str(EquiCon[iSpec+2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-MWTemperatureIC').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempRot').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempVib').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux2-TempElec').lower()] = str(EquiCon[2])




    else:
      #Init
      SimulationParameter[('Part-Species'+str(iSpec)+'-nInits').lower()] = '1'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-SpaceIC').lower()] = 'cuboid'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BaseVector1IC').lower()] = '(/0.,1.,0./)'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-BaseVector2IC').lower()] = '(/0.,0.,1./)'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-VeloVecIC').lower()] = \
            SimulationParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init1-VeloVecIC').lower()]
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-velocityDistribution').lower()] = \
            SimulationParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Init1-velocityDistribution').lower()]
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-PartDensity').lower()] = str(EquiCon[iSpec+2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-MWTemperatureIC').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempRot').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempVib').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Init1-TempElec').lower()] = str(EquiCon[2])
      #Surface Flux
      SimulationParameter[('Part-Species'+str(iSpec)+'-nSurfaceFluxBCs').lower()] = '1'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-BC').lower()] = '3'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-velocityDistribution').lower()] = 'maxwell_lpn'
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-VeloVecIC').lower()] = \
            SimulationParameter[('Part-Species'+str(Simulation.InflowSpecies[0])+'-Surfaceflux1-VeloVecIC').lower()]
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-PartDensity').lower()] = str(EquiCon[iSpec+2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-MWTemperatureIC').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempRot').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempVib').lower()] = str(EquiCon[2])
      SimulationParameter[('Part-Species'+str(iSpec)+'-Surfaceflux1-TempElec').lower()] = str(EquiCon[2])

def IntegrateSolution():
  """
  Integrate DSMC Files into one inegrated file
  """
  global Simulation
  DSMCHOFiles = sorted(glob.glob("*SMCState*"))
  DSMCHOFiles = DSMCHOFiles[len(DSMCHOFiles)-Simulation.UsedValues:]
  for i in range(len(DSMCHOFiles)):
    DSMCHOFile = h5py.File(DSMCHOFiles[i], 'r')
    if i == 0:
      ElemData = np.array(DSMCHOFile['ElemData'])
    else:
      ElemData += np.array(DSMCHOFile['ElemData'])
    DSMCHOFile.close()
  ElemData /= len(DSMCHOFiles)
  with h5py.File("Integrated.h5", 'w', libver=('earliest', 'v110')) as IntegratedFile:
    dset = IntegratedFile.create_dataset("ElemData", np.shape(ElemData), dtype='f8', data=ElemData)
    # dset = ElemData
    with h5py.File(DSMCHOFiles[-1], 'r') as DSMCHOFile:
      for attribute in list(DSMCHOFile.attrs):
        # if sys.version_info[0] == 3:
        #   attribute = attribute.decode()
        # print(attribute)
        # print(type(attribute))
        Enties = []
        for entry in DSMCHOFile.attrs[attribute]:
          # print("type(entry)",type(entry))
          if type(entry) is np.bytes_ :
            entry = entry.decode()
          else:
            entry = str(entry)
          # print("type(entry)",type(entry))
          # print("entry",entry)
          # Enties.append(entry+" "*(256-len(entry)))
          # Enties.append(entry.encode("ascii"))
          Enties.append(entry)
        #   print("Enties[-1]",Enties[-1])
        #   print("type(Enties[-1])",type(Enties[-1]))
        #   print("entry.encode('ascii')",entry.encode('ascii'))
        # print("Enties",Enties)
        # print("np.array(Enties)",np.array(Enties))
        # print("np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype)",np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype))
        # print("np.array(Enties,dtype=h5py.string_dtype('utf-8', 255))",np.array(Enties,dtype=h5py.string_dtype('utf-8', 255)))
        # IntegratedFile.attrs[attribute] = np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype)
        IntegratedFile.attrs[attribute] = np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype)
        # IntegratedFile.attrs.create(attribute,entry)
        # IntegratedFile.attrs[attribute] = np.array(Enties.encode("utf-8"),dtype=utf8_type)
        # print("DSMCHOFile.attrs[attribute].dtype",DSMCHOFile.attrs[attribute].dtype)
        # print("IntegratedFile.attrs[attribute]",IntegratedFile.attrs[attribute])
  if not h5py.is_hdf5('Integrated.h5'):
      raise ValueError('Not an hdf5 file')



def IntegrateSolution2():
  """
  Integrate DSMC Files into one inegrated file
  """
  global Simulation
  DSMCHOFiles = sorted(glob.glob("*SMCState*"))
  DSMCHOFiles = DSMCHOFiles[len(DSMCHOFiles)-Simulation.UsedValues:]
  shutil.copyfile(DSMCHOFiles[-1], "Integrated.h5")
  for i in range(len(DSMCHOFiles)):
    DSMCHOFile = h5py.File(DSMCHOFiles[i], 'r')
    if i == 0:
      ElemData = np.array(DSMCHOFile['ElemData'])
    else:
      ElemData += np.array(DSMCHOFile['ElemData'])
    DSMCHOFile.close()
  ElemData /= len(DSMCHOFiles)
  with h5py.File("Integrated.h5", 'r+') as IntegratedFile:
    #, libver=('earliest', 'v110')
    # dset = IntegratedFile.create_dataset("ElemData", np.shape(ElemData), dtype='f8', data=ElemData)
    print(IntegratedFile['ElemData'])
    IntegratedFile['ElemData'][:,:] = ElemData[:,:]
    # dset = ElemData
    # with h5py.File(DSMCHOFiles[-1], 'r') as DSMCHOFile:
    #   for attribute in list(DSMCHOFile.attrs):
    #     # if sys.version_info[0] == 3:
    #     #   attribute = attribute.decode()
    #     # print(attribute)
    #     # print(type(attribute))
    #     Enties = []
    #     for entry in DSMCHOFile.attrs[attribute]:
    #       # print("type(entry)",type(entry))
    #       if type(entry) is np.bytes_ :
    #         entry = entry.decode()
    #       else:
    #         entry = str(entry)
    #       # print("type(entry)",type(entry))
    #       # print("entry",entry)
    #       # Enties.append(entry+" "*(256-len(entry)))
    #       # Enties.append(entry.encode("ascii"))
    #       Enties.append(entry)
    #     #   print("Enties[-1]",Enties[-1])
    #     #   print("type(Enties[-1])",type(Enties[-1]))
    #     #   print("entry.encode('ascii')",entry.encode('ascii'))
    #     # print("Enties",Enties)
    #     # print("np.array(Enties)",np.array(Enties))
    #     # print("np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype)",np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype))
    #     # print("np.array(Enties,dtype=h5py.string_dtype('utf-8', 255))",np.array(Enties,dtype=h5py.string_dtype('utf-8', 255)))
    #     # IntegratedFile.attrs[attribute] = np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype)
    #     IntegratedFile.attrs[attribute] = np.array(Enties,dtype=DSMCHOFile.attrs[attribute].dtype)
    #     # IntegratedFile.attrs.create(attribute,entry)
    #     # IntegratedFile.attrs[attribute] = np.array(Enties.encode("utf-8"),dtype=utf8_type)
    #     # print("DSMCHOFile.attrs[attribute].dtype",DSMCHOFile.attrs[attribute].dtype)
    #     # print("IntegratedFile.attrs[attribute]",IntegratedFile.attrs[attribute])
  if not h5py.is_hdf5('Integrated.h5'):
      raise ValueError('Not an hdf5 file')





# Main Program

CreateTitle()
GenerateReadinValues()

args = ArgumentsParser()

# if(args.restart!=None):
#   pass

ReadShocktubeIni(args.InputFile)
Init(args.Equicon,args.restart)

if(args.restart!=None):
  InitRestart(args.restart)

if(args.Equicon):
  Path.CurrentSimDir = '%s/iter_%02i' % (Path.MasterDir,0)
  os.chdir(Path.CurrentSimDir)
  ReadMeshFile(SimulationParameter["ProjectName".lower()] + '_mesh.h5')

  ReadDSMCHOData()
  DetermineShockSpeed()

  EquiCon = GetEquilibriumCondition2()
  exit()

if not Autoadjust.Mesh:
  GenerateMeshFile()
  ReadMeshFile(SimulationParameter["ProjectName".lower()] + '_mesh.h5')

# Begin Loop
if(not RestartStanding):
  Convergence = False
  while not Convergence:
    if(not DoRestart):
      Iteration.iter += 1
      print('='*128)
      print('Iteration: %i' % Iteration.iter)

      PrepareSimulation()
      if Autoadjust.Mesh:
        GenerateMeshFile()
        ReadMeshFile(SimulationParameter["ProjectName".lower()] + '_mesh.h5')
      InitPiclasParameters()
      WriteParameterIni()

      RunPiclas()
    else:
      DoRestart = False
      print('='*128)
      print('Analyzing last iteration: %i' % Iteration.iter)
      os.chdir(Path.CurrentSimDir)
      if Autoadjust.Mesh:
        ReadMeshFile(SimulationParameter["ProjectName".lower()] + '_mesh.h5')


    ReadDSMCHOData()
    DetermineShockSpeed()
    Convergence = UpdateVariables()

  print('Matching Shock Condition found')
  print('Calculating equilibrium condition')
  EquiCon = GetEquilibriumCondition2()

  ChangeSimulationType(EquiCon)

Convergence = False
while not Convergence:
  if(not DoRestart):
    Iteration.iter += 1
    print('='*128)
    print('Standing Shock')
    print('Iteration: %i' % Iteration.iter)

    PrepareSimulation()
    if Autoadjust.Mesh:
      GenerateMeshFile()
      ReadMeshFile(SimulationParameter["ProjectName".lower()] + '_mesh.h5')
    InitPiclasParameters()
    WriteParameterIni()

    RunPiclas()

    ReadDSMCHOData()
    DetermineShockSpeed()
    Convergence = UpdateVariables()
  else:
    DoRestart = False
    print('='*128)
    print('Analyzing last standing iteration: %i' % Iteration.iter)
    os.chdir(Path.CurrentSimDir)

    if Autoadjust.Mesh:
      ReadMeshFile(SimulationParameter["ProjectName".lower()] + '_mesh.h5')

    ReadDSMCHOData()
    DetermineShockSpeed()
    Convergence = UpdateVariables()
    Convergence = False





print('FINISH')
