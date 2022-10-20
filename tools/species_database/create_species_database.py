from cmath import log
from matplotlib.rcsetup import validate_bool
import numpy as np
import h5py
from argparse import ArgumentParser
from datetime import date

parser = ArgumentParser(prog='create_species_database')
parser.add_argument("-p", "--parameter", dest="ini_filename",
                    help="DSMC.ini file to read-in parameters", metavar="FILE", default="DSMC.ini")
parser.add_argument("-e", "--electronic", dest="database_electronic",
                    help="Electronic excitation database to read-in parameters", metavar="FILE", default="Electronic-State-Database.h5")
parser.add_argument("-c", "--crosssection", dest="database_crosssection",
                    help="Crosssection database to read-in parameters", metavar="FILE", default="XSec_Database.h5")
parser.add_argument("-o", "--output", dest="database_output",
                    help="Output file name", metavar="FILE", default="Species_Database.h5")
parser.add_argument("-r", "--reference", dest="reference",
                    help="Reference for species data", default="Unknown")
parser.add_argument("-s", "--radiation", dest="rad_file_name",
                    help="Radiation file to read-in parameters", metavar="FILE", default="Rad.dat")

args = parser.parse_args()

h5_species      = h5py.File(args.database_output, 'a')
h5_electronic   = h5py.File(args.database_electronic, 'r')
h5_crosssection = h5py.File(args.database_crosssection, 'r')

def is_float(value):
  try:
    float(value)
    return True
  except:
    return False
    
# Create general structure
hdf_species_group = 'Species'
if hdf_species_group in h5_species.keys():
  print('Group Species already exists.')
  hdf_species_group = h5_species['Species']
  hdf_species_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Species does not exist, creating new group.')
  hdf_species_group = h5_species.create_group('Species')
  hdf_species_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")

hdf_xsec_group = 'Cross-Sections'
if hdf_xsec_group in h5_species.keys():
  print('Group Cross-Sections already exists.')
  hdf_xsec_group = h5_species['Cross-Sections']
  hdf_xsec_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Cross-Sections does not exist, creating new group.')
  hdf_xsec_group = h5_species.create_group('Cross-Sections')
  hdf_xsec_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")

hdf_rad_group = 'Radiation'
if hdf_rad_group in h5_species.keys():
  print('Group Radiation already exists.')
  hdf_rad_group = h5_species['Radiation']
  hdf_rad_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Radiation does not exist, creating new group.')
  hdf_rad_group = h5_species.create_group('Radiation')
  hdf_rad_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")
  
hdf_reac_group = 'Reaction'  
if hdf_reac_group in h5_species.keys():
  print('Group Reaction already exists.')
  hdf_reac_group = h5_species['Reaction']
  hdf_reac_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Reaction does not exist, creating new group.')
  hdf_reac_group = h5_species.create_group('Reaction')
  hdf_reac_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")
  
hdf_surf_group = 'Surface-Chemistry'
if hdf_surf_group in h5_species.keys():
  print('Group Surface-Chemistry already exists.')
  hdf_surf_group = h5_species['Surface-Chemistry']
  hdf_surf_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Surface-Chemistry does not exist, creating new group.')
  hdf_surf_group = h5_species.create_group('Surface-Chemistry')
  hdf_surf_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")

logical_list = ['PolyatomicMol', 'LinearMolec']

# Read-in of DSMC.ini parameters and electronic state
with open(args.ini_filename) as file:
  spec_dict = {}
  for line in file:
    if not line.startswith('!'):
      if line.startswith('Part-'):
        var_name = line.strip().replace(" ","").replace("Part-Species","").split('=')[0].split('-')
        var_value = line.strip().replace(" ","").replace("Part-Species","").split('=')[1].split('!', 1)[0]
        if var_name[1].startswith('SpeciesName'):
          hdf_species = var_value
          species_count = var_name[0]
          spec_dict[species_count] = hdf_species
          if hdf_species in hdf_species_group.keys():
            print('Species already exists: ', hdf_species)
            hdf_species = hdf_species_group[hdf_species]
          elif hdf_species == 'electron':
            print('Species added to the database: ', hdf_species)
            hdf_species = hdf_species_group.create_dataset(hdf_species,data=[0])
            hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")
          else:
            if hdf_species in h5_electronic.keys():
              hdf_input_data = h5_electronic[hdf_species]
              print('Species added to the database: ', hdf_species)
              hdf_species = hdf_species_group.create_dataset(hdf_species,data=hdf_input_data)
              hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")
            else:
              print('Species added to the database, but electronic levels are unknown: ', hdf_species)
              hdf_species = hdf_species_group.create_dataset(hdf_species,data=[0])
              hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")

with open(args.ini_filename) as file:
  for line in file:
    if not line.startswith('!'):
      if line.startswith('Part-'):
        var_name = line.strip().replace(" ","").replace("Part-Species","").split('=')[0].split('-')
        var_value = line.strip().replace(" ","").replace("Part-Species","").split('=')[1].split('!', 1)[0]

        if not var_name[1].startswith('SpeciesName'):
          species_count = var_name[0]
          hdf_species = spec_dict[species_count]
          hdf_species = hdf_species_group[hdf_species]
          if var_name[1] in hdf_species.attrs:
            print('Species parameter is already set: ', var_name[1]) #raus??
          else:
            if is_float(var_value):
              var_value = float(var_value)
            if var_name[1] not in logical_list:
              hdf_species.attrs[var_name[1]] = var_value
              print('Species parameter set: ', var_name[1]) #raus??
            else:
              if 'F' in var_value or 'false' in var_value:
                var_value = 0
              else:
                var_value = 1
              hdf_species.attrs[var_name[1]] = var_value
              print('Species parameter set: ', var_name[1]) #raus??
            # Write attributes for source and time of retrieval
            hdf_species.attrs['* Reference'] = args.reference
            hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")

# Copy cross-section data
for dataset in h5_crosssection.keys():
  if dataset in hdf_xsec_group.keys():
    print('Cross-section is already set: ', dataset)
  else:
    print('Cross-section added: ', dataset)
    hdf_xsec_group.copy(source=h5_crosssection[dataset],dest=hdf_xsec_group)
    
# Radiation data
species_list = ['N','O','NIon1','OIon1','H','Xe','XeIon1','XeIon2','Ar','ArIon1','ArIon2']
for rad_spec in species_list:
  if rad_spec in hdf_rad_group.keys():
    print('Radiative species already exists: ', rad_spec)
    rad_spec_group = hdf_rad_group[rad_spec]
  else:
    print('Radiative species added to the database: ', rad_spec)
    rad_spec_group = hdf_rad_group.create_group(rad_spec)
    rad_spec_group.attrs['* Created']   = date.today().strftime("%B %d, %Y")
    lines = rad_spec_group.create_dataset('Lines',data=[0])
    levels = rad_spec_group.create_dataset('Levels',data=[0])
  
reac_attr_list = ['Reactants', 'Products', 'NonReactiveSpecies']
  
# Read-in of gas-phase reaction data
with open(args.ini_filename) as file:
  reac_dict = {}
  for line in file:
    if not line.startswith('!'):
      if line.startswith('DSMC-Reaction'):
        var_name = line.strip().replace(" ","").replace("DSMC-Reaction","").split('=')[0].split('-',1)
        var_value = line.strip().replace(" ","").replace("DSMC-Reaction","").split('=')[1].split('!', 1)[0]
        if var_name[1].startswith('ReactionName'):
          hdf_reac = var_value
          reac_count = var_name[0]
          reac_dict[reac_count] = hdf_reac
          if hdf_reac in hdf_reac_group.keys():
            print('Reaction already exists: ', hdf_reac)
            hdf_reac = hdf_reac_group[hdf_reac]
          else:
            print('Reaction added to the database: ', hdf_reac)
            hdf_reac = hdf_reac_group.create_dataset(hdf_reac,data=[0])
            hdf_reac.attrs['* Created']   = date.today().strftime("%B %d, %Y")

with open(args.ini_filename) as file:
  for line in file:
    if not line.startswith('!'):
      if line.startswith('DSMC-Reaction'):
        var_name = line.strip().replace(" ","").replace("DSMC-Reaction","").split('=')[0].split('-',1)
        var_value = line.strip().replace(" ","").replace("DSMC-Reaction","").split('=')[1].split('!', 1)[0]

        if not var_name[1].startswith('ReactionName'):
          reac_count = var_name[0]
          hdf_reac = reac_dict[reac_count]
          hdf_reac = hdf_reac_group[hdf_reac]
          if var_name[1] in hdf_reac.attrs:
            print('Reaction parameter is already set: ', var_name[1])
          else:
            if is_float(var_value):
              var_value = float(var_value)
            if var_name[1] not in reac_attr_list:
              hdf_reac.attrs[var_name[1]] = var_value
            else: 
              spec_name_list = ''
              var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
              for val in var_value:
                spec_name_list = spec_name_list + ', ' + spec_dict[val]
              spec_name_list = spec_name_list.replace(', ', '', 1)
              hdf_reac.attrs[var_name[1]] = spec_name_list
            print('Reaction parameter set: ', var_name[1])
            # Write attributes for source and time of retrieval
            hdf_reac.attrs['* Reference'] = args.reference
            hdf_reac.attrs['* Created']   = date.today().strftime("%B %d, %Y")

surf_attr_list = ['Reactants', 'Products']
surf_wo_readin = ['Boundaries', 'NumOfBoundaries', 'Inhibition', 'Promotion']

# Read-in of gas-surface reaction data
with open(args.ini_filename) as file:
  surf_dict = {}
  for line in file:
    if not line.startswith('!'):
      if line.startswith('Surface-Reaction'):
        var_name = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[0].split('-')
        var_value = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[1].split('!', 1)[0]
        if var_name[1].startswith('SurfName'):
          hdf_surf = var_value
          surf_count = var_name[0]
          surf_dict[surf_count] = hdf_surf
          if hdf_surf in hdf_surf_group.keys():
            print('Surface reaction already exists: ', hdf_surf)
            hdf_surf = hdf_surf_group[hdf_surf]
          else:
            print('Surface reaction added to the database: ', hdf_surf)
            hdf_surf = hdf_surf_group.create_dataset(hdf_surf,data=[0])
            hdf_surf.attrs['* Created']   = date.today().strftime("%B %d, %Y")

with open(args.ini_filename) as file:
  for line in file:
    if not line.startswith('!'):
      if line.startswith('Surface-Reaction'):
        var_name = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[0].split('-')
        var_value = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[1].split('!', 1)[0]

        if not var_name[1].startswith('SurfName'):
          surf_count = var_name[0]
          hdf_surf = surf_dict[surf_count]
          hdf_surf = hdf_surf_group[hdf_surf]
          if var_name[1] in hdf_surf.attrs:
            print('Catalytic parameter is already set: ', var_name[1])
          else:
            if is_float(var_value):
              var_value = float(var_value)
            if var_name[1] not in surf_attr_list and var_name[1] not in surf_wo_readin:
              hdf_surf.attrs[var_name[1]] = var_value
              print('Catalytic parameter set: ', var_name[1]) 
            elif var_name[1] in surf_attr_list:
              spec_name_list = ''
              var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
              for val in var_value:
                spec_name_list = spec_name_list + ', ' + spec_dict[val]
              spec_name_list = spec_name_list.replace(', ', '', 1)
              hdf_surf.attrs[var_name[1]] = spec_name_list
              print('Catalytic parameter set: ', var_name[1]) 
            # Write attributes for source and time of retrieval
            hdf_surf.attrs['* Reference'] = args.reference
            hdf_surf.attrs['* Created']   = date.today().strftime("%B %d, %Y")
    

print('*** Database successfully built. ***')

h5_species.close()
h5_electronic.close()
h5_crosssection.close()
