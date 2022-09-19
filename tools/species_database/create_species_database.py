import numpy as np
import h5py
from argparse import ArgumentParser

database_output = "Species_Database.h5"
database_electronic = "Electronic-State-Database.h5"
database_crosssection = "XSec_Database_H2_Plasma.h5"

hdf = h5py.File(database_output, 'a')
h5_electronic = h5py.File(database_electronic, 'r')
h5_crosssection = h5py.File(database_crosssection, 'r')

parser = ArgumentParser(prog='create_species_database')
parser.add_argument("-f", "--file", dest="ini_filename",
                    help="DSMC.ini file to read-in parameters", metavar="FILE")

args = parser.parse_args()

def is_float(value):
  try:
    float(value)
    return True
  except:
    return False

# Create general structure
hdf_species_group = 'Species'
if hdf_species_group in hdf.keys():
  print('Group Species already exists.')
  hdf_species_group = hdf['Species']
else:
  print('Group Species does not exist, creating new group.')
  hdf_species_group = hdf.create_group('Species')

hdf_xsec_group = 'Cross-Sections'
if hdf_xsec_group in hdf.keys():
  print('Group Cross-Sections already exists.')
  hdf_xsec_group = hdf['Cross-Sections']
else:
  print('Group Cross-Sections does not exist, creating new group.')
  hdf_xsec_group = hdf.create_group('Cross-Sections')

# Read-in of DSMC.ini parameters and electronic state
with open(args.ini_filename) as file:
  for line in file:
    if not line.startswith('!'):
      if line.startswith('Part-'):
        var_name = line.strip().replace(" ","").replace("Part-Species","").split('=')[0].split('-')
        var_value = line.strip().replace(" ","").replace("Part-Species","").split('=')[1].split('!', 1)[0]
        if var_name[1].startswith('SpeciesName'):
          hdf_species = var_value
          if hdf_species in hdf_species_group.keys():
            print('Species already exists: ', hdf_species)
            species_count = var_name[0]
            hdf_species = hdf_species_group[hdf_species]
          elif hdf_species == 'electron':
            print('Species added to the database: ', hdf_species)
            hdf_species = hdf_species_group.create_dataset(hdf_species,data=[0])
          else:
            hdf_input_data = h5_electronic[hdf_species]
            print('Species added to the database: ', hdf_species)
            hdf_species = hdf_species_group.create_dataset(hdf_species,data=hdf_input_data)
        else:
          if var_name[1] in hdf_species.attrs:
            print('Species parameter is already set: ', var_name[1])
          else:
            if is_float(var_value):
              var_value = float(var_value)
            hdf_species.attrs[var_name[1]] = var_value
            print('Species parameter set: ', var_name[1])

# Copy cross-section data
for dataset in h5_crosssection.keys():
  if dataset in hdf_xsec_group.keys():
    print('Cross-section is already set: ', dataset)
  else:
    print('Cross-section added: ', dataset)
    hdf_xsec_group.copy(source=h5_crosssection[dataset],dest=hdf_xsec_group)

hdf.close()
h5_electronic.close()
h5_crosssection.close()