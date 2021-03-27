import numpy as np
import h5py

# Output database
database_output = "XSec_Database_Xe_Plasma.h5"
# Open existing file/create new
hdf = h5py.File(database_output, 'r+')
hdf.attrs['Date'] = '16.03.2021'

## User-input
# pair = 'Xe-electron'
# reaction = 'XeIon1-electron-electron'
# # Read-in input file, skip the first 8 lines of the header (NIFS output copied starting with Record No)
# data_input = np.genfromtxt('NIFS-Xe-elec-XeIon-2elec.dat', dtype=float, delimiter=',', skip_header=8, usecols=(0,1))
# source_input = 'NIFS Database: Krishnakumar, E et al. J. Phys. B 21 (1988) 1055'
########################################################################################################
## User-input
# pair = 'Xe-electron'
# reaction = 'XeIon2-electron-electron-electron'
# # Read-in input file, skip the first 8 lines of the header (NIFS output copied starting with Record No)
# data_input = np.genfromtxt('NIFS-Xe-elec-XeIon2-3elec.dat', delimiter=',', skip_header=8, usecols=(0,1))
# source_input = 'NIFS Database: Krishnakumar, E et al. J. Phys. B 21 (1988) 1055'
########################################################################################################
## User-input
pair = 'XeIon1-electron'
reaction = 'XeIon2-electron-electron'
# Read-in input file, skip the first 8 lines of the header (NIFS output copied starting with Record No)
data_input = np.genfromtxt('NIFS-XeIon-elec-XeIon2-2e.dat', delimiter=',', skip_header=8, usecols=(0,1))
source_input = 'NIFS Database: Man, K.F. et al. J. Phys. B 20 (1987) 5865'
########################################################################################################

# Sort data input by ascending energy values
print(data_input.view('f8,f8').sort(order=['f0'], axis=0))

# Convert cm2 to m2
data_input[:,1] *= 1E-4
print('Input converted from cm2 to m2.')

# Check if the collision pair (=folder) already exists, if not create it
hdf_pair = pair
if hdf_pair in hdf.keys():
  print('Collision pair already exists.')
  hdf_pair = hdf[pair]
else:
  print('Collision pair does not exist, creating new pair.')
  hdf_pair = hdf.create_group(pair)

# Check if a REACTION group (=folder) already exists, if not create it
hdf_reaction = 'REACTION'
if hdf_reaction in hdf_pair.keys():
  print('Group REACTION already exists.')
  hdf_reaction = hdf_pair['REACTION']
else:
  print('Group REACTION does not exist, creating new group.')
  hdf_reaction = hdf_pair.create_group('REACTION')

# If reaction data set already exists, delete the old set
if reaction in hdf_reaction.keys():
  del hdf_reaction[reaction]
  print('Old dataset replaced.')
hdf_reaction.create_dataset(reaction, data=data_input)

hdf_reaction[reaction].attrs['Source'] = source_input

hdf.close()