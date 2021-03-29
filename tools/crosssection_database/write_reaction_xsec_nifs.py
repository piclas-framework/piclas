import numpy as np
import h5py

### User-input
# Output database
database_output = "XSec_Database_Xe_Plasma.h5"
# Reactants and products
pair_list =     ['Xe-electron',                 'Xe-electron',                      'XeIon1-electron']
reaction_list = ['XeIon1-electron-electron',    'XeIon2-electron-electron-electron','XeIon2-electron-electron']
# Read-in input file, skip the first 8 lines of the header (NIFS output copied starting with Record No)
input_list =    ['NIFS-Xe-elec-XeIon-2elec.dat','NIFS-Xe-elec-XeIon2-3elec.dat',    'NIFS-XeIon-elec-XeIon2-2e.dat']
source_list =   ['NIFS Database: Krishnakumar, E et al. J. Phys. B 21 (1988) 1055','NIFS Database: Krishnakumar, E et al. J. Phys. B 21 (1988) 1055','NIFS Database: Man, K.F. et al. J. Phys. B 20 (1987) 5865']
########################################################################################################

hdf = h5py.File(database_output, 'r+')
for iReac, current_reaction in enumerate(reaction_list):
  print('Reactants: ' + pair_list[iReac])
  print('Products: ' + current_reaction)
  data_input = np.genfromtxt(input_list[iReac], delimiter=',', skip_header=8, usecols=(0,1))
  # Sort data input by ascending energy values
  data_input.view('f8,f8').sort(order=['f0'], axis=0)

  # Convert cm2 to m2
  data_input[:,1] *= 1E-4
  print('Input converted from cm2 to m2.')

  # Check if the collision pair (=folder) already exists, if not create it
  hdf_pair = pair_list[iReac]
  if hdf_pair in hdf.keys():
    print('Collision pair already exists.')
    hdf_pair = hdf[pair_list[iReac]]
  else:
    print('Collision pair does not exist, creating new pair.')
    hdf_pair = hdf.create_group(pair_list[iReac])

  # Check if a REACTION group (=folder) already exists, if not create it
  hdf_reaction = 'REACTION'
  if hdf_reaction in hdf_pair.keys():
    print('Group REACTION already exists.')
    hdf_reaction = hdf_pair['REACTION']
  else:
    print('Group REACTION does not exist, creating new group.')
    hdf_reaction = hdf_pair.create_group('REACTION')

  # If reaction data set already exists, delete the old set
  if current_reaction in hdf_reaction.keys():
    del hdf_reaction[current_reaction]
    print('Old dataset replaced.')
  hdf_reaction.create_dataset(current_reaction, data=data_input)

  hdf_reaction[current_reaction].attrs['Source'] = source_list[iReac]

hdf.close()