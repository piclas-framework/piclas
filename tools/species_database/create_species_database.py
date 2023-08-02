from cmath import log
from matplotlib.rcsetup import validate_bool
import numpy as np
import h5py
from argparse import ArgumentParser
from datetime import date
import re
from functions_database import custom_sort_reactants, custom_sort_products

parser = ArgumentParser(prog='create_species_database')
parser.add_argument("-p", "--parameter", dest="ini_filename",
                    help="DSMC.ini file to read-in parameters", metavar="FILE", default="DSMC.ini")
parser.add_argument("-e", "--electronic", dest="database_electronic",
                    help="Electronic excitation database to read-in parameters", metavar="FILE", default="Electronic-State-Database.h5")
parser.add_argument("-c", "--crosssection", dest="database_crosssection",
                    help="Crosssection database to read-in parameters", metavar="FILE", default="XSec_Database.h5")
parser.add_argument("-o", "--output", dest="database_output",
                    help="Output file name", metavar="FILE", default="SpeciesDatabase.h5")
parser.add_argument("-r", "--reference", dest="reference",
                    help="Reference for species data", default="Unknown")
# parser.add_argument("-s", "--radiation", dest="rad_file_name",
#                     help="Radiation file to read-in parameters", metavar="FILE", default="Rad.dat")

args = parser.parse_args()

h5_species      = h5py.File(args.database_output, 'a')
h5_electronic   = h5py.File(args.database_electronic, 'r')
h5_crosssection = h5py.File(args.database_crosssection, 'r')

# Name of the DSMC.ini is the name of the ChemicalModel
model_name = args.ini_filename[:-4]
# Split the name in its components (General structure: Topic-SpecNum-ReacNum-Source)
items = re.split('_+', model_name)
indices_to_modify = [1,2]
model_name_split = [item[:-4] if i in indices_to_modify else item for i, item in enumerate(items)]

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

# hdf_rad_group = 'Radiation'
# if hdf_rad_group in h5_species.keys():
#   print('Group Radiation already exists.')
#   hdf_rad_group = h5_species['Radiation']
#   hdf_rad_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
# else:
#   print('Group Radiation does not exist, creating new group.')
#   hdf_rad_group = h5_species.create_group('Radiation')
#   hdf_rad_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")
  
hdf_reac_group = 'Reactions'  
if hdf_reac_group in h5_species.keys():
  print('Group Reaction already exists.')
  hdf_reac_group = h5_species['Reactions']
  hdf_reac_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Reactions does not exist, creating new group.')
  hdf_reac_group = h5_species.create_group('Reactions')
  hdf_reac_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")
  
# hdf_surf_group = 'Surface-Chemistry'
# if hdf_surf_group in h5_species.keys():
#   print('Group Surface-Chemistry already exists.')
#   hdf_surf_group = h5_species['Surface-Chemistry']
#   hdf_surf_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
# else:
#   print('Group Surface-Chemistry does not exist, creating new group.')
#   hdf_surf_group = h5_species.create_group('Surface-Chemistry')
#   hdf_surf_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")

# Variables that must be treated separately
logical_list = ['PolyatomicMol', 'LinearMolec']
spec_attr_list = ['PreviousState']

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
          # Check if the species data already exists in the database
          if hdf_species in hdf_species_group.keys():
            if hdf_species in h5_electronic.keys():
              del hdf_species_group[hdf_species]
              hdf_input_data = h5_electronic[hdf_species]
              print('Electronic states added to the database: ', hdf_species)
              hdf_species = hdf_species_group.create_dataset(hdf_species,data=hdf_input_data)
              hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")
          elif hdf_species == 'electron':
            print('Species added to the database: ', hdf_species)
            hdf_species = hdf_species_group.create_dataset(hdf_species,data=[0])
            hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")
          else:
            # Add the electronic state data
            if hdf_species in h5_electronic.keys():
              hdf_input_data = h5_electronic[hdf_species]
              print('Electronic states added to the database: ', hdf_species)
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
          # Check if the species parameter is already defiend
          if var_name[1] in hdf_species.attrs:
            print('Species parameter is already set: ', var_name[1]) 
          else:
            # Float conversion
            if is_float(var_value):
              var_value = float(var_value)
            if var_name[1] not in logical_list and var_name[1] not in spec_attr_list:
              hdf_species.attrs[var_name[1]] = var_value
              print('Species parameter set: ', var_name[1]) 
              # Previous-state read in as a string
            elif var_name[1] in spec_attr_list:
              var_value = str(int(var_value))
              spec_name_list = spec_dict[var_value]
              hdf_species.attrs[var_name[1]] = np.string_(spec_name_list)
              print('Species parameter set: ', var_name[1]) 
              # Logicals
            else:
              if 'F' in var_value or 'false' in var_value:
                var_value = 0
              else:
                var_value = 1
              hdf_species.attrs[var_name[1]] = var_value
              print('Species parameter set: ', var_name[1]) 
            # Write attributes for source and time of retrieval
            hdf_species.attrs['* Reference'] = args.reference
            hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")

# Copy cross-section data if not defined already
for dataset in h5_crosssection.keys():
  if dataset in hdf_xsec_group.keys():
    print('Cross-section is already set: ', dataset)
  else:
    print('Cross-section added: ', dataset)
    hdf_xsec_group.copy(source=h5_crosssection[dataset],dest=hdf_xsec_group)
    
# # Radiation data
# species_list = ['N','O','NIon1','OIon1','H','Xe','XeIon1','XeIon2','Ar','ArIon1','ArIon2']
# for rad_spec in species_list:
#   if rad_spec in hdf_rad_group.keys():
#     print('Radiative species already exists: ', rad_spec)
#     rad_spec_group = hdf_rad_group[rad_spec]
#   else:
#     print('Radiative species added to the database: ', rad_spec)
#     rad_spec_group = hdf_rad_group.create_group(rad_spec)
#     rad_spec_group.attrs['* Created']   = date.today().strftime("%B %d, %Y")
#     lines = rad_spec_group.create_dataset('Lines',data=[0])
#     levels = rad_spec_group.create_dataset('Levels',data=[0])

# Lists for the renaming of the chemical reaction equations
educt_attr_list = ['Reactants']
product_attr_list = ['Reactants', 'Products']
educt_dict = {}
product_dict = {}
reac_attr_list = ['Reactants', 'Products', 'NonReactiveSpecies']

ReacName_dict = {}
reac_dict = {}

# Read-in of the gas-phase reaction names or create them separately if not defined
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
        if var_name[1] in educt_attr_list:
          var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
          educt_dict[var_name[0]] = var_value
        elif var_name[1] in product_attr_list:
          var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
          product_dict[var_name[0]] = var_value
        if var_name[0] not in ReacName_dict:
          ReacName_dict[var_name[0]] = []
          ReacName_dict[var_name[0]].append(var_name[1])
        else:
          ReacName_dict[var_name[0]].append(var_name[1])

# If not defined, set the reaction name: Educts1+Educt2_Product1_Product2
# Non-reactives (M or A) are added at the end of the reactants and at the second position for the products
for key in ReacName_dict:
  if 'ReactionName' not in ReacName_dict[key]:
    ReactionName = ''
    for val in educt_dict[key]:
      ReactionName = ReactionName + '+' + spec_dict[val]
    ReactionName = ReactionName + '_'
    for val in product_dict[key]:
      ReactionName = ReactionName + '+' + spec_dict[val]
    ReactionName = ReactionName[1:]
    ReactionName = ReactionName.replace('_+','_')
    reac_dict[key] = ReactionName
  
  hdf_reac = reac_dict[key]
  ReactionName = hdf_reac

  # Bring the reaction equation in the correct order
  ReactionNameSplit = re.split('_', ReactionName)
  checkeduct_list = re.split('\+', ReactionNameSplit[0])
  checkeduct_list_sorted = custom_sort_reactants(checkeduct_list)
  checkproduct_list = re.split('\+',ReactionNameSplit[1])
  checkproduct_list_sorted = custom_sort_products(checkproduct_list)

# Sorting of the species in the reaction equation and create the new name for the reaction
  ReactionName = ''
  for val in checkeduct_list_sorted:
    ReactionName = ReactionName + '+' + val
  ReactionName = ReactionName + '_'
  for val in checkproduct_list_sorted:
    ReactionName = ReactionName + '+' + val
  ReactionName = ReactionName[1:]
  ReactionName = ReactionName.replace('_+','_')
  reac_dict[key] = ReactionName

  hdf_reac = reac_dict[key]
  ReactionName = hdf_reac

  # Check if the reaction exists
  if hdf_reac in hdf_reac_group.keys():
    print('Reaction already exists: ', hdf_reac)
    hdf_reac = hdf_reac_group[hdf_reac]
  elif (hdf_reac + '#1') in hdf_reac_group.keys():
    print('Reaction already exists: ', hdf_reac)
    ReactionName = hdf_reac + '#1'
    reac_dict[key] = ReactionName
    hdf_reac = hdf_reac_group[(hdf_reac + '#1')]
  else:
    print('Reaction added to the database: ', hdf_reac)
    hdf_reac = hdf_reac_group.create_dataset(hdf_reac,data=[0])
    hdf_reac.attrs['* Created']   = date.today().strftime("%B %d, %Y")

  #Read-In of the Chemistry-model for the reaction
  if 'ChemistryModel' in hdf_reac.attrs:
    model_name_list = []
    model_attr = list(hdf_reac.attrs['ChemistryModel'])
    for val in model_attr:
      model_name_list.append(val.decode('UTF-8'))
      # Check if the reaction is already defined for the model
    if model_name in model_name_list:
      print('This Model is already defined for the reaction. No further action is taken.')
      print('To define the reaction data as a separate instance, please rename the model.')
    else:
      AddNewReaction = False
      Count = 0
      ReacNameTest = ReactionName
      # Check if the model is already defined for the reaction name
      if ReacNameTest in hdf_reac_group.keys():
        hdf_reac_test = hdf_reac_group[ReacNameTest]
        model_test_list = []
        model_attr_test = list(hdf_reac_test.attrs['ChemistryModel'])
        for val in model_attr_test:
          model_test_list.append(val.decode('UTF-8'))
          # If the model is defined, the reaction is not added again
        if model_name in model_test_list:
          print('This Model is already defined for the reaction. No further action is taken.')
          print('To define the reaction data as a separate instance, please rename the model.')
          AddNewReaction = False
          exit
        else: 
          AddNewReaction = True
    # Loop over  iterations of the reaction name and check if the model is already defined
      while ReacNameTest in hdf_reac_group.keys():
        ReacNameTest = re.sub('#\d+', '', ReactionName)
        Count = Count + 1
        ReacNameTest = ReacNameTest + '#' + str(Count)
        if ReacNameTest in hdf_reac_group.keys():
          hdf_reac_test = hdf_reac_group[ReacNameTest]
          model_test_list = []
          model_attr_test = list(hdf_reac_test.attrs['ChemistryModel'])
          for val in model_attr_test:
            model_test_list.append(val.decode('UTF-8'))
          if model_name in model_test_list:
            print('This Model is already defined for the reaction. No further action is taken.')
            print('To define the reaction data as a separate instance, please rename the model.')
            reac_dict[key] = ReacNameTest
            AddNewReaction = False
            exit
          else: 
            AddNewReaction = True
        # Add the reaction with the new model to the database, counter of the reaction name (#) is increased by 1
      if AddNewReaction:
        print('A different Model is already defined for this reaction.')
        Count = 1
        while ReactionName in hdf_reac_group.keys():
          # if no iterations of the reaction name are defined so far, move the already defined name from ReactionName to ReactionName#1
          if '#' not in ReactionName:
            ReacNameCopy = ReactionName + '#1'
            hdf_reac_group.move(ReactionName, ReacNameCopy)
            reac_dict[key] = ReacNameCopy
            # Create the new name and add the reaction to the database
          ReactionName = re.sub('#\d+', '', ReactionName)
          Count = Count + 1
          ReactionName = ReactionName + '#' + str(Count)
        reac_dict[key] = ReactionName
        hdf_reac = ReactionName
        hdf_reac = hdf_reac_group.create_dataset(hdf_reac,data=[0])
        hdf_reac.attrs['* Created']   = date.today().strftime("%B %d, %Y")
        print('The reaction with the new model was stored as ', ReactionName, '.')
        model_name_list = []
        model_name_list.append(model_name)
        hdf_reac.attrs['ChemistryModel'] = np.array(model_name_list,dtype='S255')
        exit
  else:
    # Add the attribute ChemistryModel to the database
    model_name_list = []
    model_name_list.append(model_name)
    hdf_reac.attrs['ChemistryModel'] = np.array(model_name_list,dtype='S255')

# List of variables not added to the file
exclude_list = ['NumberOfNonReactives']
str_list = ['ReactionModel']

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
          # Read-In of the attributes, check if the parameter is already defined
          if var_name[1] in hdf_reac.attrs:
            print('Reaction parameter is already set: ', var_name[1])
          else:
            if is_float(var_value):
              var_value = float(var_value)
            # Exclude certain parameters
            if var_name[1] not in exclude_list:
              # Treatment of the NonInteractiveSpecies (added in the form of the species name)
              if var_name[1] in reac_attr_list:
                spec_name_list = []
                var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
                for val in var_value:
                  spec_name_list.append(spec_dict[val])
                hdf_reac.attrs[var_name[1]] = np.array(spec_name_list,dtype='S255')
              elif var_name[1] in str_list:
                hdf_reac.attrs[var_name[1]] = np.array([var_value],dtype='S255')
              else:
                # All other parameters
                hdf_reac.attrs[var_name[1]] = var_value
                # If not defined, the standard reaction model is set to TCE
              print('Reaction parameter set: ', var_name[1]) 
              if 'ReactionModel' not in hdf_reac.attrs:
                hdf_reac.attrs['ReactionModel'] = np.array(['TCE'],dtype='S255')       
            # Write attributes for source and time of retrieval
            hdf_reac.attrs['* Reference'] = args.reference
            hdf_reac.attrs['* Created']   = date.today().strftime("%B %d, %Y")

# surf_attr_list = ['Reactants', 'Products']
# surf_wo_readin = ['Boundaries', 'NumOfBoundaries', 'Inhibition', 'Promotion']
# surf_string_list = ['Type']

# # Read-in of gas-surface reaction data
# with open(args.ini_filename) as file:
#   surf_dict = {}
#   for line in file:
#     if not line.startswith('!'):
#       if line.startswith('Surface-Reaction'):
#         var_name = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[0].split('-')
#         var_value = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[1].split('!', 1)[0]
#         if var_name[1].startswith('SurfName'):
#           hdf_surf = var_value
#           surf_count = var_name[0]
#           surf_dict[surf_count] = hdf_surf
#           if hdf_surf in hdf_surf_group.keys():
#             print('Surface reaction already exists: ', hdf_surf)
#             hdf_surf = hdf_surf_group[hdf_surf]
#           else:
#             print('Surface reaction added to the database: ', hdf_surf)
#             hdf_surf = hdf_surf_group.create_dataset(hdf_surf,data=[0])
#             hdf_surf.attrs['* Created']   = date.today().strftime("%B %d, %Y")

# with open(args.ini_filename) as file:
#   for line in file:
#     if not line.startswith('!'):
#       if line.startswith('Surface-Reaction'):
#         var_name = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[0].split('-')
#         var_value = line.strip().replace(" ","").replace("Surface-Reaction","").split('=')[1].split('!', 1)[0]

#         if not var_name[1].startswith('SurfName'):
#           surf_count = var_name[0]
#           hdf_surf = surf_dict[surf_count]
#           hdf_surf = hdf_surf_group[hdf_surf]
#           if var_name[1] in hdf_surf.attrs:
#             print('Catalytic parameter is already set: ', var_name[1])
#           else:
#             if is_float(var_value):
#               var_value = float(var_value)
#             if var_name[1] not in surf_attr_list and var_name[1] not in surf_wo_readin:
#               if var_name[1] not in surf_string_list:
#                 hdf_surf.attrs[var_name[1]] = var_value
#               else:
#                 hdf_surf.attrs[var_name[1]] = np.string_(var_value)     
#               print('Catalytic parameter set: ', var_name[1]) 
#             elif var_name[1] in surf_attr_list:
#               spec_name_list = ''
#               var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
#               for val in var_value:
#                 spec_name_list = spec_name_list + ', ' + spec_dict[val]
#               spec_name_list = spec_name_list.replace(', ', '', 1)
#               hdf_surf.attrs[var_name[1]] = spec_name_list
#               print('Catalytic parameter set: ', var_name[1]) 
#             # Write attributes for source and time of retrieval
#             hdf_surf.attrs['* Reference'] = args.reference
#             hdf_surf.attrs['* Created']   = date.today().strftime("%B %d, %Y")
    

print('*** Database successfully built. ***')

h5_species.close()
h5_electronic.close()
h5_crosssection.close()
