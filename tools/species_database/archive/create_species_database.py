from cmath import log
from matplotlib.rcsetup import validate_bool
import numpy as np
import h5py
from argparse import ArgumentParser
from datetime import date
import re
import sys
from functions_database import custom_sort_reactants, custom_sort_products

parser = ArgumentParser(prog='create_species_database')
parser.add_argument("-p", "--parameter", dest="ini_filename",
                    help="DSMC.ini file to read-in parameters", metavar="FILE", default="DSMC.ini")
parser.add_argument("-e", "--electronic", dest="database_electronic",
                    help="Electronic excitation database to read-in parameters", metavar="FILE", default="")
parser.add_argument("-c", "--crosssection", dest="database_crosssection",
                    help="Crosssection database to read-in parameters", metavar="FILE", default="")
parser.add_argument("-o", "--output", dest="database_output",
                    help="Output file name", metavar="FILE", default="SpeciesDatabase.h5")
parser.add_argument("-r", "--reference", dest="reference",
                    help="Reference for species data", default="Unknown")
parser.add_argument("-s", "--surface", dest="database_surf",
                    help="Sticking coefficient data", metavar="FILE", default="")
# parser.add_argument("-s", "--radiation", dest="rad_file_name",
#                     help="Radiation file to read-in parameters", metavar="FILE", default="Rad.dat")

args = parser.parse_args()

h5_species      = h5py.File(args.database_output, 'a')
if args.database_electronic != "":
  h5_electronic   = h5py.File(args.database_electronic, 'r')
else:
  print('No electronic state data is given as input.')
if args.database_crosssection != "":
  h5_crosssection = h5py.File(args.database_crosssection, 'r')
else:
  print('No crosssection data is given as input.')
if args.database_surf != "":
  h5_surface = h5py.File(args.database_surf, 'r')
else:
  print('No surface data is given as input.')

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

hdf_surf_group = 'Surface-Chemistry'
if hdf_surf_group in h5_species.keys():
  print('Group Surface-Chemistry already exists.')
  hdf_surf_group = h5_species['Surface-Chemistry']
  hdf_surf_group.attrs['* Last Modified'] = date.today().strftime("%B %d, %Y")
else:
  print('Group Surface-Chemistry does not exist, creating new group.')
  #hdf_surf_group = h5_species.create_group('Surface-Chemistry')
  #hdf_surf_group.attrs['* Created'] = date.today().strftime("%B %d, %Y")
  # Copy sticking coefficient data if not defined already
  if args.database_surf != "":
    for dataset in h5_surface.keys():
      print('Surface-Chemistry added: ', dataset)
      h5_species.copy(source=h5_surface[dataset],dest=h5_species)

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
            if args.database_electronic != "":
              if hdf_species in h5_electronic.keys():
                dset = np.array(hdf_species_group.get(hdf_species))
                if len(dset) == 1:
                  del hdf_species_group[hdf_species]
                  hdf_input_data = h5_electronic[hdf_species]
                  print('Electronic states added to the database: ', hdf_species)
                  hdf_species = hdf_species_group.create_dataset(hdf_species,data=hdf_input_data)
                  hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")
              else: 
                hdf_species = hdf_species_group[hdf_species]
          elif hdf_species == 'electron':
            print('Species added to the database: ', hdf_species)
            hdf_species = hdf_species_group.create_dataset(hdf_species,data=[0])
            hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")
          else:
            # Add the electronic state data
            if args.database_electronic != "":
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
          if species_count not in spec_dict:
            print('Error: the species', species_count, 'has no species name in the input file.')
            sys.exit()
          hdf_species = spec_dict[species_count]
          hdf_species = hdf_species_group[hdf_species]
          # Check if the species parameter is already defiend
          if var_name[1] not in hdf_species.attrs:
            # Float conversion
            if is_float(var_value):
              var_value = float(var_value)
            if var_name[1] not in logical_list and var_name[1] not in spec_attr_list:
              hdf_species.attrs[var_name[1]] = var_value
              # Previous-state read in as a string
            elif var_name[1] in spec_attr_list:
              var_value = str(int(var_value))
              spec_name_list = spec_dict[var_value]
              hdf_species.attrs[var_name[1]] = np.string_(spec_name_list)
              # Logicals
            else:
              if 'F' in var_value or 'false' in var_value:
                var_value = 0
              else:
                var_value = 1
              hdf_species.attrs[var_name[1]] = var_value
            # Write attributes for source and time of retrieval
            hdf_species.attrs['* Reference'] = args.reference
            hdf_species.attrs['* Created']   = date.today().strftime("%B %d, %Y")

# Copy cross-section data if not defined already
if args.database_crosssection != "":
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
      if val not in spec_dict:
        print('Error: the reactant ', val, 'is not defined as a species.')
        sys.exit()
      ReactionName = ReactionName + '+' + spec_dict[val]
    ReactionName = ReactionName + '_'
    for val in product_dict[key]:
      if val not in spec_dict:
        print('Error: the product ', val, 'is not defined as a species.')
        sys.exit()
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
          if var_name[1] not in hdf_reac.attrs:
            if is_float(var_value):
              var_value = float(var_value)
            # Exclude certain parameters
            if var_name[1] not in exclude_list:
              # Treatment of the NonInteractiveSpecies (added in the form of the species name)
              if var_name[1] in reac_attr_list:
                spec_name_list = []
                var_value = var_value.replace(',0', '').replace('(/', '').replace('/)', '').split(',')
                for val in var_value:
                  if val not in spec_dict:
                    print('Error: the species ', val, 'is not defined in the input file.')
                    sys.exit()
                  spec_name_list.append(spec_dict[val])
                hdf_reac.attrs[var_name[1]] = np.array(spec_name_list,dtype='S255')
              elif var_name[1] in str_list:
                hdf_reac.attrs[var_name[1]] = np.array([var_value],dtype='S255')
              else:
                # All other parameters
                hdf_reac.attrs[var_name[1]] = var_value
                # If not defined, the standard reaction model is set to TCE
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
if args.database_electronic != "":
  h5_electronic.close()
if args.database_crosssection != "":
  h5_crosssection.close()


#--------------------------------------- combine reactions ---------------------------------------#
###################################################################################################
#----------------------------------------- quick changes -----------------------------------------#
# Disclaimer: script only works on 'old' database from a27fcdae1f4d312a7f629273dc17ed8eade95769
###################################################################################################
# Name of Attribute of the Reaction #
reaction_attribute_name = 'ChemistryModel,NonReactiveSpecies'
# this attribute is an array with dimension:    (X,1), where in the first column the Name of the Chemisty Model is stored 
#                                               or
#                                               (X,2), where in the first column the Name of the Chemisty Model is stored and in the second column the nonreactive species as a string seperated by ','
###################################################################################################
# Name of Attribute of Reactions (the group containing all reactions) #
all_reactions_attribute_name = 'AvailableChemistryModels,NonReactiveSpecies,Reference'
###################################################################################################
# relative tolerance for arrhenius comparison to check if two reactions should be put together (keeping the firsst arrhenius coefficients)
Rtol = 0.003
###################################################################################################
# relative path to species database
relative_path = '../../SpeciesDatabase.h5'
###################################################################################################

import h5py
import re
import numpy as np

def create_dataset(hdf_unified_data, chem_model_list, non_reac_species_list, reference_list):
    # array_of_strings is used for the single reactions to combine chemisty models and Full_Chem_Model_List is for saving the Models for the attribute of the reactions group
    if non_reac_species_list != []:
        array_of_strings = np.empty((len(chem_model_list),2), dtype='S255')
        for i in range(len(chem_model_list)):

            Full_Chem_Model_List.append(chem_model_list[i][0])
            Full_Ref_List.append(str(reference_list[i]).encode('ascii', 'ignore').decode('utf-8'))

            array_of_strings[i,0] = chem_model_list[i][0]

            non_reac_species_str = ''
            for species in non_reac_species_list[i]:
                non_reac_species_str += species + ','

            non_reac_species_str = non_reac_species_str[:-1]
            array_of_strings[i,1] = non_reac_species_str
            Full_Non_Reac_List.append(non_reac_species_str)
    
    else:
        array_of_strings = np.empty((len(chem_model_list),1), dtype='S255')
        for i in range(len(chem_model_list)):
            Full_Chem_Model_List.append(chem_model_list[i][0])
            Full_Ref_List.append(str(reference_list[i]).encode('ascii', 'ignore').decode('utf-8'))
            Full_Non_Reac_List.append('-')
            array_of_strings[i,0] = chem_model_list[i][0]

    # delete old attributes
    hdf_unified_data['Reactions'][combined_reaction_name].attrs.create(reaction_attribute_name, array_of_strings)
    if '* Reference' in hdf_unified_data["Reactions"][combined_reaction_name].attrs:
        del hdf_unified_data["Reactions"][combined_reaction_name].attrs['* Reference']
    if 'NonReactiveSpecies' in hdf_unified_data["Reactions"][combined_reaction_name].attrs:
        del hdf_unified_data["Reactions"][combined_reaction_name].attrs['NonReactiveSpecies']
    if 'ChemistryModel' in hdf_unified_data["Reactions"][combined_reaction_name].attrs:
        del hdf_unified_data["Reactions"][combined_reaction_name].attrs['ChemistryModel']


# Read in Species Database
hdf_unified_data = h5py.File(relative_path, 'a')

# get list of all reactions
reactions_list = list(hdf_unified_data["Reactions"].keys())

Full_Chem_Model_List = []
Full_Ref_List = []
Full_Non_Reac_List = []

# Loop over species list
for i,reaction in enumerate(reactions_list):
    
    # create new reaction dataset if current name differs from name before
    if (re.sub(r'#\d+','',reaction) != re.sub(r'#\d+','',reactions_list[i-1])) or reaction==reactions_list[0]:
        # set counter for new combined name later if Arrhenius changes (if not reaction with different arrhenius coefficientes would be stored under the same name)
        j = 0
        # set name for reaction
        combined_reaction_name = re.sub(r'#\d+','',reaction)

        # Store existing attributes
        attrs = hdf_unified_data["Reactions"][reaction].attrs

        # create lists for created reaction
        chem_model_list = []
        non_reac_species_list = []
        reference_list = []

        # check its the first iteration and the reaction name exists already
        if combined_reaction_name == reaction:
            combined_reaction_name = reaction
            dataset = hdf_unified_data["Reactions"][reaction]
        else:
            dataset = hdf_unified_data["Reactions"].create_dataset(combined_reaction_name, data=0)
        
        #  Restore attributes
        for attr_name, attr_value in attrs.items():
            dataset.attrs[attr_name] = attr_value
    
    # enter loop to check arr coefficients if the current name is equal to the name before (after 'or' just to insure loop is entered even with new changed name from j counter)
    elif re.sub(r'#\d+','',reaction) == combined_reaction_name or re.sub(r'#\d+','',reaction) == re.sub(r'#0\d+','',combined_reaction_name):
        # Check if both reactions have the necessary attributes before comparing them
        if "Arrhenius-Powerfactor" in hdf_unified_data["Reactions"][reaction].attrs and \
            "Arrhenius-Prefactor" in hdf_unified_data["Reactions"][reaction].attrs and \
            "Activation-Energy_K" in hdf_unified_data["Reactions"][reaction].attrs:

            # check the three Arrhenius coefficients to decide if new reaction has to be created
            if (not np.isclose(hdf_unified_data["Reactions"][reaction].attrs['Arrhenius-Powerfactor'], Arrhenius_Powerfactor, rtol=Rtol, atol=0.0)) or \
            (not np.isclose(hdf_unified_data["Reactions"][reaction].attrs['Arrhenius-Prefactor'], Arrhenius_Prefactor, rtol=Rtol, atol=0.0)) or \
            (not np.isclose(hdf_unified_data["Reactions"][reaction].attrs['Activation-Energy_K'], Activation_Energy_Kr, rtol=Rtol, atol=0.0)) or \
            (not (hdf_unified_data["Reactions"][reaction].attrs['ReactionModel'][0] == ReactionModel)):

                # add old chem_models to attribute to create new reaction before resetting lists to store new data
                create_dataset(hdf_unified_data, chem_model_list, non_reac_species_list, reference_list)
                # counter plus one to ensure other name of reaction
                j += 1
                # empty lists for new reaction
                chem_model_list = []
                non_reac_species_list = []
                reference_list = []

                # Store existing attributes
                attrs = hdf_unified_data["Reactions"][reaction].attrs 
                combined_reaction_name = re.sub(r'#\d+','',reaction) + '#0%i' % j
                dataset = hdf_unified_data["Reactions"].create_dataset(combined_reaction_name, data=0)
                #  Restore attributes
                for attr_name, attr_value in attrs.items():
                    dataset.attrs[attr_name] = attr_value

        elif ("Arrhenius-Powerfactor" not in hdf_unified_data["Reactions"][reaction].attrs or \
            "Arrhenius-Prefactor"  not in hdf_unified_data["Reactions"][reaction].attrs or \
            "Activation-Energy_K" not in hdf_unified_data["Reactions"][reaction].attrs) and \
            (Arrhenius_Powerfactor!=0 or Arrhenius_Prefactor!=0 or Activation_Energy_Kr!=0):
            # counter plus one to ensure other name of reaction
            j += 1
            # empty lists for new reaction
            chem_model_list = []
            non_reac_species_list = []
            reference_list = []

            # Store existing attributes
            attrs = hdf_unified_data["Reactions"][reaction].attrs 
            combined_reaction_name = re.sub(r'#\d+','',reaction) + '#0%i' % j
            dataset = hdf_unified_data["Reactions"].create_dataset(combined_reaction_name, data=0)
            #  Restore attributes
            for attr_name, attr_value in attrs.items():
                dataset.attrs[attr_name] = attr_value
    
    # store values for loop in next iter since reaction will be deleted before the if loop
    ReactionModel = hdf_unified_data["Reactions"][reaction].attrs['ReactionModel'][0]
    if "Arrhenius-Powerfactor" in hdf_unified_data["Reactions"][reaction].attrs and \
        "Arrhenius-Prefactor" in hdf_unified_data["Reactions"][reaction].attrs and \
        "Activation-Energy_K" in hdf_unified_data["Reactions"][reaction].attrs:
        
        Activation_Energy_Kr = hdf_unified_data["Reactions"][reaction].attrs['Activation-Energy_K']
        Arrhenius_Prefactor= hdf_unified_data["Reactions"][reaction].attrs['Arrhenius-Prefactor']
        Arrhenius_Powerfactor = hdf_unified_data["Reactions"][reaction].attrs['Arrhenius-Powerfactor']
    else:
        Activation_Energy_Kr    = 0
        Arrhenius_Prefactor     = 0
        Arrhenius_Powerfactor   = 0

    # if current reaction has the same name as the combined reaction and arr are equal since elif was skipped -> add the data to the lists (after 'or' just to insure loop is entered even with new changed name from j counter)
    if re.sub(r'#\d+','',reaction) == combined_reaction_name or re.sub(r'#\d+','',reaction) == re.sub(r'#0\d+','',combined_reaction_name):
        # Store existing attributes
        attrs = hdf_unified_data["Reactions"][reaction].attrs

        chem_model_list.append(attrs['ChemistryModel'])
        reference_list.append(attrs['* Reference'])

        if 'NonReactiveSpecies' in hdf_unified_data["Reactions"][reaction].attrs:
            # get strings in byte format
            byte_strings = attrs['NonReactiveSpecies']
            # Convert each byte string to a regular string
            species_list = [byte_string.decode('utf-8') for byte_string in byte_strings]
            non_reac_species_list.append(species_list)

        # delete old reaction
        if combined_reaction_name != reaction:
            del hdf_unified_data["Reactions"][reaction]

    # if last reaction in list or current reaction is not equal to next reaction (after second 'or' just to insure loop is entered even with new changed name from j counter)
    if reaction == reactions_list[-1] or re.sub(r'#\d+','',reaction) != re.sub(r'#\d+','',reactions_list[i+1]) or re.sub(r'#\d+','',reaction) != re.sub(r'#0\d+','',reactions_list[i+1]):
        
        create_dataset(hdf_unified_data, chem_model_list, non_reac_species_list, reference_list)

# convert Full lists to arrays to store in attribute
Full_Non_Reac_array = np.array(Full_Non_Reac_List, dtype='S255')
Full_Chem_Model_array = np.array(Full_Chem_Model_List, dtype='S255')
Full_Ref_array = np.array(Full_Ref_List, dtype='S255')
Full_Array = np.vstack((Full_Chem_Model_array, Full_Non_Reac_array, Full_Ref_array))
Full_Array = Full_Array.T
unique_array = np.unique(Full_Array, axis=0)
hdf_unified_data['Reactions'].attrs.create(all_reactions_attribute_name, unique_array)


###################################################################################################
#----------------------------------------- second part -----------------------------------------#
###################################################################################################
# checking all new reactions if arr coefficients are equal and reactions werent adjacent since this wasnt considered before

# Read in Species Database
hdf_unified_data = h5py.File(relative_path, 'a')

# get list of all reactions )stripped from '#01',etc.
unique_reaction_list = []
reaction_list = list(hdf_unified_data["Reactions"].keys())
for item in hdf_unified_data["Reactions"].keys():
    if re.sub(r'#0\d+','',item) not in unique_reaction_list:
        unique_reaction_list.append(re.sub(r'#0\d+','',item))
# counter for reactions in full reaction list (with '#01')
j = 0
for i,reaction in enumerate(unique_reaction_list):

    Activation_Energy_Kr = []
    Arrhenius_Prefactor = []
    Arrhenius_Powerfactor = []
    Reaction_Model = []
    Reaction = []

    # appending all coefficients and reaction names as long as the unique name is equal to the full name when removed '#01'
    while j < len(reaction_list) and reaction == re.sub(r'#0\d+','',reaction_list[j]):
        if "Arrhenius-Powerfactor" in hdf_unified_data["Reactions"][reaction_list[j]].attrs and \
            "Arrhenius-Prefactor" in hdf_unified_data["Reactions"][reaction_list[j]].attrs and \
            "Activation-Energy_K" in hdf_unified_data["Reactions"][reaction_list[j]].attrs:

            Activation_Energy_Kr.append(hdf_unified_data["Reactions"][reaction_list[j]].attrs['Activation-Energy_K'])
            Arrhenius_Prefactor.append(hdf_unified_data["Reactions"][reaction_list[j]].attrs['Arrhenius-Prefactor'])
            Arrhenius_Powerfactor.append(hdf_unified_data["Reactions"][reaction_list[j]].attrs['Arrhenius-Powerfactor'])
            Reaction_Model.append(hdf_unified_data["Reactions"][reaction_list[j]].attrs['ReactionModel'][0])
            Reaction.append(reaction_list[j])
        else:
            Activation_Energy_Kr.append('-')
            Arrhenius_Prefactor.append('-')
            Arrhenius_Powerfactor.append('-')
            Reaction_Model.append(hdf_unified_data["Reactions"][reaction_list[j]].attrs['ReactionModel'][0])
            Reaction.append(reaction_list[j])
        j += 1

    Activation_Energy_Array = np.array(Activation_Energy_Kr)
    Arrhenius_Prefactor_Array = np.array(Arrhenius_Prefactor)
    Arrhenius_Powerfactor_Array = np.array(Arrhenius_Powerfactor)
    Reaction_Model_Array = np.array(Reaction_Model)
    Reaction_Array = np.array(Reaction)
    Full_Array = np.vstack((Reaction_Array,Arrhenius_Powerfactor_Array,Arrhenius_Prefactor_Array,Activation_Energy_Array, Reaction_Model_Array)).T
   # Find rows where the last three columns match
    matching_rows = []
    for i in range(len(Full_Array)):
        for k in range(i + 1, len(Full_Array)):
            if all(Full_Array[i][-4:] == Full_Array[k][-4:]):
                matching_rows.append((i, k))
                
    if matching_rows != 0:
        for match in matching_rows:
            # get matching idices
            index1, index2 = match
            # Store existing attributes
            attrs = hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs

            current_models_array = hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[reaction_attribute_name]
            append_models_array = hdf_unified_data["Reactions"][Full_Array[index2,0]].attrs[reaction_attribute_name]

            # stacked array for attribute creation
            new_models_array = np.vstack((current_models_array,append_models_array))

            # keeping first match so only deleting attribute
            del hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[reaction_attribute_name]
            # deleting whole reaction of second match
            del hdf_unified_data["Reactions"][Full_Array[index2,0]]
            
            hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs.create(reaction_attribute_name, new_models_array)

print('Output successful.')
hdf_unified_data.close()