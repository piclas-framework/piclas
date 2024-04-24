import h5py
import re
import argparse
from general_functions import *
from edit_species import *
from edit_crosssections import *
from edit_reactions import *
from edit_surfchem import *

###################################################################################################
# - Program starts here
###################################################################################################
# arguments for ionization and species setting
parser = argparse.ArgumentParser(description='Create electronic database')
parser.add_argument('--species', nargs='+', metavar='', help='Select species to add into Species Database (also possible directly in script, line 160). Enter species like H Fe C')
parser.add_argument('--setting', type=str, metavar='', help='Select which setting should be run: "a" for adding new dataset to Species database or "c" for checking existing database for updates')
parser.add_argument('--reactions', nargs='+', metavar='', help='Select reactions to add into Database (also possible directly in script). Enter reactions like C+N_CNIon1+electron C2+M_C+M+C')
args = parser.parse_args()

# Read in Species Database
relative_path = '../../SpeciesDatabase.h5'
hdf_unified_data = h5py.File(relative_path, 'a')
# select ATcT Version (script build for 1.130)
version = '1.130'
ATcT_URL = f'https://atct.anl.gov/Thermochemical%20Data/version%20{version}'
# Base URL of the query
URL_base = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
# flag for reading in data from ATcT only once
species_dict = []
species_dict_flag = False

###################################################################################################
# variables to set in code if not by user input or arguments, leave empty if not used
###################################################################################################
# set reaction list in code like ['C+N_CNIon1+electron','C2+M_C+M+C']
code_reactions = []
code_species = []

# Prompt user to choose which dataset to keep
print("\n\nThis script is meant to edit and maintain the " +yellow("unified species database") + " of PICLas")
user_input = get_valid_input(create_prompt('to maintain/edit species', 'to maintain/edit chemical reactions', 'to maintain/edit cross section data', 'to maintain/edit surface chemistry'), lambda x: x == '1' or x == '2' or x =='3' or x =='4' or x =='5')


###################################################################################################
# SPECIES
###################################################################################################
if user_input == "1":
    user_input = get_valid_input(create_prompt('check existing species', 'add new species'), lambda x: x == '1' or x == '2' or x =='3')
    if user_input == "1":
        print("Checking electronic levels "+red("currently only for atoms")+ " from "+blue(URL_base)+" and from custom csv file for molecules if set\nand attributes for all species from "+ blue(ATcT_URL))
        existing_species_list = list(hdf_unified_data["Species"].keys())
        species_list = []
        # check if species to check was provided
        if args.species == None:
            if len(code_species) == 0:
                # run through all atoms species from species database
                print('Checking all existing atom species in database and species with provided custom csv electronic levels')
                for i, species in enumerate(existing_species_list):
                    if sum(1 for c in species.replace('Ion','') if c.isupper()) == 1:
                        if bool(re.search(r'[A-Za-z]*\d', re.sub(r'Ion\d+','',species))):
                            continue
                        species_list.append(species)
                
                # Directory path
                current_directory = os.getcwd()
                for root, dirs, files in os.walk(current_directory):
                    for filename in files:
                        if 'custom_electronic_levels' in filename and 'SPECIES' not in filename:        # exclude template
                            species_list.append(filename.split('_')[-1].split('.')[0])

            else:
                species_list = code_species
    
        else:
            species_list = args.species
            existing_species_list = args.species
        
        print('Checking species: ' + ', '.join(species_list))
        
        # check electronic levels for atoms and custom data
        for current_species in species_list:
            if check_datasets(current_species, hdf_unified_data, relative_path) == -1:
                break

        # check attributes for all species in database
        for current_species in existing_species_list:
            species_dict_flag, break_flag = check_ATcT(current_species, hdf_unified_data, species_dict_flag, species_dict, ATcT_URL)
            if break_flag == -1:
                break

    elif user_input == "2":
        # check if species list was provided as argument, if not use species list provided in python
        if args.species == None:
            user_input_species = input(bold('\nPlease enter species list as comma separated string') + ', e.g. Fe,Ar,H,CIon1,CIon2,C\n')
            species_list = user_input_species.split(',')
        else:
            species_list = args.species
        
        for current_species in species_list:
            # still fails if species not found in species dict -> edit_species.py line 370
            add_dataset(hdf_unified_data,current_species,relative_path,ATcT_URL,species_dict,species_dict_flag)

    elif user_input == "3":
        own_exit()

###################################################################################################
# REACTIONS
###################################################################################################

elif user_input == "2":
    user_input = get_valid_input(create_prompt('add new reactions', 'delete reactions'), lambda x: x == '1' or x == '2' or x =='3')
    if user_input == '1':
        function = create_reaction
    elif user_input == '2':
        function = delete_reaction
    elif user_input == '3':
     own_exit()
    
    # operation
    if args.reactions == None:
        if len(code_reactions) == 0:
            user_input_reaction = input(bold('\nPlease enter reactions list as comma separated string') + ', e.g. C+N_CNIon1+electron,C2+M_C+M+C\n')
            reaction_list = user_input_reaction.split(',')
        else:
            reaction_list = code_reactions
    else:
        reaction_list = args.reactions

    for reaction in reaction_list:
        function(hdf_unified_data, reaction)

###################################################################################################
# CROSS SECTION DATA
###################################################################################################

elif user_input == "3":
    print("Not implemented yet")
    
###################################################################################################
# SURFACE CHEMISTRY
###################################################################################################

elif user_input == "4":
    print("Not implemented yet")

###################################################################################################
# EXIT
###################################################################################################
elif user_input == "5":
    own_exit()

print('Done.')
hdf_unified_data.close()