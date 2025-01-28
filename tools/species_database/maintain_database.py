import re
import argparse
from pathlib import Path
from general_functions import *
from edit_species import *
from edit_crosssections import *
from edit_reactions import *
from edit_surfchem import *
from edit_diffusion_coefficients import *
from config import *

###################################################################################################
# - Program starts here
###################################################################################################
# arguments for ionization and species setting
parser = argparse.ArgumentParser(description='Create electronic database')
parser.add_argument('--species', nargs='+', metavar='', help='Select species to add into Species Database (also possible directly in script, line 160). Enter species like H Fe C')
args = parser.parse_args()

###################################################################################################
# variables to set in code if not by user input or arguments, leave empty if not used
###################################################################################################
# set reaction list in code like ['C+N_CNIon1+electron','C2+M_C+M+C']
code_reactions = []
###################################################################################################
# Initial prompt
print("\n\nThis script is meant to edit and maintain the " +yellow("unified species database") + " of PICLas")
user_input = get_valid_input(create_prompt('to maintain/edit species',
                                           'to maintain/edit chemical reactions',
                                           'to maintain/edit cross section data',
                                           'to maintain/edit surface chemistry',
                                           'to maintain/edit diffusion coefficients',),
                                           lambda x: x == '1' or x == '2' or x =='3' or x =='4' or x =='5' or x =='6')
###################################################################################################
# SPECIES
###################################################################################################
if user_input == "1":
    # create list of existing species
    existing_species_list = list(hdf_unified_data["Species"].keys())
    # create list of store classes for each species
    list_of_species_classes = []
    # create class for each species
    for species in existing_species_list:
        instance = create_instance_from_data(species)
        list_of_species_classes.append(instance)
    # list contains the species where the excitation states are checked
    check_levels_list = []
    # list contains the species where the attributes are checked
    check_attribute_list = []

    user_input = get_valid_input(create_prompt('check existing species',
                                               'add new data to existing species',
                                               'add new species'),
                                               lambda x: x == '1' or x == '2' or x =='3' or x =='4')
    if user_input == "1":
        # create dictionary from ATcT data for faster access
        species_dict = create_species_dict()
        # check if species to check was provided by argument or in code, if not check all atom species
        if args.species == None:
            print("Checking electronic levels "+red("currently only for atoms")+ " from "+blue(URL_base)+" and from custom csv file for molecules if set\nand attributes for all species from "+ blue(ATcT_URL))
            # run through all atoms species from species database
            print('Checking all existing atom species in database and species with provided custom csv electronic levels')
            for i, species in enumerate(existing_species_list):
                if sum(1 for c in species.replace('Ion','') if c.isupper()) == 1:
                    if bool(re.search(r'[A-Za-z]*\d', re.sub(r'Ion\d+','',species))):
                        continue
                    check_levels_list.append(species)

            # also append species with custom electronic levels
            current_directory = os.getcwd()
            for root, dirs, files in os.walk(current_directory):
                for filename in files:
                    if 'custom_electronic_levels' in filename and 'SPECIES' not in filename:        # exclude template
                        check_levels_list.append(filename.split('_')[-1].split('.')[0])
                    if 'custom_vibrational_levels' in filename and 'SPECIES' not in filename:        # exclude template
                        check_levels_list.append(filename.split('_')[-1].split('.')[0])
            # check attributes for all species
            check_attribute_list = existing_species_list

        else:
            check_levels_list = args.species
            check_attribute_list = args.species

        print('Checking species: ' + ', '.join(check_attribute_list))

        # loop over all classes (equals all species in database)
        for species_class in list_of_species_classes:
            # check electronic levels for atoms and custom data
            if species_class.name in check_levels_list:
                # read in inner energy levels (for atoms only electronic)
                read_datasets_from_existing_species(species_class)
                species_class.add_or_update_all_possible_datasets()

            if species_class.name in check_attribute_list:
                # check attributes for species in database
                # use data from ATcT, see edit_species.py for creation of species_dict dictionary
                atct_data_for_instance={
                    'HeatOfFormation_K':species_dict[species_class.name][0],
                    'MassIC':species_dict[species_class.name][1],
                    'ChargeIC':species_dict[species_class.name][2],
                    'InteractionID':species_dict[species_class.name][3],
                    }
                check_attributes_from_actc(species_class, atct_data_for_instance)
                write_instance_to_database(species_class)


    elif user_input == "2":
        # check if species list was provided as argument, if not use species list provided in python
        if args.species == None:
            user_input_species = input(bold('\nPlease enter species list as comma separated string') + ', e.g. Fe,Ar,H,CIon1,CIon2,C\n')
            check_species_list = user_input_species.split(',')
        else:
            check_species_list = args.species

        print(green('\nAdding data to species: ' + ', '.join(check_species_list)))
        hit = False
        # loop over all classes (equals all species in database)
        for species_class in list_of_species_classes:
            if species_class.name in check_species_list:
                hit = True
                # get inner energy levels and degenracy, function also allows to add new datasets if none are found
                read_datasets_from_existing_species(species_class)
                edit_attributes(species_class)
                write_instance_to_database(species_class)
        if not hit:
            print('No species found in database for given species list!')

    elif user_input == "3":
        # create dictionary from ATcT data to access it faster
        species_dict = create_species_dict()
        # check if species list was provided as argument, if not use species list provided in python
        if args.species == None:
            user_input_species = input(bold('\nPlease enter species list as comma separated string') + ', e.g. Fe,Ar,H,CIon1,CIon2,C\n')
            check_species_list = user_input_species.split(',')
        else:
            check_species_list = args.species

        # loop over all classes (equals all species in database)
        for new_species in check_species_list:
            # skip if species exists in database
            if new_species in existing_species_list:
                print(red('Error:')+' Species ' + new_species + ' is already stored in database! It can only be edited.')
                continue

            # create empty instance for given species which results in right class depening on species name
            species_class = create_empty_instance(new_species)
            # add inner energy levels and degenracy
            species_class.add_or_update_all_possible_datasets()
            # try to add attributes from ATcT data
            try:
                # write all possible attributes for given species
                atct_data_for_instance={
                    'HeatOfFormation_K':species_dict[new_species][0],
                    'MassIC':species_dict[new_species][1],
                    'ChargeIC':species_dict[new_species][2],
                    'InteractionID':species_dict[new_species][3],
                    }
                for attr_name, attr_value in atct_data_for_instance.items():
                    print("Adding attribute " + attr_name + " with value " + str(attr_value))
                    species_class.attributes[attr_name] = attr_value
            except Exception as e:
                # pass since missing data is added via user input anyway
                pass

            # get remaining attributes from user input
            for attr_name, attr_value in species_class.attributes.items():
                if attr_value is None:
                    # Set the attribute to a default value based on its type
                    if attr_name in attribute_types:
                        attr_value = input(f"Please enter a value for {attr_name}: ")
                        species_class.attributes[attr_name] = attribute_types[attr_name](attr_value)
                print("Adding attribute " + attr_name + ": " + str(attr_value))
            write_instance_to_database(species_class)

    elif user_input == "4":
        own_exit()

###################################################################################################
# REACTIONS
###################################################################################################

elif user_input == "2":
    user_input = get_valid_input(create_prompt('add new reactions',
                                               'delete reactions'),
                                               lambda x: x == '1' or x == '2' or x =='3')
    if user_input == '1':
        function = create_reaction
    elif user_input == '2':
        function = delete_reaction
    elif user_input == '3':
     own_exit()

    # operation
    if len(code_reactions) == 0:
        user_input_reaction = input(bold('\nPlease enter reactions list as comma separated string') + ', e.g. C+N_CNIon1+electron,C2+M_C+M+C\n')
        reaction_list = user_input_reaction.split(',')
    else:
        reaction_list = code_reactions

    for reaction in reaction_list:
        function(reaction)

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
# Diffusion coefficients
###################################################################################################

elif user_input == "5":
    txt_files = [f for f in os.listdir('.') if f.endswith('.txt')]
    if txt_files == []:
        # check templates directory
        txt_files = [f for f in os.listdir('templates') if f.endswith('.txt')]
        for i, file in enumerate(txt_files):
            txt_files[i] = 'templates/' + file
        # exclude template files
        # if txt_files == ['swarm_data_EXAMPLE.txt']:
        if False:
            dir_input = input('No .txt files found in current directory and templates directory (except swarm_data_EXAMPLE.txt). Please supply the path to the .txt files with diffusion coefficients.\n').strip()
            path = Path(dir_input).expanduser()
            if not path.is_absolute():
                # Try the path as-is first
                if not path.exists():
                    # Only resolve if the direct path doesn't exist
                    path = path.resolve()
                    if not path.exists():
                        s = f"Path '{path}' does not exist"
                        raise ValueError(s)
                    if not path.is_dir():
                        s = f"'{path}' is not a directory"
                        raise ValueError(s)
                    if not any(f.suffix == '.txt' for f in path.iterdir()):
                        s = f'No .txt files found in {path}'
                        raise ValueError(s)
            txt_files = [f for f in os.listdir(path) if f.endswith('.txt')]
    print('\n------------------------------------------------------------------------')
    # Display the .txt files with numbers
    for i, file in enumerate(txt_files):
        print(f" [{i + 1}] - {file}")
    print('------------------------------------------------------------------------')
    InputFILES = read_file_by_numbers(txt_files)
    for file in InputFILES:
        bolsig_bool = check_for_bolsig(file)
        if bolsig_bool == True:
            append_to_database_bolsig(relative_path, file)
        else:
            append_to_database(relative_path, file)

###################################################################################################
# EXIT
###################################################################################################
elif user_input == "5":
    own_exit()

print('Done.')
hdf_unified_data.close()