import requests
import pandas as pd
import io
import re
import numpy as np
import os
from tabulate import tabulate
from datetime import date
from requests_html import HTMLSession
from bs4 import BeautifulSoup
from general_functions import *
from collections import defaultdict
from config import *

# General workflow:
# Different species are handled by different classes - Atom, DiatomicMolecule, PolyatomicMolecule
# This is done to define the right attributes and methods for different types of species since different types of attributes are needed, e.g. characteristic temperatures for molecules
# If all attributes and inner energy levels are set for a species, the instance of the class is written to the database
#
# To extend the functionality and to handle new attributes, simply add them to the attributes dictionary in the class file and initialize them with None,
# so they will be created and need to be set if a new species is created
#
# If the value can be determined automatically, the value can be set in the __init__ function of the class, e.g. ChargeIC and InteractionID
# The key of the dictionary will be the name of the attribute in the database!!!
# For already existing species with existing attributes the __init__ functions will just loop over all attributes so the addition must only be made for the case of creating and adding a new species
# In addition the type of the new attribute must be added to the attribute_types dictionary in the config.py file or needs to be set if a new attribute is added using the script
# Generally the attributes dictionary should be used to store all attributes which are written to the database
# Other variables of the class like the name of the species or the electronic levels are be stored separately to simplify the handling of a different number of attributes in loops
# The new references are also stored in a separate dictionary so they can be added to the existing references

######################################################################################################################################################################################################
# The classes are defined in the following code snippet
######################################################################################################################################################################################################
class Atom:
    def __init__(self, species, read_species=True, ChargeIC=0.0, InteractionID=1.0):
        # read_species: if True, read species from database, if False, create new species and write values later
        # create dictionary to store attributes with keys corresponding to the right name in the database
        self.attributes = {}
        # storing the species name and electronic levels separately (not into the attributes dict) since it should not be written to the database
        self.name = species
        self.new_references = {}
        # read species from database into class variables
        if read_species:
            attrs = hdf_unified_data["Species"][species].attrs
            # Store attributes in class variables
            for attr_name, attr_value in attrs.items():
                self.attributes[attr_name] = attr_value
        # create new class to store new data
        else:
            # initialize the attributes of the class with the configuration from the config file
            for attr_name, attr_value in atom_attributes.items():
                self.attributes[attr_name] = attr_value
            if 'ion' in self.name.lower():
                Ion_Number = int(re.search(r'\d+$', self.name).group())
                self.attributes['ChargeIC'] = -Ion_Number * electron_charge
                self.attributes['InteractionID'] = 10.0
            else:
                self.attributes['ChargeIC'] = ChargeIC
                self.attributes['InteractionID'] = InteractionID
            self.ElectronicLevels = np.array([])

    def add_or_update_all_possible_datasets(self):
        """
        Add or Check the electronic levels of the species (since this class is only for atoms)

        If the electronic levels are already stored in the database, the user is prompted to choose which dataset to keep.
        Otherwise a new dataset is created.
        """
        # read in data to compare existing data or set new data
        try: # try to get electronic data from the nist database
            # Get the ionization number as the last digits and 1 to comply with the NIST standard
            ion_level = int(re.sub('.*?([0-9]*)$',r'\1',self.name) or 0) + 1

            # get data from Nist database
            elec_levels_nist, current_species_NIST = get_data_from_NIST(self.name, ion_level)
            if type(elec_levels_nist) == int:
                return

            # convert data from website to fit specifications (e.g. filter unnecessary information )
            elec_levels_nist = convert_electronic_data(elec_levels_nist)
            # data from URL into numpy array
            elec_levels_nist_array = elec_levels_nist.to_numpy()
            # set reference
            data_origin = 'Levels and Degeneracy - Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2020). NIST Atomic Spectra Database (ver. 5.8), [Online]. Available: https://physics.nist.gov/asd. National Institute of Standards and Technology, Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F. Retrieved on ' + date.today().strftime("%B %d, %Y") + '.'

        except Exception as e:
            print(red('Error'), "reading data from the NIST database for ", green(self.name), ":", e)
            user_input_csv = input(bold('\nPlease enter the csv file name containing the electronic levels: ') )
            custom_dataset = pd.read_csv(f'{user_input_csv}.csv',skipinitialspace=True,delimiter=",",usecols=['J','Levelcm-1','Term'], na_values=['---'], dtype=str)
            elec_levels_nist_array = custom_dataset.astype(float)
            print(f'Using custom data from {user_input_csv}!')
            # set reference
            data_origin = 'user input from' + date.today().strftime("%B %d, %Y") + '.'

        # add new reference
        self.new_references['ElectronicLevels and degeneracy'] = data_origin
        # Check if the electronic levels were read in from the database or new species is created
        if self.ElectronicLevels.size == 0:
            print("Creating electronic levels for species ", green(self.name), "\n")
            self.ElectronicLevels = elec_levels_nist_array
        # electronic levels are stored in the class and will be compared with other data to select which to keep
        else:
            elec_levels_database = self.ElectronicLevels
            # sanity check to see if read-in data from species database has wrong format (might be empty)
            if elec_levels_database.ndim == 1 or elec_levels_database.shape[1] != 2:
                print(red('Error'), ": Read-in from ", cyan(os.path.abspath(relative_path))," failed or dataset has unexpected dimensions (needs to be (X,2)! Please check the dataset of species ", green(self.name), " - species will be skipped\n")
                return
            # Compare the arrays to check if they are already the same
            if np.array_equal(elec_levels_nist_array, elec_levels_database):
                print("\nDatasets of species "+ green(self.name) + " are equal, so dataset from " + yellow('unified species database') + " will be kept.\n")
                return

            # get user input to determine which dataset should be saved
            print("It seems like the species " + green(self.name) + " is already stored in the SpeciesDatabase at " + cyan(os.path.abspath(os.path.abspath(relative_path)))+' with the following data:\n')
            print_diffs(elec_levels_nist_array, elec_levels_database, 'Electronic')

            # Prompt user to choose which dataset to keep
            user_input = get_valid_input(create_prompt('to keep data from '+yellow('unified species database'),
                                                       'to save only electronic level data from '+blue(URL_base)),
                                                       lambda x: x == '1' or x == '2' or x == '3')
            if user_input == '1':
                print("Keeping electronic level dataset and attributes for species ",green(self.name), "\n")
                return
            elif user_input == '2':
                self.ElectronicLevels = elec_levels_nist_array
                print("Saving electronic level dataset from " + blue(data_origin) + " for species ",green(self.name), " but keeping attributes\n")
            elif user_input == '3':
                print(bold(red("Exiting program")))
                exit(1)

    def __str__(self):
        return f"Atom(name={self.name}, attributes={self.attributes})"

# Functions only needed for molecule classes
class MoleculeHelper:
    """Helper class for molecules for shared methods"""
    def add_or_update_dataset(instance, level_type, create_new=False):
        """
        Check the selected levels of the species in the unified species database and the NIST database and prompt the user to choose which dataset to keep
        Currently only supports custom csv data for molecules

        Parameters:
        instance: Instance of the class
        level_type: Type of levels to add or update (Rotational, Vibrational, Electronic)
        create_new: If True, create a new dataset, if False, compare the existing dataset with the new dataset
        """
        if not create_new:
            if level_type.lower() == 'rotational':
                level_type == 'Rotational'
                levels_database = instance.RotationalLevels
            elif level_type.lower() == 'vibrational':
                level_type == 'Vibrational'
                levels_database = instance.VibrationalLevels
            else:
                level_type = 'Electronic'
                levels_database = instance.ElectronicLevels
        else:
            # set up dummy
            levels_database = np.array([])
        # try to read in data from csv file needed to compare and/or store
        try:
            # Attempt to convert data to pandas dataframe from custom csv file
            if level_type == 'Rotational':
                custom_dataset = pd.read_csv(f'custom_rotational_levels_{instance.name}.csv',skipinitialspace=True,delimiter=",",usecols=['J','RotLevelJ'], na_values=['---'], dtype=str)
                custom_dataset = custom_dataset.to_numpy()
                print(f'Using custom data from custom_rotational_levels_{instance.name}.csv')
            elif level_type == 'Vibrational':
                custom_dataset = pd.read_csv(f'custom_rotational_levels_{instance.name}.csv',skipinitialspace=True,delimiter=",",usecols=['VibLevelJ'], na_values=['---'], dtype=str)
                custom_dataset = custom_dataset.astype(float)
                print(f'Using custom data from custom_vibrational_levels_{instance.name}.csv')
            else:
                try:
                    # Convert custom electronic data to format which fits nist data
                    custom_dataset = pd.read_csv(f'custom_electronic_levels_{instance.name}.csv',skipinitialspace=True,delimiter=",",usecols=['J','Levelcm-1'], na_values=['---'], dtype=str)
                    custom_dataset.insert(0, 'Term', '---')
                    last_row = custom_dataset.iloc[-1].copy()
                    last_row['Term'] = 'Limit'
                    last_row['J'] = '---'
                    custom_dataset.iloc[-1] = last_row
                    print(f'Using custom data from custom_electronic_levels_{instance.name}.csv')
                    custom_dataset = convert_electronic_data(custom_dataset)
                    custom_dataset = custom_dataset.to_numpy()
                except:
                    print(red('Caution')+": Please make sure you provide a atom species. Automatic read in from online database of electronic levels for molecules not possible yet! Current species: "+green(instance.name)+" It is possible to provide a custom csv file containing the electronic levels. To use this option please use the template in this directory, insert the necessary information and add the species to the end like this: 'custom_electronic_levels_H2.csv'")

        except Exception as e:
            print(red('Error')+e+"\nCurrent species: "+green(instance.name))
            return

        instance.new_references[f'{level_type}Levels and degeneracy'] = 'user input from' + date.today().strftime("%B %d, %Y") + '.'
        # Check if the levels are already stored in the database/were read in to the class
        if levels_database.size == 0:
            print(f"Creating {level_type} levels for species ", green(instance.name), "\n")
            if level_type == 'Rotational':
                instance.RotationalLevels  = custom_dataset
            elif level_type == 'Vibrational':
                instance.VibrationalLevels = custom_dataset
            else:
                instance.ElectronicLevels  = custom_dataset
        # electronic levels are stored in the class and will be compared with other data to select which to keep
        else:
            # check if read-in data from species database has wrong format (might be empty)
            if levels_database.ndim == 1 or levels_database.shape[1] != 2:
                print(red('Error'), ": Read-in from ", cyan(os.path.abspath(relative_path))," failed or dataset has unexpected dimensions (needs to be (X,2)! Please check the dataset of species ", green(instance.name), " - species will be skipped\n")
                return
            # Compare the arrays to check if they are already the same
            if np.array_equal(custom_dataset, levels_database):
                print("\nDatasets of species "+ green(instance.name) + " are equal, so dataset from " + yellow('unified species database') + " will be kept.\n")
                return

            # get user input to determine which dataset should be saved
            print("It seems like the species " + green(instance.name) + " is already stored in the SpeciesDatabase at " + cyan(os.path.abspath(os.path.abspath(relative_path)))+' with the following data:\n')
            print_diffs(custom_dataset, levels_database, level_type)

            # Prompt user to choose which dataset to keep
            user_input = get_valid_input(create_prompt('to keep data from '+yellow('unified species database'),
                                                       f'to save only {level_type} level data from '+blue(f'custom_{level_type}_levels_{instance.name}.csv')),
                                                       lambda x: x == '1' or x == '2' or x == '3')
            if user_input == '1':
                print(f"Keeping {level_type} level dataset and attributes for species ",green(instance.name), "\n")
                return
            elif user_input == '2':
                if level_type == 'Rotational':
                    instance.RotationalLevels  = custom_dataset
                elif level_type == 'Vibrational':
                    instance.VibrationalLevels = custom_dataset
                else:
                    instance.ElectronicLevels  = custom_dataset
                print(f"Saving {level_type} level dataset from " + blue(f'custom_{level_type}_levels_{instance.name}.csv') + " for species ",green(instance.name), "\n")
            elif user_input == '3':
                print(bold(red("Exiting program")))
                exit(1)

    def get_num_of_atoms(species):
        """Returns number of atoms in a polyatomic molecule"""
        # Regular expression to match elements and their counts
        pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
        matches = pattern.findall(species)

        # Dictionary to store the count of each element
        element_counts = defaultdict(int)

        for element, count in matches:
            if count == '':
                count = 1
            else:
                count = int(count)
            element_counts[element] += count

        # Sum the counts to get the total number of atoms
        total_atoms = sum(element_counts.values())
        return total_atoms

class DiatomicMolecule:
    def __init__(self, species, read_species=True, ChargeIC=0.0, InteractionID=2.0):
        # read_species: if True, read species from database, if False, create new species and write values later
        # create dictionary to store attributes with keys corresponding to the right name in the database
        self.attributes = {}
        # storing the species name and electronic levels separately (not into the attributes dict) since it should not be written to the database
        self.name = species
        self.new_references = {}
        # read species from database into class variables
        if read_species:
            attrs = hdf_unified_data["Species"][species].attrs
            # Store attributes in class variables
            for attr_name, attr_value in attrs.items():
                self.attributes[attr_name] = attr_value
        # create new class to store new data
        else:
            # initialize the attributes of the class with the configuration from the config file
            for attr_name, attr_value in diatomic_attributes.items():
                self.attributes[attr_name] = attr_value
            if 'ion' in self.name.lower():
                Ion_Number = int(re.search(r'\d+$', self.attributes['name']).group())
                self.attributes['ChargeIC'] = -Ion_Number * electron_charge
                self.attributes['InteractionID'] = 20.0
            else:
                self.attributes['ChargeIC'] = ChargeIC
                self.attributes['InteractionID'] = InteractionID
            self.ElectronicLevels = np.array([])
            self.RotationalLevels = np.array([])
            self.VibrationalLevels = np.array([])

    def add_or_update_all_possible_datasets(self, level_type=None):
        if level_type is not None:
            print(f"Reading {level_type} data")
            MoleculeHelper.add_or_update_dataset(self, level_type)
        else:
            for level_type in ['ElectronicLevels', 'RotationalLevels', 'VibrationalLevels']:
                if level_type in list(vars(self).keys()):
                    print(f"Reading {level_type} data")
                    MoleculeHelper.add_or_update_dataset(self, level_type)

    def __str__(self):
        return f"DiatomicMolecule(name={self.name}, attributes={self.attributes})"

class PolyatomicMolecule:
    def __init__(self, species, read_species=True, ChargeIC=0.0, InteractionID=2.0):

        # read_species: if True, read species from database, if False, create new species and write values later
        # create dictionary to store attributes with keys corresponding to the right name in the database
        self.attributes = {}
        # storing the species name and electronic levels separately (not into the attributes dict) since it should not be written to the database
        self.name = species
        self.new_references = {}
        # read species from database into class variables
        if read_species:
            attrs = hdf_unified_data["Species"][species].attrs
            # Store attributes in class variables
            for attr_name, attr_value in attrs.items():
                self.attributes[attr_name] = attr_value
        # create new class to store new data
        else:
            # initialize the attributes of the class with the configuration from the config file
            for attr_name, attr_value in polyatomic_attributes.items():
                self.attributes[attr_name] = attr_value
            # gets updated in function
            self.attributes['NumOfAtoms'] = MoleculeHelper.get_num_of_atoms(self.name)
            if 'ion' in self.name.lower():
                Ion_Number = int(re.search(r'\d+$', self.attributes['name']).group())
                self.attributes['ChargeIC'] = -Ion_Number * electron_charge
                self.attributes['InteractionID'] = 20.0
            else:
                self.attributes['ChargeIC'] = ChargeIC
                self.attributes['InteractionID'] = InteractionID
            # Get LinearMolec to set right number of charac temps
            while True:
                self.attributes['LinearMolec'] = int(input("Please enter '1' if species is a linear molecule and '0' if not: "))
                if self.attributes['LinearMolec'] in [0, 1]:
                    break
                else:
                    print("Please enter 1 for linear or 0 for non-linear.")
            # set correct number of characteristic temperatures
            if self.attributes['LinearMolec'] == 1:
                for i in range(1, 3*self.attributes['NumOfAtoms']-5):
                    self.attributes[f'CharaTempVib{i}'] = None
                self.attributes['CharaTempRot'] = None
                self.attributes['MomentOfInertia'] = None
            else:
                for i in range(1, 3*self.attributes['NumOfAtoms']-6):
                    self.attributes[f'CharaTempVib{i}'] = None
                for i in range(1, 3):
                    self.attributes[f'CharaTempRot{i}'] = None
                self.attributes['MomentOfInertia'] = np.array([])
            # set inner energy levels
            self.ElectronicLevels = np.array([])
            self.RotationalLevels = np.array([])
            self.VibrationalLevels = np.array([])

    def add_or_update_all_possible_datasets(self, level_type=None):
        if level_type is not None:
            print(f"Reading {level_type} data")
            MoleculeHelper.add_or_update_dataset(self, level_type)
        else:
            for level_type in ['ElectronicLevels', 'RotationalLevels', 'VibrationalLevels']:
                if level_type in list(vars(self).keys()):
                    print(f"Reading {level_type} data")
                    MoleculeHelper.add_or_update_dataset(self, level_type)

    def __str__(self):
        return f"PolyatomicMolecule(name={self.name}, attributes={self.attributes})"

######################################################################################################################################################################################################
# Functions to handle the classes, e.g. create a new class depending on the species, edit attributes or write an instance to the database
######################################################################################################################################################################################################

def get_interaction_id(species_name):
    if sum((1 for c in species_name.replace('Ion', '') if c.isupper())) == 1 and (not bool(re.search('\\d+', re.sub('Ion\\d+', '', species_name)))):
        if not bool(re.search('[A-Za-z]*\\d', re.sub('Ion\\d+', '', species_name))):
            if 'Ion' in species_name:
                interactionID = 10
            elif 'Ion' not in species_name:
                interactionID = 1
    elif bool(re.search('\\d+', re.sub('Ion\\d+', '', species_name))) or sum((1 for c in species_name.replace('Ion', '') if c.isupper())) != 1:
        if 'Ion' in species_name:
            interactionID = 20
        elif 'Ion' not in species_name:
            interactionID = 2
    elif 'electron' in species_name:
        interactionID = 4
    return interactionID

def create_empty_instance(species):
    """
    Function to create an empty instance of the class based on the interaction ID of the species

    The class will have the right attributes corresponding to the interaction ID of the species.
    """
    interaction_id = get_interaction_id(species)
    if interaction_id == 1 or interaction_id == 10:
        cls = Atom
    elif interaction_id == 2 or interaction_id == 20:
        # get number of atoms in the molecule
        num_of_atoms = MoleculeHelper.get_num_of_atoms(species)
        if num_of_atoms == 2:
            cls = DiatomicMolecule
        else:
            cls = PolyatomicMolecule
    else:
        raise ValueError("Unknown interaction ID for the given species.")

    # Create an empty instance of the class and add some default attributes
    instance = cls(species, read_species=False)
    instance.attributes['* Created'] = date.today().strftime("%B %d, %Y")
    instance.attributes['* Reference'] = 'This species was created manually.'
    return instance

def determine_class_with_given_data(attrs):
    """Determines the class of the species based on the attributes of the species"""
    if 'PolyatomicMol' in attrs:
        return PolyatomicMolecule
    else:
        if 'CharaTempRot' in attrs or 'CharaTempVib' in attrs or 'CharaTempRot1' in attrs or 'CharaTempVib1' in attrs:
            return DiatomicMolecule
        else:
            return Atom

def create_instance_from_data(species):
    """Creates an instance of the correct class based on the attributes of the species when reading in species data from the database"""
    attrs = hdf_unified_data["Species"][species].attrs
    cls = determine_class_with_given_data(attrs.keys())
    if cls is None:
        raise ValueError("Unknown class type for the given data.")

    # Create an instance of the determined class
    instance = cls(species)
    return instance

def read_datasets_from_existing_species(instance):
    """Read inner energy data from the existing species database and store it in the class variables"""
    if isinstance(instance, Atom):
        level_types = ['Electronic']
        class_variables = ['ElectronicLevels']
    else:
        level_types = ['Electronic', 'Rotational', 'Vibrational']
        class_variables = ['ElectronicLevels', 'RotationalLevels', 'VibrationalLevels']

    for level_type, class_variable in zip(level_types, class_variables):
        try:
            setattr(instance, class_variable, np.array(hdf_unified_data["Species"][instance.name][f'{level_type}Level']))
            print(f"Reading {level_type} data from the database for species ", green(instance.name))
        except Exception as e:
            print(red('Error'), f"reading {level_type} data from the database for ", green(instance.name), ":", e)
            user_input = get_valid_input(create_prompt('to add new data from csv file',
                                                       'to skip data for this species'),
                                                       lambda x: x == '1' or x == '2' or x == '3')
            if user_input == '1':
                MoleculeHelper.add_or_update_dataset(instance, level_type, create_new=True)
            elif user_input == '2':
                print(f"Skipping {class_variable} for species ", green(instance.name), "\n")
            elif user_input == '3':
                own_exit()

def check_attributes_from_actc(instance, atct_dict):
    """Function to check the attributes of the species with data from the ATcT database

    This only happens for 4 attributes: HeatOfFormation_K, MassIC, ChargeIC, and InteractionID.
    Parameters:
    instance: Instance of the species class
    atct_dict: Dictionary containing the data from the ATcT database, e.g.
               atct_data_for_instance={
                'HeatOfFormation_K' : example_value,
                'MassIC'            : example_value,
                'ChargeIC'          : example_value,
                'InteractionID'     : example_value,
                }
                with this dictionary the attributes to check are defined! So if more data is retrieved from ATcT, the dictionary has to be updated, but the function just loops over the dictionary!
    """
    print(f"\nChecking attributes: {list(atct_dict.keys())} from {ATcT_URL} for species ", green(instance.name))
    # create list to see if user input is necessary
    diffs_list = []
    for attr_name, attr_value in atct_dict.items():
        try:
            # heat of formation is not set for ions in the database, use ions as sanity check together with electronic data
            if attr_name == 'HeatOfFormation_K' and 'HeatOfFormation_K' not in instance.attributes.keys():
                Ion_Number =  int(re.search(r'\d+$', instance.name).group())
                ground_state = re.sub(r'Ion\d+','',instance.name)
                HeatOfFormation_Sum = hdf_unified_data['Species'][ground_state].attrs['HeatOfFormation_K']+hdf_unified_data['Species'][ground_state]['ElectronicLevel'][:][-1,1]
                for i in range(1, Ion_Number):
                    prev_state = ground_state + 'Ion' + str(i)
                    HeatOfFormation_Sum = HeatOfFormation_Sum + hdf_unified_data['Species'][prev_state]['ElectronicLevel'][:][-1,1]
                if not np.isclose(atct_dict['HeatOfFormation_K'], HeatOfFormation_Sum, rtol=1e-03, atol=0.0):          # HeatOfFormation_K being value from ATcT and HeatOfFormation_Sum cal value with relative tolerance of 1e-05
                    print(red('Caution:')+' Please check the heat of formation of the current species: '+green(instance.name)+', since it did not matcht with the value from '+blue(ATcT_URL)+'!\nIt is calculated in the PICLas source code like this: HeatOfFormation(CurrentSpecies) = HeatOfFormation(previous state) + Ionization Energy(last entry in electronic level dataset of previous state)')
                    print(yellow("Calculated HeatOfFormation: ")+str(HeatOfFormation_Sum)+"\n"+blue("HeatOfFormation form ATcT: ")+str(atct_dict['HeatOfFormation_K']))
            else:   # compare the attributes specified by the directory
                if np.isclose(instance.attributes[attr_name], attr_value, rtol=1e-03, atol=1e-08):
                    print('%s for ' % (attr_name) +green(instance.name)+' is equal so will be kept')
                else:
                    print(underlinE(attr_name)+' of species '+green(instance.name)+'\n'+yellow('unified species database')+f': {instance.attributes[attr_name]}'+"\n"+blue(ATcT_URL) +f': {attr_value}\n')
                    diffs_list.append(attr_name)
        except KeyError:
            if attr_name in ['HeatOfFormation_K','MassIC']:
                print("It seems like %s is not set for %s so it will be set as %s with data from "% (attr_name,instance.name,attr_value)+blue(ATcT_URL))
            elif attr_name in ['ChargeIC','InteractionID']:
                print("It seems like %s is not set for %s so it will be set as %s"% (attr_name,instance.name,attr_value))

    if diffs_list == []:    # all attributes are equal
        pass
    else:
        user_input = get_valid_input(create_prompt('to keep attributes from '+yellow('unified species database'),
                                                   'to save attributes from '+blue(ATcT_URL),
                                                   'to exit program here'),
                                                   lambda x: x == '1' or x == '2' or x == '3')
        if user_input == '1':
            print("Keeping attributes for species ",green(instance.name), "\n")
        elif user_input == '2':
            for attr_name, attr_value in atct_dict.items():
                instance.attributes[attr_name] = attr_value
            instance.new_references[attr_name] = str([atct_dict.keys()] + 'from ATcT Version: ' + version + ',' + date.today().strftime("%B %d, %Y"))
        elif user_input == '3':
            own_exit()

def edit_attributes(instance):
    """
    Function to edit attributes of existing species in class

    Parameters:
    instance: Instance of the species class
    """
    # list existing attributes
    print("\nExisting attributes for species ", green(instance.name))
    for attr_name, attr_value in instance.attributes.items():
        if isinstance(attr_value, np.ndarray):
            # used for reference data
            if attr_value.ndim == 2:
                print(f"{attr_name}:")
                for i in range(attr_value.shape[0]):
                    print(f"\t{attr_value[i,0]} : {attr_value[i,1]}")
            # used for array data, e.g. CharaTempVib or CharaTempRot
            elif attr_value.ndim == 1:
                print(f"{attr_name}:")
                for i in range(attr_value.shape[0]):
                    print(f"\tMode {i+1} : {attr_value[i]}")
        else:
            print(f"{attr_name}: {attr_value}")
    user_input_attrs = input(bold('\nPlease enter attribute list as comma separated string. It is also possible to add new attributes.') + ', e.g. Tref, HeatOfFormation, ... or None to skip this species\n')
    attrs_list = user_input_attrs.split(',')
    if user_input_attrs.lower() == 'none':
        return
    # check if attribute exists, if yes user input to overwrite or not
    for attr_name in attrs_list:
        attr_name = attr_name.replace(" ", "")
        if attr_name in list(instance.attributes.keys()):
            print("\n\nThe attribute " + attr_name + " already exists. Do you want to overwrite the value?")
            user_input = get_valid_input(create_prompt('yes',
                                                       'no'),
                                                       lambda x: x == '1' or x == '2' or x =='3')
            if user_input == "1":
                new_data = input(bold('\nPlease enter ' + attr_name + ' of current species (arrays can be entered as comma separated strings, e.g. 1.2,3.4,5.9) %s\n-->') % instance.name)
                if ',' in new_data:
                    new_data = [np.float64(val) for val in new_data.split(',')]
                    # sanity check if array has same length as before
                    if len(new_data) != len(instance.attributes[attr_name]):
                        print(red('Error')+': New array has not the same length as the old one! Please check the input!')
                        own_exit()
                # Replace the existing value
                instance.attributes[attr_name] = attribute_types[attr_name](new_data)
            elif user_input == "3":
                own_exit()
        else:
            if attr_name in attribute_types.keys():
                new_data = input(bold('\nPlease enter ' + attr_name + ' of current species (arrays can be entered as comma separated strings, e.g. 1.2,3.4,5.9) %s\n-->') % instance.name)
                if ',' in new_data:
                    new_data = [np.float64(val) for val in new_data.split(',')]
            else:
                # User input for attribute type
                attr_type = input(f"Please enter the datatype for {attr_name} (str, int, float or array): ").strip().lower()
                if attr_type == 'str':
                    attribute_types[attr_name] = str
                    attr_type_str = 'str'
                elif attr_type == 'int':
                    attribute_types[attr_name] = int
                    attr_type_str = 'int'
                elif attr_type == 'float':
                    attribute_types[attr_name] = float
                    attr_type_str = 'float'
                elif attr_type == 'array':
                    attribute_types[attr_name] = np.array
                    attr_type_str = 'np.array'
                else:
                    print("Invalid datatype entered!")
                    own_exit()
                new_data = (input(bold('\nPlease enter ' + attr_name + ' of current species (arrays can be entered as comma separated strings, e.g. 1.2,3.4,5.9) %s\n-->') % instance.name))
                if ',' in new_data:
                    new_data = [np.float64(val) for val in new_data.split(',')]

                # Add the new attribute type to the config file
                config_path = os.path.join(os.path.dirname(__file__), 'config.py')
                # decide which class to append to
                if instance.attributes['InteractionID'] == 1 or instance.attributes['InteractionID'] == 10:
                    line_catch = 'atom_attributes = {'
                elif instance.attributes['InteractionID'] == 2 or instance.attributes['InteractionID'] == 20:
                    if 'PolyatomicMol' in instance.attributes.keys():
                        line_catch = 'polyatomic_attributes = {'
                    else:
                        line_catch = 'diatomic_attributes = {'
                # Write the new attribute type to config.py
                with open(config_path, 'r+') as config_file:
                    lines = config_file.readlines()
                    config_file.seek(0)
                    for line in lines:
                        config_file.write(line)
                        if line.strip() == 'attribute_types = {':
                            config_file.write(f"    '{attr_name}'         : {attr_type_str},\n")
                        if line.strip() == line_catch:
                            config_file.write(f"    '{attr_name}'         : None,\n")
                    config_file.truncate()

            # sanity check in both cases for characteristic temperatures
            if 'LinearMolec' in instance.attributes.keys():
                if attr_name == 'CharaTempRot' or attr_name == 'MomentOfInertia':
                    print(f'Current species {instance.name} is a linear molecule, so only one characteristic temperature/moment of inertia is needed!')
                    own_exit()
                if attr_name == 'CharaTempVib':
                    if len(new_data) != 3*instance.attributes['NumOfAtoms']-5:
                        print(red('Error')+f': New array for characteristic vibrational temperatures does not fit expected value for a linear molecule with {instance.attributes["NumOfAtoms"]} atoms!')
                        own_exit()
            else:
                if attr_name == 'CharaTempRot' or attr_name == 'MomentOfInertia':
                    if len(new_data) != 3:
                        print(red('Error')+': New array for characteristic rotational temperatures/moment of inertia does not fit expected value for a non-linear molecule!')
                        own_exit()
                    if attr_name == 'CharaTempVib':
                        if len(new_data) != 3*instance.attributes['NumOfAtoms']-6:
                            print(red('Error')+f': New array for characteristic vibrational temperatures does not fit expected value for a non-linear molecule with {instance.attributes["NumOfAtoms"]} atoms!')
                            own_exit()
            # if passed write to class
            instance.attributes[attr_name] = attribute_types[attr_name](new_data)

def write_instance_to_database(instance):
    """
    This function writes the attributes and inner energy levels/degeneracies of an instance to the database

    Parameters:
    instance: Instance of the species class
    """
    # write attributes to the database
    deleted = False
    if instance.name in list(hdf_unified_data["Species"].keys()):
        # if species exists in database, save the inner energy levels and degeneracies to temporary variables
        try:
            if not isinstance(instance, Atom):
                temp_rot = np.array(hdf_unified_data["Species"][instance.name]['RotationalLevel'])
                temp_vib = np.array(hdf_unified_data["Species"][instance.name]['VibrationalLevel'])
        except Exception as e:
            # catch for missing datasets
            temp_rot = None
            temp_vib = None
            pass
        temp_elec = np.array(hdf_unified_data["Species"][instance.name]['ElectronicLevel'])
        # delete species from database
        del hdf_unified_data["Species"][instance.name]
        deleted = True
    else: # default values for temporary variables to ensure correct functionality in zip loop
        temp_rot, temp_vib = None, None
    hdf_unified_data["Species"].create_group(instance.name)
    for attr_name, attr_value in instance.attributes.items():
        # check if attribute of class is not set yet and give user the option to set it
        if attr_value is None:
            print(f"Attribute {attr_name} is not set for species {instance.name}!")
            user_input = get_valid_input(create_prompt('to set attribute',
                                                       'to skip attribute'),
                                                       lambda x: x == '1' or x == '2' or x == '3')
            if user_input == '1':
                instance.attributes[attr_name] = attribute_types[attr_name](input(bold('\nPlease enter ' + attr_name + ' of current species %s\n-->') % instance.name))
            elif user_input == '2':
                continue
        hdf_unified_data["Species"][instance.name].attrs[attr_name] = attr_value

    # write inner energy levels
    if not deleted:
        # species was not in database before, is added
        hdf_unified_data["Species"][instance.name].create_dataset('ElectronicLevel', data=instance.ElectronicLevels)
    else:
        if 'ElectronicLevel' in list(vars(instance).keys()):
            # check if electronic level was set in class, if yes use this data
            hdf_unified_data["Species"][instance.name].create_dataset('ElectronicLevel', data=instance.ElectronicLevels)
        else:
            # if not, restore data with temporary variable
            hdf_unified_data["Species"][instance.name].create_dataset('ElectronicLevel', data=temp_elec)
    # equivalent for rotational and vibrational levels
    if not isinstance(instance, Atom):
        for level_type, temp_level in zip(['Rotational', 'Vibrational'], [temp_rot, temp_vib]):
            if not deleted: # new species
                hdf_unified_data["Species"][instance.name].create_dataset(f'{level_type}Level', data=instance.RotationalLevels if level_type == 'Rotational' else instance.VibrationalLevels)
            else:
                level_type_matching = level_type + 'Levels'
                if level_type_matching in list(vars(instance).keys()): # new data for existing species
                    hdf_unified_data["Species"][instance.name].create_dataset(f'{level_type}Level', data=instance.RotationalLevels if level_type == 'Rotational' else instance.VibrationalLevels)
                else: # restore
                    if temp_level is not None:
                        hdf_unified_data["Species"][instance.name].create_dataset(f'{level_type}Level', data=temp_level)
                    else: # dataset did not exist before
                        pass

    # update references
    num_of_new_references = len(instance.new_references)
    for attr_name, attr_value in instance.new_references.items():
        if attr_name == '* Reference':
            if np.isscalar(attr_value):
                if isinstance(attr_value, str):
                    strings = attr_value.split(':')
                else:
                    strings = attr_value.astype(str).split(':')
                new_refs = np.empty((len(strings)//2+num_of_new_references,2), dtype=datatype_h5)
                for j in range(len(strings)//2):
                    new_refs[j,0] = strings[2*j]
                    new_refs[j,1] = strings[2*j+1]
            elif isinstance(attr_value, np.ndarray):
                attr_value = attr_value.astype(str)
                new_refs = np.empty((attr_value.shape[0]+num_of_new_references,2), dtype=datatype_h5)
                for j in range(attr_value.shape[0]):
                    new_refs[j,0] = attr_value[j,0]
                    new_refs[j,1] = attr_value[j,1]

            for j, attr_name, attr_value in enumerate(instance.new_references.items()):
                new_refs[len(strings)//2+j,0] = attr_name
                new_refs[len(strings)//2+j,1] = attr_value

            del hdf_unified_data["Species"][instance.name].attrs['* Reference']
            hdf_unified_data["Species"][instance.name].attrs.create('* Reference', new_refs)


######################################################################################################################################################################################################
#   functions for getting electronic levels data from NIST database
######################################################################################################################################################################################################
def int_to_Roman(num):
    val = (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
    syb = ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
    roman_num = ''
    for i in range(len(val)):
        count = int(num / val[i])
        roman_num += syb[i] * count
        num -= val[i] * count
    return roman_num

def get_data_from_NIST(CURRENT_SPECIES, ION_LEVEL):
    """Function to get the electronic data from the NIST database"""
    # Build the species in the NIST database query format: Xe+I
    species = CURRENT_SPECIES.split('Ion')[0]
    current_species_NIST = species + '+' + int_to_Roman(ION_LEVEL)
    print(bold('Converted '+CURRENT_SPECIES+' to '+current_species_NIST+' for download from NIST database.'))
    # Build the request URL
    URL_request = URL_base + '?de=0&spectrum=' + current_species_NIST + '&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&term_out=on&level_out=on&j_out=on&temp=&submit=Retrieve+Data'
    response = requests.get(URL_request)
    # Check the website
    if response.status_code == 200:
        print('Download from NIST successful.\n')
    # Get the data from the NIST website
    data = response.content.decode('utf-8')
    # Clean the data
    data = data.replace('[','')
    data = data.replace(']','')
    data = data.replace('(','')
    data = data.replace(')','')
    data = data.replace('= ','')
    data = data.replace('=','')
    data = data.replace(' ','')
    data = data.replace('"','')
    try:
        # Attempt to convert data to pandas dataframe
        data = pd.read_csv(io.StringIO(data),skipinitialspace=True,delimiter=",",usecols=['J','Prefix','Levelcm-1','Suffix','Term'], na_values=['---'])
        # Drop 'Prefix' and 'Suffix' columns
        data = data.drop(columns=['Prefix', 'Suffix'])
        return data, current_species_NIST
    except Exception as e:
        # check if species is not found because it is fully ionized in database
        if "FullyIonized" in hdf_unified_data["Species"][CURRENT_SPECIES].attrs:
            return int(-1), current_species_NIST
        # Print an error message if there's an exception
        print(red('Error'), "converting data from \n", blue(URL_request), " for ", green(CURRENT_SPECIES), "might not be available on requested URL. Tried to access data with ", green(current_species_NIST), "Please check this combination on the URL. \nError:", e ,"\n")
        # get user input to determine which dataset should be saved
        user_input_create = get_valid_input(create_prompt('to create '+green(CURRENT_SPECIES)+' without the electronic level',
                                                       'to exit program here'),
                                                       lambda x: x == '1' or x == '2')
        if user_input_create == '1':
            column_names = ['J', 'Levelcm-1', 'Term']
            data = pd.DataFrame(columns=column_names)
            first_row_data = {'J': '0', 'Levelcm-1': '0', 'Term': 'Limit'}
            data = data.append(first_row_data, ignore_index=True)
        elif user_input_create == '2':
            return int(-1), current_species_NIST

# converting electronic data in useful format
def convert_electronic_data(DATA):
    """Converts the electronic data from the NIST database to a useful format"""
    # get index of J column
    J     = DATA.columns.get_loc("J")
    # find the first limit in term column to drop every entry after ionization
    index_of_limit_ionization = DATA[DATA['Term'].str.contains('Limit', na=False)].index[0]
    # drop every entry after limit index
    DATA = DATA.drop(DATA.index[index_of_limit_ionization+1:])
    # set J of limit row to 0 so its not lost in nan computation and it later becomes g=1
    DATA.iloc[index_of_limit_ionization,J] = 0
    max_level = len(DATA)-1
    # dropping term column for nan filtering
    DATA = DATA.drop(columns=['Term'])
    # find rows containing nan
    rows_with_nan = [index for index, row in DATA.iterrows() if row.isnull().any()]
    # remove rows with nan
    if len(rows_with_nan) == 1:
        max_level = rows_with_nan[0]
        drop_to_end = 1 - (len(DATA)-max_level)
        if drop_to_end < 0:
            DATA = DATA.iloc[:drop_to_end]
        elif drop_to_end == 0:
            pass
        else:
            print("ERROR: drop_to_end must be negative!")
            exit(1)
    else:
        for val in rows_with_nan:
            DATA = DATA.drop(val)
        max_level = len(DATA)-1
    # get J again because 'Term' was deleted
    J     = DATA.columns.get_loc("J")
    # set last J to 0 so ionization is later (1,UNIFIED_DATA_ARRAY)
    if DATA['J'].dtype == 'float64':
        DATA.iloc[max_level,J] = 0.0
    else:
        DATA.iloc[max_level,J] = "0.0"
    # Check the DATAtype: if its a float, then all non-numerical characters have already been removed
    if DATA['Levelcm-1'].dtype != 'float64':
        # Drop rows with a question mark ("This level/line may not be real.")
        DATA.drop(DATA[DATA['Levelcm-1'].str.contains(r'[?]')].index,inplace=True)
        # Drop rows with a +x ("The relative positions of the levels within such a system are accurate within experimental uncertainties, but no experimental connection between this system and the other levels of the spectrum has been made.")
        DATA.drop(DATA[DATA['Levelcm-1'].str.contains(r'[+x]')].index,inplace=True)
    # Execute fractions and convert J to g
    if DATA['J'].dtype != 'float64':
        for i in range(len(DATA['J'])):
            # Convert fraction entries to floats (5/2 -> 2.5)
            #print(type(DATA.iloc[i,J]))
            found = re.search(r'\d+/\d+', DATA.iloc[i,J])
            if found:
                numbers = DATA.iloc[i,J].split("/")
                DATA.iloc[i,J] = float(numbers[0])/float(numbers[1])
            else:
                DATA.iloc[i,J] = float(DATA.iloc[i,J])
    # Convert J to g
    DATA.iloc[:,J] = 2 * DATA.iloc[:,J] + 1.0
    # Convert type for HDF5 output
    DATA['J'] = DATA['J'].astype(float)
    # make sure first level is (1,0)
    DATA.iloc[0,J] = 1
    # Convert 1/cm to K
    DATA['Levelcm-1'] = DATA['Levelcm-1'].astype(float)
    DATA['Levelcm-1'] = 100 * DATA['Levelcm-1'] * 1.986E-025 / 1.38065E-023             # 1/cm * 100cm/m * (J m) / (J/K) = K
    # check for energy limits if they are increasing with every level
    x_old=0.
    Level = DATA.columns.get_loc("Levelcm-1")
    for i in range(len(DATA['Levelcm-1'])):
            x = DATA.iloc[i,Level]
            if x < x_old:
                print('Error in level %s: the energy is not increasing with the levels E2=%s < E1=%s' % (i,x,x_old))
                exit(1)
            else:
                x_old = x
    # Write to hdf: If DATA set already exists, delete the old set first
    DATA.columns=['Degeneracy', 'Level(K)']
    return DATA

# function to print differences in datasets
def print_diffs(new_data, ref_data, level_type='Electronic'):
    '''Print the differneces in the two data sets onto terminal\n Inputs: new data array from url/custom data and data array from species database'''
    try:
        diff_indices_degeneracy = np.where(~np.isclose(new_data[:,0], ref_data[:,0], rtol=1e-05, atol=0))
        diff_indices_levels = np.where(~np.isclose(new_data[:,1], ref_data[:,1], rtol=1e-05, atol=0))
        # special case for vibrational levels since no degeneracy is given
        if level_type == 'Vibrational':
            degeneracy = np.ones(len(new_data[:,0]))
            new_data = np.column_stack((degeneracy,new_data))
            ref_data = np.column_stack((degeneracy,ref_data))
        # Create DataFrame for levels
        diff_df_levels = pd.DataFrame({
            'Index Level': diff_indices_levels[0],
            'URL Array Level': [new_data[idx, 1] for idx in diff_indices_levels[0]],
            'Unified Data Array Level': [ref_data[idx, 1] for idx in diff_indices_levels[0]],
        })
        # Create DataFrame for degeneracy
        diff_df_degeneracy = pd.DataFrame({
            'Index Degeneracy': diff_indices_degeneracy[0],
            'URL Array Degeneracy': [new_data[idx, 0] for idx in diff_indices_degeneracy[0]],
            'Unified Data Array Degeneracy': [ref_data[idx, 0] for idx in diff_indices_degeneracy[0]],
        })
        # Set display option to show all rows
        with pd.option_context('display.max_rows', None):
            # Print DataFrames separately with a clear separator
            print("\n" + underlinE('Differnce for Levels:'))
            print(tabulate(diff_df_levels, headers='keys', tablefmt='psql', numalign='center'))
            print("\n" + underlinE('Differnce for Degeneracy:'))
            print(tabulate(diff_df_degeneracy, headers='keys', tablefmt='psql', numalign='center'))
    except Exception as e:
        if 'operands could not be broadcast together' in str(e):
            # print Error message with colors if expected errors occurs
            print(bold('\033[91mError:\033[0m'),"operands could not be broadcast together with shapes \033[33m%s\033[0m and \033[34m%s\033[0m!" % (ref_data.shape,new_data.shape)+"\n")
        else:
            # print error message if not expected error occurs
            print(e)
    return

######################################################################################################################################################################################################
#   functions for getting data from the ATcT
######################################################################################################################################################################################################
def create_species_dict():
    """
    Function to create a dictionary with species data from the ATcT database

    Returns:
    A dictionary with species names as keys and a tuple with the following values as values:
    - Heat of formation at 298.15K in J/mol
    - Mass in kg
    - Charge of the species in C
    - Interaction ID of the species
    """
    print('Requesting thermodynamic data from ATcT database...')
    species_dict = {}
    # Create an HTML session
    session = HTMLSession()
    # Make a GET request to the URL
    response = session.get(ATcT_URL)
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        print('Request at ATcT successful.\n')
    # Render the HTML content (execute JavaScript)
    response.html.render()
    # Parse the HTML content
    soup = BeautifulSoup(response.html.html, 'html.parser')
    # Navigating to wanted data
    grand_table = soup.find(class_ = 'bodyProperInner')
    tbody = grand_table.find('tbody')
    tr = tbody.find('tr')
    td = tr.find('td')
    inner_table = td.find('table')
    inner_tbody = inner_table.find_all('tbody')
    # data not easily accesable due to missing ids and many tables,bodies,etc -> 5 must equal the number of the tbody containing the species names, heat, ...
    inner_tbody = inner_tbody[5]
    rows = inner_tbody.find_all('tr')
    # slice first entry since it only shows the structure
    rows = rows[1:]
    for row in rows:
        # get species name
        row_str = str(row)
        start_index = row_str.index('type="button">') + len('type="button">')
        end_index = row_str.index(' </button>')
        formula = row_str[start_index:end_index]
        formula = formula.replace('[','')
        formula = formula.replace(']','')
        formula = formula.replace(' ','')
        if formula.endswith("(g)") and '(' not in formula[:-3]: # only save data of species in gas phase
            formula = re.sub(r'\([a-zA-Z0-9,\/]+\)','',formula)
        else:
            continue
        # get heat of formation at 298.15K
        bkgDHf298_element = row.find('td', class_='bkgDHf298')
        if bkgDHf298_element:
            # Extract the text content from the <span> element inside <td> and convert to kelvin
            delta_f_H_298K = float(bkgDHf298_element.span.text) * 120.27311595710411
        # get molecular mass of species
        bkgMass_element = row.find('td', class_='bkgMass')
        if bkgMass_element:
            # Extract the text content from the <span> element inside <td>
            mol_mass = bkgMass_element.span.text
            # clean data
            index_of_plus = mol_mass.find('±')
            if index_of_plus != -1:
                mol_mass = float(mol_mass[:index_of_plus-1])
                # convert with avogadro const and to kg
                mass = mol_mass/(6.02214076e23 * 1000)
        # convert species name to fit database
        if '+' in formula:
            Ion_Number = re.search(r'\d+$', formula)
            if Ion_Number is None:
                species_name = re.sub(r'\+','Ion1',formula)
            else:
                Ion_Number = int(re.search(r'\d+$', formula).group())
                species_name = re.sub(r'\+\d+$',f'Ion{Ion_Number}',formula)
        else:
            species_name = formula
        # get charge of current species
        if 'Ion' in species_name:
            Ion_Number =  int(re.search(r'\d+$', species_name).group())
            species_charge = - Ion_Number * electron_charge
        elif 'electron' in species_name:
            species_charge = electron_charge
        else:
            species_charge = 0
        # get interactionID from current species name
        interactionID = get_interaction_id(species_name)
        species_dict[species_name] = (float(f"{delta_f_H_298K:.8e}"), float(f"{mass:.8e}"), species_charge, interactionID)
    # add electron to species_dict since it is not in ATcT database and would cause a keyerror later
    species_dict['electron'] = (0.0, 9.11e-31, electron_charge, 4)
    return species_dict