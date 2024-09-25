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

###################################################################################################
#   functions for electronic levels Data handeling 
###################################################################################################
def get_data_from_NIST(CURRENT_SPECIES, SPECIES, ION_LEVEL, URL_BASE, HDF_UNIFIED_DATA):
    # Build the species in the NIST database query format: Xe+I
    current_species_NIST = SPECIES + '+' + int_to_Roman(ION_LEVEL)
    print(bold('Converted '+CURRENT_SPECIES+' to '+current_species_NIST+' for download from NIST database.'))
    # Build the request URL
    URL_request = URL_BASE + '?de=0&spectrum=' + current_species_NIST + '&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&term_out=on&level_out=on&j_out=on&temp=&submit=Retrieve+Data'
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
        if "FullyIonized" in HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES].attrs:
            return int(-1), current_species_NIST
        # Print an error message if there's an exception
        print(red('Error'), "converting data from \n", blue(URL_request), " for ", green(CURRENT_SPECIES), "might not be available on requested URL. Tried to access data with ", green(current_species_NIST), "Please check this combination on the URL. \nError:", e ,"\n")
        # get user input to determine which dataset should be saved
        user_input_create = get_valid_input("Do you want to create "+green(CURRENT_SPECIES)+" without the electronic level?\n" + bold('Please enter\n ') + purple('1') + " for yes or \n " + purple('2') + " for no\n-->", lambda x: x == '1' or x == '2')
        if user_input_create == '1':
            column_names = ['J', 'Levelcm-1', 'Term']
            data = pd.DataFrame(columns=column_names)
            first_row_data = {'J': '0', 'Levelcm-1': '0', 'Term': 'Limit'}
            data = data.append(first_row_data, ignore_index=True)
        elif user_input_create == '2':
            return int(-1), current_species_NIST

# converting electronic data in useful format
def convert_electronic_data(DATA):
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

# create dataset and save attributes
def create_dataset(HDF_UNIFIED_DATA, CURRENT_SPECIES, DATA_ARRAY):
    '''Create new dataset while saving old attributes and adding a new reference\n Inputs: database where dataset is saved and species name'''
    # Store existing attributes
    attrs = HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES].attrs
    del HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES]['ElectronicLevel']
    dataset = HDF_UNIFIED_DATA['Species'][CURRENT_SPECIES].create_dataset('ElectronicLevel', data=DATA_ARRAY)
    # Restore attributes
    for attr_name, attr_value in attrs.items():
        dataset.attrs[attr_name] = attr_value
    # Get current references -> will be an array?!
    current_references = [attrs['* Reference']]
    levels_reference = 'Levels and Degeneracy - Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2020). NIST Atomic Spectra Database (ver. 5.8), [Online]. Available: https://physics.nist.gov/asd. National Institute of Standards and Technology, Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F. Retrieved on ' + date.today().strftime("%B %d, %Y") + '.'
    # Append the new reference
    current_references.append(levels_reference)
    # Update the attribute with the references
    del dataset.attrs['* Reference']
    dataset.attrs.create('* Reference', current_references)

# function to print differences in datasets
def print_diffs(DATA_ARRAY, UNIFIED_DATA_ARRAY, CURRENT_SPECIES_NIST):
    '''Print the differneces in the two data sets onto terminal\n Inputs: new data array from url and data array from species database'''
    try:
        diff_indices_degeneracy = np.where(~np.isclose(DATA_ARRAY[:,0], UNIFIED_DATA_ARRAY[:,0], rtol=1e-05, atol=0))
        diff_indices_levels = np.where(~np.isclose(DATA_ARRAY[:,1], UNIFIED_DATA_ARRAY[:,1], rtol=1e-05, atol=0))
        # Create DataFrame for levels
        diff_df_levels = pd.DataFrame({
            'Index Level': diff_indices_levels[0],
            'URL Array Level': [DATA_ARRAY[idx, 1] for idx in diff_indices_levels[0]],
            'Unified Data Array Level': [UNIFIED_DATA_ARRAY[idx, 1] for idx in diff_indices_levels[0]],
        })
        # Create DataFrame for degeneracy
        diff_df_degeneracy = pd.DataFrame({
            'Index Degeneracy': diff_indices_degeneracy[0],
            'URL Array Degeneracy': [DATA_ARRAY[idx, 0] for idx in diff_indices_degeneracy[0]],
            'Unified Data Array Degeneracy': [UNIFIED_DATA_ARRAY[idx, 0] for idx in diff_indices_degeneracy[0]],
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
            print(bold('\033[91mError:\033[0m'),"operands could not be broadcast together with shapes \033[33m%s\033[0m from \033[33munified species database\033[0m and \033[34m%s\033[0m with data from \033[34mURL\033[0m accessed as \033[32m%s\033[0m" % (UNIFIED_DATA_ARRAY.shape,DATA_ARRAY.shape,CURRENT_SPECIES_NIST)+"\n")
        else:
            # print error message if not expected error occurs
            print(e)
    return

def check_datasets(current_species, HDF_UNIFIED_DATA, RELATIVE_PATH):
    # Base URL of the query
    URL_base = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
    # Get the ionization number as the last digits and 1 to comply with the NIST standard
    ion_level = int(re.sub('.*?([0-9]*)$',r'\1',current_species) or 0) + 1
    # Get the species (everything before Ion)
    species = current_species.split('Ion')[0]
    # break for molecules
    if bool(re.search(r'[A-Za-z]*\d', species)) or sum(1 for c in species.replace('Ion','') if c.isupper()) != 1:
        try:
            # Attempt to convert data to pandas dataframe from custom csv file
            data = pd.read_csv(f'custom_electronic_levels_{current_species}.csv',skipinitialspace=True,delimiter=",",usecols=['J','Levelcm-1','Term'], na_values=['---'], dtype=str)
            print(f'Using custom data from custom_electronic_levels_{current_species}.csv')
        except:
            print(red('Caution')+": Please make sure you provide a atom species. Automatic read in from online database of electronic levels for molecules not possible yet! Current species: "+green(current_species)+" It is possible to provide a custom csv file containing the electronic levels. To use this option please use the template in this directory, insert the necessary information and add the species to the end like this: 'custom_electronic_levels_H2.csv'")
            return
    else:
        # get data from Nist database
        data, current_species_NIST = get_data_from_NIST(current_species, species, ion_level, URL_base, HDF_UNIFIED_DATA)
        if type(data) == int:
            return
    
    # convert data from website to fit specifications (e.g. filter unnecessary information )
    data = convert_electronic_data(data)
    
    # check if species is already in Species database
    if current_species in HDF_UNIFIED_DATA['Species'].keys():
        # data from URL into numpy array
        data_array = data.to_numpy()
        # data from species database into numpy array
        unified_data_array = np.array(HDF_UNIFIED_DATA["Species"][current_species]['ElectronicLevel'])

        # check if read-in data from species database has wrong format (might be empty)
        if unified_data_array.ndim == 1 or unified_data_array.shape[1] != 2:
            print(red('Error'), ": Read-in from ", cyan(os.path.abspath(RELATIVE_PATH))," failed or dataset has unexpected dimensions (needs to be (X,2)! Please check the dataset of species ", green(current_species), " - species will be skipped\n")
            # print("What should be done in this case? Dimension do not match")
            # old_attrs = HDF_UNIFIED_DATA["Species"][current_species].attrs
            # del HDF_UNIFIED_DATA["Species"][current_species]
            # dataset = HDF_UNIFIED_DATA['Species'].create_dataset(current_species, data=data_array)
            return
        # Compare the arrays to check if they are already the same
        if np.array_equal(data_array, unified_data_array):
            print("\nDatasets of species "+ green(current_species) + " are equal, so dataset from " + yellow('unified species database') + " will be kept.\n")
            return

        # get user input to determine which dataset should be saved
        print("It seems like the species " + green(current_species) + " is already stored in the SpeciesDatabase at " + cyan(os.path.abspath(os.path.abspath(RELATIVE_PATH)))+' with the following data:\n')
        print_diffs(data_array, unified_data_array,current_species_NIST)

        # Prompt user to choose which dataset to keep
        user_input = get_valid_input(create_prompt('to keep data and attributes from '+yellow('unified species database'), 'to save only electronic level data from '+blue(URL_base), 'to skip all electronic levels and continue with only attributes'), lambda x: x == '1' or x == '2' or x == '3' or x == '4')
        if user_input == '1':
            print("Keeping electronic level dataset and attributes for species ",green(current_species), "\n")
            return
        elif user_input == '2':
            create_dataset(HDF_UNIFIED_DATA, current_species, data_array)
            print("Saving electronic level dataset from " + blue('URL') + " for species ",green(current_species), " but keeping attributes\n")
        elif user_input == '3':
            return -1           # to break loop in maintain_database.py
        elif user_input == '4':
            print(bold(red("Exiting program")))
            exit(1)

    else:
        print("Species "+green(current_species)+" not found in database. Please create new species.")

def add_dataset(HDF_UNIFIED_DATA,current_species,RELATIVE_PATH,ATCT_URL,SPECIES_DICT,SPECIES_DICT_FLAG):
    # Base URL of the query
    URL_base = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
    # Get the ionization number as the last digits and 1 to comply with the NIST standard
    ion_level = int(re.sub('.*?([0-9]*)$',r'\1',current_species) or 0) + 1
    # Get the species (everything before Ion)
    species = current_species.split('Ion')[0]
    # break for molecules
    if bool(re.search(r'[A-Za-z]*\d', species)) or sum(1 for c in species.replace('Ion','') if c.isupper()) != 1:
        try:
            # Attempt to convert data to pandas dataframe from custom csv file
            data = pd.read_csv(f'custom_electronic_levels_{current_species}.csv',skipinitialspace=True,delimiter=",",usecols=['J','Levelcm-1','Term'], na_values=['---'], dtype=str)
            print(f'Using custom data from custom_electronic_levels_{current_species}.csv')
        except:
            print(red('Caution')+": Please make sure you provide a atom species. Automatic read in from online database of electronic levels for molecules not possible yet! Current species: "+green(current_species)+" It is possible to provide a custom csv file containing the electronic levels. To use this option please use the template in this directory, insert the necessary information and add the species to the end like this: 'custom_electronic_levels_H2.csv'")
            return


    else:
        # get data from Nist database
        data, current_species_NIST = get_data_from_NIST(current_species, species, ion_level, URL_base, HDF_UNIFIED_DATA)
        if type(data) == int:
            print("Error")
            return
    
    # convert data from website to fit specifications (e.g. filter unnecessary information )
    data = convert_electronic_data(data)
    # check if exists
    if current_species not in HDF_UNIFIED_DATA['Species'].keys():
        print("Species "+ green(current_species) + " not found in Species Database - new dataset is being generated at "+ cyan(os.path.abspath(RELATIVE_PATH)))
        dataset = HDF_UNIFIED_DATA['Species'][current_species].create_dataset('ElectronicLevel', data=data)
        try:
            SPECIES_DICT_FLAG, HeatOfFormation_K, MassIC_kg, current_species_charge, interactionID = get_attr_values(SPECIES_DICT_FLAG, SPECIES_DICT, current_species, ATCT_URL)
        except:
            return SPECIES_DICT_FLAG
        create_attributes(dataset, current_species, ATCT_URL, MassIC_kg, HeatOfFormation_K, current_species_charge, interactionID)
    else:
        print("Species already exists in database\n")


###################################################################################################
#   functions for ATcT data handeling 
###################################################################################################
def create_species_dict(SPECIES_DICT, ATCT_URL):
    # Create an HTML session
    session = HTMLSession()
    # Make a GET request to the URL
    response = session.get(ATCT_URL)
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
        SPECIES_DICT.append({
        "Species": formula,
        "HeatOfFormation_K": delta_f_H_298K,
        "MassIC_kg": mass
    })
    # print(SPECIES_DICT) 
    return SPECIES_DICT

# Function to search for a species in the dictionary
def search_species_value(SPECIES_LABEL, SPECIES_DICT, ATCT_URL):
    # convert species name to fit SPECIES_DICT
    if 'Ion' in SPECIES_LABEL:
        Ion_Number =  int(re.search(r'\d+$', SPECIES_LABEL).group())
        if Ion_Number == 1:
            SPECIES_LABEL = re.sub('Ion1','+', SPECIES_LABEL)
        elif Ion_Number != None:
            SPECIES_LABEL = re.sub('Ion','+', SPECIES_LABEL)
    for species in SPECIES_DICT:
        try:
            if species["Species"] == SPECIES_LABEL:
                return species["HeatOfFormation_K"], species['MassIC_kg']
        except:
            return None  # Return None if species is not found

# function to get values of species out of dictionary
def get_attr_values(SPECIES_DICT_FLAG, SPECIES_DICT, CURRENT_SPECIES, ATCT_URL):
    # read-in heat of formation and mass from ATcT at 298.15K
    if not SPECIES_DICT_FLAG:
        # create dict with heat of formation and mass for faster read in for second species
        SPECIES_DICT_FLAG = True
        SPECIES_DICT = create_species_dict(SPECIES_DICT, ATCT_URL)
    try:
        HeatOfFormation_K, MassIC_kg = search_species_value(CURRENT_SPECIES, SPECIES_DICT, ATCT_URL)
    except:
        if CURRENT_SPECIES != 'electron':
            print(red('Error:')+' It seems like the species %s was not found at %s! Please check for the species on website\n' % (CURRENT_SPECIES, ATCT_URL))
        return SPECIES_DICT_FLAG
    # read-in ChargeIC from name of species
    if 'Ion' in CURRENT_SPECIES:
        Ion_Number =  int(re.search(r'\d+$', CURRENT_SPECIES).group())
        electron_charge = -1.60217653E-19 # coulombs
        current_species_charge = - Ion_Number * electron_charge
    elif 'electron' in CURRENT_SPECIES:
        current_species_charge = electron_charge
    else:
        current_species_charge = 0
    # get interactionID from current species name
    if sum(1 for c in CURRENT_SPECIES.replace('Ion','') if c.isupper()) == 1 and not bool(re.search(r'\d+', re.sub(r'Ion\d+','',CURRENT_SPECIES))):
        if not bool(re.search(r'[A-Za-z]*\d', re.sub(r'Ion\d+','',CURRENT_SPECIES))):
            if 'Ion' in CURRENT_SPECIES:
                interactionID = 10
            elif 'Ion' not in CURRENT_SPECIES:
                interactionID = 1
    elif bool(re.search(r'\d+', re.sub(r'Ion\d+','',CURRENT_SPECIES))) or sum(1 for c in CURRENT_SPECIES.replace('Ion','') if c.isupper()) != 1:
        if 'Ion' in CURRENT_SPECIES:
            interactionID = 20
        elif 'Ion' not in CURRENT_SPECIES:
            interactionID = 2
    elif 'electron' in CURRENT_SPECIES:
        interactionID = 4
    return SPECIES_DICT_FLAG, HeatOfFormation_K, MassIC_kg, current_species_charge, interactionID

# creating attrubutes for current species with data from dict
def create_attributes(DATASET, CURRENT_SPECIES, ATCT_URL, MassIC_kg, HeatOfFormation_K, current_species_charge, interactionID,old_attrs=[]):
    # get input from user -> maybe also database read in? where?
    if 'Tref' in list(old_attrs):
        DATASET.attrs['Tref'] = old_attrs['Tref']
    else:
        DATASET.attrs['Tref'] = float(input(bold('\nPlease enter Tref in [K] of current species %s\n-->') % CURRENT_SPECIES))
    
    if 'dref' in list(old_attrs):
        DATASET.attrs['dref'] = old_attrs['dref']
    else:
        DATASET.attrs['dref'] = float(input(bold('\nPlease enter dref in [m] of current species %s\n-->') % CURRENT_SPECIES))

    if 'omega' in list(old_attrs):
        DATASET.attrs['omega'] = old_attrs['omega']
    else:
        DATASET.attrs['omega'] = float(input(bold('\nPlease enter omega in [] of current species %s\n-->') % CURRENT_SPECIES))

    # set attributes
    DATASET.attrs['MassIC'] = MassIC_kg
    DATASET.attrs['HeatOfFormation_K'] = HeatOfFormation_K
    DATASET.attrs['ChargeIC'] = current_species_charge
    DATASET.attrs['InteractionID'] = interactionID

    array_of_strings = np.empty((6,2), dtype='S500')
    array_of_strings[0,0] = 'Electronic Energy Levels'
    array_of_strings[1,0] = 'HeatOfFormation_K'
    array_of_strings[2,0] = 'MassIC'
    array_of_strings[3,0] = 'Tref'
    array_of_strings[4,0] = 'dref'
    array_of_strings[5,0] = 'omega'

    array_of_strings[0,1] = 'Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2020). NIST Atomic Spectra Database (ver. 5.8), [Online]. Available: https://physics.nist.gov/asd. National Institute of Standards and Technology, Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F. Retrieved on ' + date.today().strftime("%B %d, %Y") + '.'
    array_of_strings[1,1] = 'B. Ruscic and D. H. Bross, Active Thermochemical Tables (ATcT) values based on ver. 1.130 of the Thermochemical Network. Argonne National Laboratory, Lemont, Illinois 2023; available at ATcT.anl.gov. DOI: https://doi.org/10.17038/CSE/1997229. Retrieved on ' + date.today().strftime("%B %d, %Y") + ' from ' + ATCT_URL + '.'
    array_of_strings[2,1] = 'B. Ruscic and D. H. Bross, Active Thermochemical Tables (ATcT) values based on ver. 1.130 of the Thermochemical Network. Argonne National Laboratory, Lemont, Illinois 2023; available at ATcT.anl.gov. DOI: https://doi.org/10.17038/CSE/1997229. Retrieved on ' + date.today().strftime("%B %d, %Y") + ' from ' + ATCT_URL + '.'
    array_of_strings[3,1] = 'User input at ' + date.today().strftime("%B %d, %Y")
    array_of_strings[4,1] = 'User input at ' + date.today().strftime("%B %d, %Y")
    array_of_strings[5,1] = 'User input at ' + date.today().strftime("%B %d, %Y")

    DATASET.attrs.create('* Reference', array_of_strings)
    DATASET.attrs['* Created']   = date.today().strftime("%B %d, %Y")

# checking attribute in existing database
def check_single_attr(HDF_UNIFIED_DATA, CURRENT_SPECIES, WANTED_VALUE, VALUE_NAME, ATCT_URL):
    try:
        if np.isclose(HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES].attrs[VALUE_NAME], WANTED_VALUE, rtol=1e-03, atol=1e-08):
            print('%s for ' % (VALUE_NAME) +green(CURRENT_SPECIES)+' is equal so will be kept')
        else:
            print(underlinE(VALUE_NAME)+' of species '+green(CURRENT_SPECIES)+'\n'+yellow('unified species database')+f': {HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES].attrs[VALUE_NAME]}'+"\n"+blue(ATCT_URL) +f': {WANTED_VALUE}\n')
            return -1
    except KeyError:
        if VALUE_NAME in ['HeatOfFormation_K','MassIC']:
            print("It seems like %s is not set for %s so it will be set as %s with data from "% (VALUE_NAME,CURRENT_SPECIES,WANTED_VALUE)+blue(ATCT_URL))
        elif VALUE_NAME in ['ChargeIC','InteractionID']:
            print("It seems like %s is not set for %s so it will be set as %s"% (VALUE_NAME,CURRENT_SPECIES,WANTED_VALUE))
        return VALUE_NAME
        

# checking all attributes in existing database
def check_attributes(HDF_UNIFIED_DATA, SPECIES_DICT_FLAG, SPECIES_DICT, CURRENT_SPECIES, ATCT_URL):
    try:
        SPECIES_DICT_FLAG, HeatOfFormation_K, MassIC_kg, current_species_charge, interactionID = get_attr_values(SPECIES_DICT_FLAG, SPECIES_DICT, CURRENT_SPECIES, ATCT_URL) 
    except:
        return SPECIES_DICT_FLAG, 1
    
    return_list = []

    if 'HeatOfFormation_K' not in HDF_UNIFIED_DATA['Species'][CURRENT_SPECIES].attrs:                                   # sanity check: calculate heat of formation for ions and compare with data from ATcT
        Ion_Number =  int(re.search(r'\d+$', CURRENT_SPECIES).group())
        ground_state = re.sub(r'Ion\d+','',CURRENT_SPECIES)
        HeatOfFormation_Sum = HDF_UNIFIED_DATA['Species'][ground_state].attrs['HeatOfFormation_K']+HDF_UNIFIED_DATA['Species'][ground_state]['ElectronicLevel'][:][-1,1]
        for i in range(1, Ion_Number):
            prev_state = ground_state + 'Ion' + str(i)
            HeatOfFormation_Sum = HeatOfFormation_Sum + HDF_UNIFIED_DATA['Species'][prev_state]['ElectronicLevel'][:][-1,1]

        if not np.isclose(HeatOfFormation_K, HeatOfFormation_Sum, rtol=1e-03, atol=0.0):          # HeatOfFormation_K being value from ATcT and HeatOfFormation_Sum cal value with relative tolerance of 1e-05
            print(red('Caution:')+' Please check the heat of formation of the current species: '+green(CURRENT_SPECIES)+', since it did not matcht with the value from '+blue(ATCT_URL)+'!\nIt is calculated in the PICLas source code like this: HeatOfFormation(CurrentSpecies) = HeatOfFormation(previous state) + Ionization Energy(last entry in electronic level dataset of previous state)')
            print(yellow("Calculated HeatOfFormation: ")+str(HeatOfFormation_Sum)+"\n"+blue("HeatOfFormation form ATcT: ")+str(HeatOfFormation_K))

    else:        # dont check heat of formation of ions because it is not stored
        return_list.append(check_single_attr(HDF_UNIFIED_DATA, CURRENT_SPECIES, HeatOfFormation_K, 'HeatOfFormation_K', ATCT_URL))
    
    return_list.append(check_single_attr(HDF_UNIFIED_DATA, CURRENT_SPECIES, MassIC_kg, 'MassIC', ATCT_URL))
    return_list.append(check_single_attr(HDF_UNIFIED_DATA, CURRENT_SPECIES, current_species_charge, 'ChargeIC', ATCT_URL))
    return_list.append(check_single_attr(HDF_UNIFIED_DATA, CURRENT_SPECIES, interactionID, 'InteractionID', ATCT_URL))
    if -1 in return_list:
        user_input = get_valid_input(bold('\nPlease enter') + "\n " + purple('1') + " to keep attributes from " + yellow('unified species database') + " or \n " + purple('2') + " to save attributes from " + blue(ATCT_URL) + " or \n " + purple('3') + " to exit program here\n-->", lambda x: x == '1' or x == '2' or x == '3')
        if user_input == '1':
            print("Keeping attributes for species ",green(CURRENT_SPECIES), "\n")
            return SPECIES_DICT_FLAG, 1
        elif user_input == '2':
            old_attrs = HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES].attrs
            data_array = HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES][:]
            del HDF_UNIFIED_DATA["Species"][CURRENT_SPECIES]
            dataset = HDF_UNIFIED_DATA['Species'].create_dataset(CURRENT_SPECIES, data=data_array)
            try:
                SPECIES_DICT_FLAG, HeatOfFormation_K, MassIC_kg, current_species_charge, interactionID = get_attr_values(SPECIES_DICT_FLAG, SPECIES_DICT, CURRENT_SPECIES, ATCT_URL)
            except TypeError:
                return SPECIES_DICT_FLAG, -1
            create_attributes(dataset, CURRENT_SPECIES, ATCT_URL, MassIC_kg, HeatOfFormation_K, current_species_charge, interactionID,old_attrs)
        elif user_input == '3':
            own_exit()

    elif any(item != -1 and item is not None for item in return_list):
        value_dict={
            'HeatOfFormation_K':HeatOfFormation_K,
            'MassIC':MassIC_kg,
            'ChargeIC':current_species_charge,
            'InteractionID':interactionID,
        }

        for VALUE_NAME in return_list:
            if VALUE_NAME != None:
                HDF_UNIFIED_DATA['Species'][CURRENT_SPECIES].attrs.create(VALUE_NAME, value_dict.get(VALUE_NAME) ,dtype='float')

    print("\n")
    return SPECIES_DICT_FLAG, 1

def check_ATcT(CURRENT_SPECIES, HDF_UNIFIED_DATA, SPECIES_DICT_FLAG, SPECIES_DICT, ATCT_URL):
    # check if species is already in Species database
    if CURRENT_SPECIES in HDF_UNIFIED_DATA['Species'].keys():
        SPECIES_DICT_FLAG, break_flag = check_attributes(HDF_UNIFIED_DATA, SPECIES_DICT_FLAG, SPECIES_DICT, CURRENT_SPECIES, ATCT_URL)
        if break_flag == -1:
            return SPECIES_DICT_FLAG,-1
        return SPECIES_DICT_FLAG, 1
    else:
        print("Species "+green(CURRENT_SPECIES)+" not found in database. Please create new species.")