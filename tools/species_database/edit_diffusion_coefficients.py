import logging
import os
import io
import pandas as pd
import argparse
import h5py
import numpy as np
import re

from general_functions import *
from config import datatype_h5

def check_for_bolsig(file):
    with open(file, "r") as f:
        line = f.readline()
        while line:
            if "BOLSIG" in line.upper():
                return True
            line = f.readline()
        return False

def parse_selection(selection, total_files):
    selected_indices = set()
    parts = selection.split(',')
    for part in parts:
        if '-' in part:
            start, end = part.split('-')
            start, end = int(start.strip()), int(end.strip())
            selected_indices.update(range(start, end + 1))
        else:
            selected_indices.add(int(part.strip()))
    return [idx - 1 for idx in selected_indices if 0 <= idx - 1 < total_files]

def read_file_by_numbers(txt_files):
    while True:
        try:
            file_numbers = input("Enter the numbers corresponding to the files you want to read (comma separated, ranges allowed): ")
            selected_indices = parse_selection(file_numbers, len(txt_files))
            if selected_indices:
                selected_files = [txt_files[i] for i in selected_indices]
                print(f"Selected files: {', '.join(selected_files)}")
                break
            else:
                print("Invalid numbers. Please try again.")
        except ValueError:
            print("Invalid input. Please enter numbers separated by commas or ranges.")
    return selected_files

class DiffusionCoefficientsSetLXCatBolsig:
    """A class containing a set of diffusion coefficients.

    Attributes:
        species: combination of species as a str seperated by -, example: N2-electron
        config: parameters of diffusion coefficients
        datasetname: name of the dataset for species database
        reference: reference of database
        diff_coef: list of diffusion coefficients"""

    def __init__(self, input_file=None):
        """
        Reads a set of diffusion coefficients from a file.

        By default, reads the first diffusion coefficients set found in the input file.
        The input file should be compatible with the LXCat diffusion coefficients data format.
        """

        self.diff_coef = []
        coefficient_dict = {}
        if input_file is not None:
            with open(input_file, "r") as f:
                logging.info(f"Starting to read the contents "
                             f"of {os.path.basename(input_file)}")
                found_data = False
                reference = ''
                line = f.readline()
                while line:
                    if 'RECOMMENDED REFERENCE FORMAT' in line.upper():
                        line = f.readline()
                        while line.startswith('-'):
                            if "DATABASE" in line.upper():
                                self.database = line.split(' ')[1]
                            reference = reference + line.split('-')[1].strip()
                            line = f.readline()
                    if "CONFIG" in line.upper():
                        line = f.readline()
                        names = []
                        values = []
                        names.append(line.replace(' ','').split('=')[0])
                        values.append(line.replace(' ','').split('=')[1])
                        self.config = line.strip()
                        line = f.readline()
                        while '=' in line:
                            names.append(line.replace(' ','').split('=')[0])
                            values.append(line.replace(' ','').split('=')[1])
                            line = f.readline()
                        names_values_arr = np.zeros((len(names),2), dtype=datatype_h5)
                        names_values_arr[:,0] = names
                        names_values_arr[:,1] = values
                        self.config = names_values_arr
                    # find available data columns
                    if "transport coefficients" in line.lower():
                        line = f.readline()
                        while line.startswith('A'):
                            parts = line.split(maxsplit=1)
                            # check if line was split correctly
                            if len(parts) == 2:
                                coefficient_dict[parts[0]] = parts[1].strip()
                            line = f.readline()
                        # get index of columns
                        for col_name, col_value in coefficient_dict.items():
                            if col_value == 'Mobility *N (1/m/V/s)':
                                mobility_col_name = col_name
                            elif col_value == 'Diffusion coefficient *N (1/m/s)':
                                diffusion_col_name = col_name
                            elif col_value == 'Townsend ioniz. coef. alpha/N (m2)':
                                ionization_col_name = col_name
                    # check if data starts
                    if 'R# ' in line:
                        found_data = True
                    if found_data:
                        plots_list = []
                        mobility_coef = []
                        diffusion_coef = []
                        townsend_coef = []
                        parts_header = re.sub(r'\s+',' ', re.sub(r'\(.*?\)', '', line)).split(' ')
                        for i in range(len(parts_header)):
                            # used column names from above
                            if "E/N" in parts_header[i].strip():
                                plot_index = i
                            elif mobility_col_name == parts_header[i].strip():
                                mobility_index = i
                            elif diffusion_col_name == parts_header[i].strip():
                                diffusion_index = i
                            elif ionization_col_name == parts_header[i].strip():
                                townsend_index = i
                        line = f.readline()
                        while line.strip() != "":
                            if 'NaN' in line:
                                print("NaN in calculated values!")
                                exit()
                            parts = re.sub(r'\s+',' ', line).split(' ')
                            plots_list.append(parts[plot_index])
                            mobility_coef.append(parts[mobility_index])
                            diffusion_coef.append(parts[diffusion_index])
                            townsend_coef.append(parts[townsend_index])
                            line = f.readline()
                        # columns:  E/N (Td) | Mobility *N (1/m/V/s) | Diffusion coefficient *N (1/m/s) | Townsend ioniz. coef. alpha/N (m2)
                        data_array = np.zeros((len(plots_list),4))
                        data_array[:,0] = plots_list
                        data_array[:,1] = mobility_coef
                        data_array[:,2] = diffusion_coef
                        data_array[:,3] = townsend_coef
                        self.diff_coef = data_array
                        self.reference = reference
                        found_data = False
                    # for a mixture of species this part might be relevant and must be changed to read in parameters
                    if "RATE " in line.upper():
                        line = f.readline()
                        specList = []
                        specList.append(re.sub(r'\s+',' ', line.strip()).split(' ')[1])
                        while line.strip().startswith('C'):
                            spec = re.sub(r'\s+',' ', line.strip()).split(' ')[1]
                            if spec not in specList:
                                specList.append(spec)
                            line = f.readline()
                    line = f.readline()
                self.species = specList[0]
                # build name for PICLas Read-in
                if len(specList) != 1:
                    for i in range(len(specList)-1):
                        self.species = self.species + '-' + specList[i+1]
                self.species = self.species + '-electron'

                # add the format of plot axis (linear, exponential or quadratic)
                spacing_type = np.zeros((1,2), dtype=datatype_h5)
                spacing_type[0,0] = 'Point distribution'
                diff1 = np.diff(self.diff_coef[:,0])
                # Check linear spacing (constant differences)
                if np.allclose(diff1, diff1[0], rtol=1e-2):
                    spacing_type[0,1] = 'Linear'
                # Check exponential spacing (constant ratios)
                ratios = self.diff_coef[:,0][1:] / self.diff_coef[:,0][:-1]
                if np.allclose(ratios, ratios[0], rtol=1e-2):
                    spacing_type[0,1] = 'Exponential'
                # Check quadratic spacing (linear difference of differences)
                diff2 = np.diff(diff1)
                if np.allclose(diff2, diff2[0], rtol=1e-2):
                    spacing_type[0,1] = 'Quadratic'
                # append the spacing type to the config
                self.config = np.vstack((self.config, spacing_type))
            logging.info(f"Initialized {self}")


class DiffusionCoefficientsSetLXCat:
    """A class containing a set of diffusion coefficients.

    Attributes:
        species: combination of species as a str seperated by -, example: N2-electron
        config: parameters of diffusion coefficients
        database: name of the database
        datasetname: name of the dataset for species database
        reference: reference of database
        diff_coef: list of diffusion coefficients
        description: description of diffusion coefficients"""

    def __init__(self, input_file=None):
        """
        Reads a set of diffusion coefficients from a file.

        By default, reads the first diffusion coefficients set found in the input file.
        The input file should be compatible with the LXCat diffusion coefficients data format.
        """
        print(red('Caution! Functionality for datasets not calculated with BOLSIG+ is not properly tested yet!'))
        self.diff_coef = []

        if input_file is not None:
            with open(input_file, "r") as f:
                logging.info(f"Starting to read the contents "
                             f"of {os.path.basename(input_file)}")
                line = f.readline()
                while line:
                    if 'RECOMMENDED REFERENCE FORMAT' in line.upper():
                        line = f.readline()
                        reference = line.split('-')[1].strip()
                    # find the name of the database (optional)
                    if line.startswith("DATABASE:"):
                        self.database = line[9:].strip().replace(' ','-')
                        line = f.readline()
                    # find a line starting with one of the cross_section_types
                    found_species = line.startswith('SPECIES')
                    if found_species:
                        current_species = line.split(':')[1]
                        line = f.readline()
                        coef_type = line.split(':')[1]
                        if current_species.split('/')[0].strip() == 'e':
                            self.species = current_species.split('/')[1].strip() + '-electron'
                        else:
                            self.species = current_species.split('/')[1].strip() + '-' + current_species.split('/')[0].strip()
                        # the next lines may contain optional, additional
                        # information on the cross section with the format
                        # KEY: information
                        other_info = {}
                        line = f.readline()
                        while not line.startswith("-----"):
                            s = line.split(":")
                            key = s[0].strip()
                            other_info[key] = line[len(key) + 1:].strip()
                            line = f.readline()
                        # "-----" mars the start of the tabulated data
                        # put the data into an ioString
                        datasetname = coef_type.strip()
                        if 'IONIZATION' in datasetname.upper():
                            self.datasetname = 'IONIZATION'
                            self.description = 'Reduced Townsend coefficient (alpha/N) in m^2'
                        elif'MOBILITY' in datasetname.upper():
                            self.datasetname = 'MOBILITY'
                            self.description = 'Mobility x gas density (mu*N)'
                        elif'ENERGY' in datasetname.upper():
                            self.datasetname = 'ENERGY'
                            self.description = ' Characteristic energy (D/μ)'
                        elif'DIFFUSION' in datasetname.upper():
                            self.datasetname = 'DIFFUSION'
                            self.description = ' Diffusion x gas density (DN)'

                        if 'PARAM.' in other_info:
                            self.config = other_info['PARAM.'].strip()
                        else:
                            self.config = 'dummy'
                        data_stream = io.StringIO()
                        line = f.readline()
                        while not line.startswith("-----"):
                            data_stream.write(line)
                            line = f.readline()
                        data_stream.seek(0)
                        # "-----" marks the end of the tabulated data
                        # read the data into a pandas DataFrame
                        data = pd.read_csv(data_stream, sep="\t",
                                            names=[other_info['COLUMNS'].split('|')[0], other_info['COLUMNS'].split('|')[1]])
                        data_array = np.zeros((len(data),2))
                        data_array[:,0] = data[other_info['COLUMNS'].split('|')[0]].to_numpy()
                        data_array[:,1] =data[other_info['COLUMNS'].split('|')[1]].to_numpy()
                        self.diff_coef = data_array
                        red_field = data[other_info['COLUMNS'].split('|')[0]]
                        self.reference = reference #+ f'from {red_field.iloc[0]} to {red_field.iloc[-1]}'
                    line = f.readline()
        else:
            logging.info(f"Initialized {self}")

def append_to_database(database_path, Input):
    data = DiffusionCoefficientsSetLXCat(Input)
    # Open the HDF5 file in read/write mode
    with h5py.File(database_path, 'r+') as f:
        if 'Diffusion-Coefficients' not in f:
            f.create_group('Diffusion-Coefficients')
        diff_coef_group = f['Diffusion-Coefficients']

        speciesgroups = [ds for ds in diff_coef_group if isinstance(diff_coef_group[ds], h5py.Group)]
        if data.species not in speciesgroups:
            diff_coef_group.create_group(data.species)

        # List all datasets in the diff_coef_group
        databasegroups = [ds for ds in diff_coef_group[data.species] if isinstance(diff_coef_group[data.species][ds], h5py.Group)]

        if data.database not in databasegroups:
            diff_coef_group[data.species].create_group(data.database)
        else:
            datasets = [ds for ds in diff_coef_group[data.species][data.database] if isinstance(diff_coef_group[data.species][data.database][ds], h5py.Dataset)]
            if data.datasetname in datasets:
                print(f'Caution! Current dataset is already stored in species database for species: {data.species}')
                # maybe select what should be done -> override or keep

                exit()

        diff_coef_group[data.species][data.database].create_dataset(data.datasetname,data=data.diff_coef)
        diff_coef_group[data.species][data.database][data.datasetname].attrs.create('* Reference', data.reference)
        diff_coef_group[data.species][data.database][data.datasetname].attrs.create('* Description', data.description)
        if data.config != None:
            diff_coef_group[data.species][data.database].attrs.create('* Parameter', data.config)


def append_to_database_bolsig(database_path, Input):
    data = DiffusionCoefficientsSetLXCatBolsig(Input)
    # Open the HDF5 file in read/write mode
    with h5py.File(database_path, 'r+') as f:
        if 'Diffusion-Coefficients' not in f:
            f.create_group('Diffusion-Coefficients')
        diff_coef_group = f['Diffusion-Coefficients']

        # speciesgroups is list of all groups of Diffusion-Coefficients, e.g. 'Ar-electron', 'N2-electron'
        speciesgroups = [ds for ds in diff_coef_group if isinstance(diff_coef_group[ds], h5py.Group)]
        if data.species not in speciesgroups:
            diff_coef_group.create_group(data.species)

        # List all datasets in the diff_coef_group (name of the database where the swarm data was retrieved from), e.g. 'Biagi-v7.1', 'Phelps'
        databasegroups = [ds for ds in diff_coef_group[data.species] if isinstance(diff_coef_group[data.species][ds], h5py.Group)]

        # if not existent, create new group for database
        if data.database not in databasegroups:
            diff_coef_group[data.species].create_group(data.database)
        else:
            # check attributes of existing databases
            for database in databasegroups:
                if diff_coef_group[data.species][database].attrs['* Parameters'].all() == data.config.all():
                    print(red('Caution!')+f' Current dataset is already stored in species database for species: {data.species}, with parameters:')
                    for i in range(data.config.shape[0]):
                        print(f'{data.config[i, 0]}: {data.config[i, 1]}')
                    user_input_ = get_valid_input(create_prompt('to keep the current datasets',
                                                       'to overwrite the datasets',),
                                                       lambda x: x == '1' or x == '2' or x == '3')
                    if user_input_ == '1':
                        return
                    elif user_input_ == '2':
                        del diff_coef_group[data.species][data.database]
                    elif user_input_ == '3':
                        own_exit()

        DSetNameList = ['MOBILITY', 'DIFFUSION', 'IONIZATION']
        DescriptionList = ['Mobility x gas density (mu*N)', ' Diffusion x gas density (DN)', 'Reduced Townsend coefficient (alpha/N) in m^2']

        for i in range(3):
            data_array = np.vstack((data.diff_coef[:,0],data.diff_coef[:,i+1]))
            diff_coef_group[data.species][data.database].create_dataset(DSetNameList[i],data=data_array.T)
            diff_coef_group[data.species][data.database][DSetNameList[i]].attrs.create('* Reference', data.reference)
            diff_coef_group[data.species][data.database][DSetNameList[i]].attrs.create('* Description', DescriptionList[i])

        diff_coef_group[data.species][data.database].attrs.create('* Parameters', data.config)


if __name__ == "__main__":
     # read in all arguments
    parser = argparse.ArgumentParser(description='LxCat Diffusion Coefficients')
    parser.add_argument('InputFile',nargs='+', type=str, help='')
    args = parser.parse_args()
    database_PATH = '../../SpeciesDatabase.h5'
    for file in args.InputFile:
        bolsig_bool = check_for_bolsig(file)
        if bolsig_bool == True:
            append_to_database_bolsig(database_PATH, file)
        else:
            append_to_database(database_PATH, file)
    print("Output sucessful!")
