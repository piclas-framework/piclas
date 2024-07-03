import logging
import os
import io
import pandas as pd
import argparse
import h5py
import numpy as np
import re

from general_functions import *

def check_for_bolsig(file):
    with open(file, "r") as f:
        line = f.readline()
        while line:
            if "BOLSIG" in line.upper():
                return True
            line = f.readline()
        return False


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

        if input_file is not None:
            with open(input_file, "r") as f:
                logging.info(f"Starting to read the contents "
                             f"of {os.path.basename(input_file)}")
                found_data = False
                reference = ''
                line = f.readline()
                while line:
                    if 'RECOMMENDED REFERENCE FORMAT' in line:
                        line = f.readline()
                        while line.startswith('-'):
                            if "DATABASE" in line.upper():
                                self.database = line.split(' ')[1]
                            reference = reference + line.split('-')[1].strip()
                            line = f.readline()
                    # find the name of the database (optional)
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
                        names_values_arr = np.zeros((len(names),2), dtype='S255')
                        names_values_arr[:,0] = names
                        names_values_arr[:,1] = values
                        self.config = names_values_arr
                    if 'A29' in line:
                        found_data = True
                    # line = f.readline()
                    if found_data:
                        plots_list = []
                        A1 = []
                        A11 = []
                        A12 = []
                        A18 = []
                        line = f.readline()
                        line = f.readline()
                        while line.strip() != "":
                            if 'NaN' in line:
                                print("NaN in calculated values!")
                                exit()
                            parts = re.sub(r'\s+',' ', line).split(' ')
                            plots_list.append(parts[1])
                            A1.append(parts[2])
                            A11.append(parts[5])
                            A12.append(parts[6])
                            A18.append(parts[10])
                            line = f.readline()
                        # columns:  E/N (Td) | Mean energy (eV) | Energy mobility *N (1/m/V/s) | Energy diffusion coef. D*N (1/m/s) | Townsend ioniz. coef. alpha/N (m2)
                        data_array = np.zeros((len(plots_list),4))
                        data_array[:,0] = plots_list
                        # data_array[:,1] = A1
                        data_array[:,1] = A11
                        data_array[:,2] = A12
                        data_array[:,3] = A18
                        self.diff_coef = data_array
                        self.reference = reference #+ f'from {red_field.iloc[0]} to {red_field.iloc[-1]}'
                        found_data = False
                    if "RATE " in line.upper():
                        line = f.readline()
                        specList = []
                        specList.append(re.sub(r'\s+',' ', line).split(' ')[2])
                        while line.strip().startswith('C'):
                            line = f.readline()
                            spec = re.sub(r'\s+',' ', line).split(' ')[2]
                            if spec not in specList:
                                specList.append(spec)
                    line = f.readline()
                if "/" in str(specList[-1]):
                    specList = specList[:-1]
                if len(specList) == 1:
                    self.species = specList[0] + '-electron'
                    return
                self.species = specList[0]
                for i in range(len(specList)-1):
                    print(specList[i+1])
                    self.species = self.species + '-' + specList[i+1]
        else:
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

        self.diff_coef = []

        if input_file is not None:
            with open(input_file, "r") as f:
                logging.info(f"Starting to read the contents "
                             f"of {os.path.basename(input_file)}")
                line = f.readline()
                while line:
                    if 'RECOMMENDED REFERENCE FORMAT' in line:
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

        speciesgroups = [ds for ds in diff_coef_group if isinstance(diff_coef_group[ds], h5py.Group)]
        if data.species not in speciesgroups:
            diff_coef_group.create_group(data.species)

        # List all datasets in the diff_coef_group
        databasegroups = [ds for ds in diff_coef_group[data.species] if isinstance(diff_coef_group[data.species][ds], h5py.Group)]
        
        if data.database not in databasegroups:
            diff_coef_group[data.species].create_group(data.database)
        else:
            datasets = [ds for ds in diff_coef_group[data.species][data.database] if isinstance(diff_coef_group[data.species][data.database][ds], h5py.Dataset)]
            # check if dataset already in database
            if data.database in datasets:
                print(f'Caution! Current dataset is already stored in species database for species: {data.species}')
                # maybe select what should be done -> override or keep
                exit()

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
    