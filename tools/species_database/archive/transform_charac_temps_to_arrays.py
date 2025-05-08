import h5py
import re
import numpy as np

import sys
sys.path.append('../')
from edit_species import *
from general_functions import *

############################################################################################################################
# This script is not used any more, but it is kept for reference
# It is used to transform the attributes CharaTempRot, CharaTempVib, and MomentOfInertia from numbered attributes to arrays
############################################################################################################################

# Read in Species Database
relative_path = '../../../SpeciesDatabase.h5'
hdf_unified_data = h5py.File(relative_path, 'a')
# datatype for strings in database file
datatype_h5 = h5py.string_dtype(encoding='utf-8')

def transform_numbered_attrs_to_array(list_of_species_classes,base_name):
    # loop over all classes (equals all species in database)
    for species_class in list_of_species_classes:
        if isinstance(species_class,PolyatomicMolecule):
            if base_name == 'CharaTempRot' or base_name == 'MomentOfInertia':
                if species_class.attributes['LinearMolec'] == 1:
                    continue
            # Find all matching attributes
            pattern = f"^{base_name}\\d+$"
            matching_attrs = [attr for attr in species_class.attributes.keys()
                            if re.match(pattern, attr)]
            if not matching_attrs:
                print(f"No matching attributes found for {species_class.name} with {base_name}\n")
                continue

            # Sort by number
            matching_attrs.sort(key=lambda x: int(re.findall(r'\d+', x)[0]))

            # Collect values into array
            values = [species_class.attributes[attr] for attr in matching_attrs]
            array_data = np.array(values)
            # Create new attribute and delete old ones
            species_class.attributes[base_name] = array_data
            for attr in matching_attrs:
                del hdf_unified_data["Species"][species_class.name].attrs[attr]
            hdf_unified_data["Species"][species_class.name].attrs[base_name] = array_data
            print(f'Deleting {matching_attrs} for {species_class.name}, Added {base_name} to {species_class.name} as array')

if __name__ == '__main__':
    existing_species_list = list(hdf_unified_data["Species"].keys())
    list_of_species_classes = []
    # create class for each species
    for species in existing_species_list:
        instance = create_instance_from_data(species)
        list_of_species_classes.append(instance)

    change_attributes =['CharaTempRot','CharaTempVib','MomentOfInertia']
    for attr in change_attributes:
        print(green('Working on'),attr,'\n')
        transform_numbered_attrs_to_array(list_of_species_classes,attr)
        print(f'Done with {attr}!\n')