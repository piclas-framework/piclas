import re

import h5py
import lxcat_data_parser as ldp

# Input database
database_input = "Example.txt"
# Output database
database_output = "Example.h5"
# Species list to be included in the output database
species_list = ["Xe","Ar"]
# Reference of the utilized database
reference = "Biagi-v7.1 database, www.lxcat.net, retrieved on April 04, 2022. LXCat is an open-access website with databases contributed by members of the scientific community."

# Output of the HDF5 database
hdf = h5py.File(database_output, 'w')
hdf.attrs['Info'] = 'Cross-section database: First column is the collision energy in [eV], second column is the cross-section in [m^2]'

def sort_by_threshold(elem):
    return elem.threshold

# Read-in of data and output of the species containers
for current_species in species_list:
    print(f'\nCurrent species: {current_species}')
    grp_spec = hdf.create_group(current_species + "-electron")
    data_spec = ldp.CrossSectionSet(database_input,imposed_species=current_species)
    result_rot = 0
    result_vib = 0
    result_elec = 0
    total = 0
    for cross_section in data_spec.cross_sections:
        if cross_section.type == ldp.CrossSectionTypes.EXCITATION: #and cross_section.type != ldp.CrossSectionTypes.IONIZATION:
            # Electronic levels: Look for a star in the process description or scan the info for electronic
            if cross_section.info.get('PROCESS').find('*') != -1 or re.search(r'[A-Za-z]+\(([0-9]+[SPDF][0-9]*!*|HIGH)\)', cross_section.info.get('PROCESS')) !=1 or str(cross_section.info).casefold().find('electronic') != -1:
                result_elec += 1
            # Rotational excitation: Look for a j or rot in the process description or scan the info for rotation
            if cross_section.info.get('PROCESS').casefold().find('j') != -1 or cross_section.info.get('PROCESS').casefold().find('rot') != -1 or str(cross_section.info).casefold().find('rotation') != -1:
                result_rot += 1
            # Vibrational excitation: Look for v (but check that its not the v from eV) or scan the info for vibration
            if cross_section.info.get('PROCESS').casefold().find('v') != -1:
                if cross_section.info.get('PROCESS').casefold().find('v')-1 != cross_section.info.get('PROCESS').casefold().find('ev') or str(cross_section.info).casefold().find('vibration') != -1:
                    result_vib += 1
            elif str(cross_section.info).casefold().find('vibration') != -1:
                result_vib += 1
            total += 1
    if sum([result_rot,result_vib,result_elec]) != total:
        print('WARNING: '+ str(total-sum([result_rot,result_vib,result_elec])) +' excitation cross-section(s) could not be assigned!')
        print('WARNING: Cross-sections will be written in the UNDEFINED group.')
        grp_undef = grp_spec.create_group("UNDEFINED")
        grp_undef.attrs['Reference'] = reference
    else:
        print('Found:')
        print(str(result_rot)+' rotational cross-section(s).')
        print(str(result_vib)+' vibrational cross-section(s).')
        print(str(result_elec)+' electronic cross-section(s).')
    if result_rot > 0:
        grp_rot = grp_spec.create_group("ROTATION")
        grp_rot.attrs['Reference'] = reference
    if result_vib > 0:
        grp_vib = grp_spec.create_group("VIBRATION")
        grp_vib.attrs['Reference'] = reference
    if result_elec > 0:
        grp_elec = grp_spec.create_group("ELECTRONIC")
        grp_elec.attrs['Reference'] = reference
    for cross_section in data_spec.cross_sections:
        # Read-in of the effective cross-sections
        if cross_section.type == ldp.CrossSectionTypes.EFFECTIVE:
            ## Get the string of the cross section type (ELASTIC = 0, EFFECTIVE = 1, EXCITATION = 2, ATTACHMENT = 3, IONIZATION = 4)
            type_spec = ldp.CrossSectionTypes(cross_section.type).name
            ## Print name of the current species in console
            print('Found EFFECTIVE cross-section.')
            ## Write cross-section dataset of the current species in the HDF5 database
            dataset = grp_spec.create_dataset(type_spec, data=cross_section.data)
            ## Save the type of cross-section
            dataset.attrs['Type'] = type_spec
            ## Save the additional information and reference
            dataset.attrs['Info'] = str(cross_section.info)
            dataset.attrs['Reference'] = reference
        # Read-in of the elastic cross-sections
        if cross_section.type == ldp.CrossSectionTypes.ELASTIC:
            ## Get the string of the cross section type (ELASTIC = 0, EFFECTIVE = 1, EXCITATION = 2, ATTACHMENT = 3, IONIZATION = 4)
            type_spec = ldp.CrossSectionTypes(cross_section.type).name
            ## Print name of the current species in console
            print('Found ELASTIC cross-section.')
            ## Write cross-section dataset of the current species in the HDF5 database
            dataset = grp_spec.create_dataset(type_spec, data=cross_section.data)
            ## Save the type of cross-section
            dataset.attrs['Type'] = type_spec
            ## Save the additional information
            dataset.attrs['Info'] = str(cross_section.info)
            dataset.attrs['Reference'] = reference
    for cross_section in sorted([cross_section for cross_section in data_spec.cross_sections if cross_section.type == ldp.CrossSectionTypes.EXCITATION],key=sort_by_threshold):
        ## Write cross-section dataset of the current species in the HDF5 database
        test = 0
        # Electronic levels: Look for a star in the process description or scan the info for electronic
        if cross_section.info.get('PROCESS').find('*') != -1 or re.search(r'[A-Za-z]+\(([0-9]+[SPDF][0-9]*!*|HIGH)\)', cross_section.info.get('PROCESS')) !=1 or str(cross_section.info).casefold().find('electronic') != -1:
            dataset = grp_elec.create_dataset(str(cross_section.threshold), data=cross_section.data)
            test += 1
        # Rotational excitation: Look for a j or rot in the process description or scan the info for rotation
        if cross_section.info.get('PROCESS').casefold().find('j') != -1 or cross_section.info.get('PROCESS').casefold().find('rot') != -1 or str(cross_section.info).casefold().find('rotation') != -1:
            dataset = grp_rot.create_dataset(str(cross_section.threshold), data=cross_section.data)
            test += 1
        # Vibrational excitation: Look for v (but check that its not the v from eV) or scan the info for vibration
        if cross_section.info.get('PROCESS').casefold().find('v') != -1:
            if cross_section.info.get('PROCESS').casefold().find('v')-1 != cross_section.info.get('PROCESS').casefold().find('ev') or str(cross_section.info).casefold().find('vibration') != -1:
                dataset = grp_vib.create_dataset(str(cross_section.threshold), data=cross_section.data)
                test += 1
        elif str(cross_section.info).casefold().find('vibration') != -1:
            dataset = grp_vib.create_dataset(str(cross_section.threshold), data=cross_section.data)
            test += 1
        if test == 0:
            dataset = grp_undef.create_dataset(str(cross_section.threshold), data=cross_section.data)
        elif test > 1:
            print('WARNING: Cross-section data was added to multiple excitation types: '+str(cross_section.info))
        ## Get the string of the cross section type (ELASTIC = 0, EFFECTIVE = 1, EXCITATION = 2, ATTACHMENT = 3, IONIZATION = 4)
        type_spec = ldp.CrossSectionTypes(cross_section.type).name
        ## Save the type of cross-section
        dataset.attrs['Type'] = type_spec
        ## Save the additional information
        dataset.attrs['Info'] = str(cross_section.info)
        dataset.attrs['Threshold [eV]'] = cross_section.threshold
    for cross_section in sorted([cross_section for cross_section in data_spec.cross_sections if cross_section.type == ldp.CrossSectionTypes.IONIZATION],key=sort_by_threshold):
        ## Write cross-section dataset of the current species in the HDF5 database
        grp_reaction = grp_spec.create_group("REACTION")
        # split ionization into educts and products
        products = cross_section.info.get('PROCESS').split('->')[1]
        # cut off description at the end
        products = products.split(',')[0]
        # replace + sign (with space to avoid removing + from charge) with due to expected format
        products = products.replace(' + ','-')
        # replace single E with electron due to expected format
        products = re.sub(r'\bE\b', 'electron', products)
        if '+' in products:
            # search for digits in products to find charge of particle
            Ion_Number = re.search(r'\d+$', products)
            if Ion_Number is None:
                products = re.sub(r'\+','Ion1',products)
            else:
                Ion_Number = int(re.search(r'\d+$', products).group())
                products = re.sub(r'\+\d+$',f'Ion{Ion_Number}',products)
        # reorder the products to place each electron at the end
        product_parts = [p.strip() for p in products.split('-')]
        non_electrons = [p for p in product_parts if p != 'electron']
        electrons = [p for p in product_parts if p == 'electron']
        products = '-'.join(non_electrons + electrons)

        dataset = grp_reaction.create_dataset(str(products), data=cross_section.data)
        ## Get the string of the cross section type (ELASTIC = 0, EFFECTIVE = 1, EXCITATION = 2, ATTACHMENT = 3, IONIZATION = 4)
        type_spec = ldp.CrossSectionTypes(cross_section.type).name
        ## Save the type of cross-section
        dataset.attrs['Type'] = type_spec
        ## Save the additional information
        dataset.attrs['Info'] = str(cross_section.info)
        dataset.attrs['Threshold [eV]'] = cross_section.threshold
        print('Found IONIZATION cross-section.')

hdf.close()
