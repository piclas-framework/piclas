import numpy as np
import h5py
import lxcat_data_parser as ldp

hdf = h5py.File('LXCat_Database_Phelps_Electron_Scattering_EFFECTIVE.h5', 'w')
hdf.attrs['Info'] = 'Effective electron-neutrals collision cross-section database for read-in with PICLas. First column is the collision energy in [eV], second column is the cross-section in [m^2]'
hdf.attrs['Source'] = 'Phelps database, www.lxcat.net, retrieved on February 18, 2020.'

species_list = ["Ar","CO","CO2","H2","H2O","He","Mg","N2","NO","Na","Ne","O2","SF6"]

for current_species in species_list:
    print(current_species)
    data_spec = ldp.CrossSectionSet("Database_Phelps.txt",imposed_species=current_species)
    for cross_section in data_spec.cross_sections:
        # For now, only considering effective cross-sections in PICLas
        if cross_section.type == ldp.CrossSectionTypes.EFFECTIVE:
            ## Get the string of the cross section type (ELASTIC = 0, EFFECTIVE = 1, EXCITATION = 2, ATTACHMENT = 3, IONIZATION = 4)
            type_spec = ldp.CrossSectionTypes(cross_section.type).name
            ## Print name of the current species in console
            print(ldp.CrossSectionTypes(cross_section.type).name)
            ## Write cross-section dataset of the current species in the HDF5 database
            dataset = hdf.create_dataset(current_species + "-electron", data=cross_section.data)
            ## Save the type of cross-section
            dataset.attrs['Type'] = type_spec
            ## Save the additional information
            dataset.attrs['Info'] = str(cross_section.info)
            ## FUTURE DEVELOPMENTS
            ## Access a specific key of the info array
            # print(cross_section.info.get('SPECIES'))
            ## Save the mass ratio or threshold, depending on the type of the cross-section
            # if cross_section.type == ldp.CrossSectionTypes.EFFECTIVE or cross_section.type == ldp.CrossSectionTypes.ELASTIC:
            #     dataset.attrs['Mass ratio'] = cross_section.mass_ratio
            # if cross_section.type == ldp.CrossSectionTypes.IONIZATION or cross_section.type == ldp.CrossSectionTypes.EXCITATION:
            #     dataset.attrs['Threshold [eV]'] = cross_section.threshold

hdf.close()