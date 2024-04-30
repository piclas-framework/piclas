###################################################################################################
#----------------------------------------- quick changes -----------------------------------------#
# Disclaimer: script only works on 'old' database from a27fcdae1f4d312a7f629273dc17ed8eade95769
###################################################################################################
# Name of Attribute of the Reaction #
reaction_attribute_name = 'ChemistryModel'
# this attribute is an array with dimension:    (X,1), where in the first column the Name of the Chemisty Model is stored 
#                                               or
#                                               (X,2), where in the first column the Name of the Chemisty Model is stored and in the second column the nonreactive species as a string seperated by ','
###################################################################################################
# Name of Attribute of Reactions (the group containing all reactions) #
all_reactions_attribute_name = 'AvailableChemistryModels,Reference'
###################################################################################################
# relative tolerance for arrhenius comparison to check if two reactions should be put together (keeping the first arrhenius coefficients)
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
    
    else:
        array_of_strings = np.empty((len(chem_model_list),1), dtype='S255')
        for i in range(len(chem_model_list)):
            Full_Chem_Model_List.append(chem_model_list[i][0])
            Full_Ref_List.append(str(reference_list[i]).encode('ascii', 'ignore').decode('utf-8'))
            array_of_strings[i,0] = chem_model_list[i][0]

    # delete old attributes
    if '* Reference' in hdf_unified_data["Reactions"][combined_reaction_name].attrs:
        del hdf_unified_data["Reactions"][combined_reaction_name].attrs['* Reference']
    if 'NonReactiveSpecies' in hdf_unified_data["Reactions"][combined_reaction_name].attrs:
        del hdf_unified_data["Reactions"][combined_reaction_name].attrs['NonReactiveSpecies']
    if 'ChemistryModel' in hdf_unified_data["Reactions"][combined_reaction_name].attrs:
        del hdf_unified_data["Reactions"][combined_reaction_name].attrs['ChemistryModel']
    hdf_unified_data['Reactions'][combined_reaction_name].attrs.create(reaction_attribute_name, array_of_strings)



# Read in Species Database
hdf_unified_data = h5py.File(relative_path, 'a')

# get list of all reactions
reactions_list = list(hdf_unified_data["Reactions"].keys())

# initialize variables
Full_Chem_Model_List = []
Full_Ref_List = []
Activation_Energy_Kr    = 0
Arrhenius_Prefactor     = 0
Arrhenius_Powerfactor   = 0

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
            
            else:
                for ArrString,OldValue in zip(["Arrhenius-Powerfactor", "Arrhenius-Prefactor" , "Activation-Energy_K"],[Arrhenius_Powerfactor, Arrhenius_Prefactor, Activation_Energy_Kr]):
                    ArrValue = hdf_unified_data["Reactions"][reaction].attrs[ArrString]
                    precision1 = len(str(ArrValue).split('.')[1]) if '.' in str(ArrValue) else 0
                    precision2 = len(str(OldValue).split('.')[1]) if '.' in str(OldValue) else 0
                    if precision1 > precision2:
                        hdf_unified_data["Reactions"][reaction].attrs[ArrString] = ArrValue
                    else:
                        hdf_unified_data["Reactions"][reaction].attrs[ArrString] = OldValue

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
        # fix for wrong value in old database
        if reaction == 'N+electron_NIon1+electron+electron#1' or reaction == 'N+electron_NIon1+electron+electron#2':
            del hdf_unified_data["Reactions"][combined_reaction_name].attrs['ReactionModel']
            hdf_unified_data['Reactions'][combined_reaction_name].attrs.create('ReactionModel', ['QK'],dtype='S255')
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
Full_Chem_Model_array = np.array(Full_Chem_Model_List, dtype='S255')
Full_Ref_array = np.array(Full_Ref_List, dtype='S255')
Full_Array = np.vstack((Full_Chem_Model_array, Full_Ref_array))
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
delete_list = []
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
            Activation_Energy_Kr.append(-1)
            Arrhenius_Prefactor.append(-1)
            Arrhenius_Powerfactor.append(-1)
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
            if Full_Array[i,-1] == Full_Array[k,-1]:
                if np.all(np.isclose(Full_Array[i][1:3].astype(float), Full_Array[k][1:3].astype(float), rtol=Rtol, atol=0.0)):
                    matching_rows.append((i, k))
                
    if matching_rows != 0:
        for match in matching_rows:
            # get matching idices
            index1, index2 = match
            current_models_array = hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[reaction_attribute_name]
            append_models_array = hdf_unified_data["Reactions"][Full_Array[index2,0]].attrs[reaction_attribute_name]

            if "Arrhenius-Powerfactor" in hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs and \
            "Arrhenius-Prefactor" in hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs and \
            "Activation-Energy_K" in hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs:
                
                for ArrString in ["Arrhenius-Powerfactor", "Arrhenius-Prefactor" , "Activation-Energy_K"]:
                    ArrValue1 = hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[ArrString]
                    ArrValue2 = hdf_unified_data["Reactions"][Full_Array[index2,0]].attrs[ArrString]
                    precision1 = len(str(ArrValue).split('.')[1]) if '.' in str(ArrValue1) else 0
                    precision2 = len(str(OldValue).split('.')[1]) if '.' in str(ArrValue2) else 0
                    if precision1 > precision2:
                        hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[ArrString] = ArrValue1
                    else:
                        hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[ArrString] = ArrValue2
            # stacked array for attribute creation
            new_models_array = np.vstack((current_models_array,append_models_array))

            # Store existing attributes
            attrs = hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs

            # keeping first match so only deleting attribute
            del hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs[reaction_attribute_name]
            # deleting whole reaction of second match
            delete_list.append(Full_Array[index2,0])
            
            hdf_unified_data["Reactions"][Full_Array[index1,0]].attrs.create(reaction_attribute_name, new_models_array)

delete_list = np.unique(delete_list)
for reaction in delete_list:
    del hdf_unified_data["Reactions"][reaction]

# get list of all reactions )stripped from '#01',etc.
unique_reaction_list = []
reaction_list = list(hdf_unified_data["Reactions"].keys())
for item in hdf_unified_data["Reactions"].keys():
    if re.sub(r'#0\d+','',item) not in unique_reaction_list:
        unique_reaction_list.append(re.sub(r'#0\d+','',item))
# counter for reactions in full reaction list (with '#01')
j = 0
for i,reaction in enumerate(unique_reaction_list):
    k = 1
    if j != len(unique_reaction_list) or not (re.sub(r'#0\d+','',reaction_list[j]) != re.sub(r'#0\d+','',reaction_list[j+1])):
        while j < len(reaction_list) and (reaction == re.sub(r'#0\d+','',reaction_list[j]) or reaction == re.sub(r'#\d+','',reaction_list[j])):
            attrs = hdf_unified_data["Reactions"][reaction_list[j]].attrs
            del hdf_unified_data["Reactions"][reaction_list[j]]
            new_name = re.sub(r'#0\d+','',reaction_list[j]) + f'#{k}'
            dataset = hdf_unified_data["Reactions"].create_dataset(new_name, data=0)
            # Restore attributes
            for attr_name, attr_value in attrs.items():
                dataset.attrs[attr_name] = attr_value
            j += 1
            k += 1

    else:
        attrs = hdf_unified_data["Reactions"][reaction].attrs
        del hdf_unified_data["Reactions"][reaction]
        new_name = re.sub(r'#0\d+','',reaction) + '#1'
        dataset = hdf_unified_data["Reactions"].create_dataset(new_name, data=0)
        # Restore attributes
        for attr_name, attr_value in attrs.items():
            dataset.attrs[attr_name] = attr_value
        j += 1

reaction_list = list(hdf_unified_data["Reactions"].keys())
for i,reaction in enumerate(reaction_list):
    if reaction in ['O2+electron_O2Ion1+electron+electron#1','O+electron_OIon1+electron+electron#1','CO+electron_COIon1+electron+electron#1']:
        hdf_unified_data["Reactions"][reaction].attrs['ReactionModel'] = 'QK'

print('Output successful.')
hdf_unified_data.close()