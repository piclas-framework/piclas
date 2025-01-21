from general_functions import *
from datetime import date
import re
import numpy as np
from config import *

def display_reaction(CURRENT_REACTION):
    print("\n"+red("Reaction: ") + green(CURRENT_REACTION))
    for AttrName, AttrValue in list(hdf_unified_data["Reactions"][CURRENT_REACTION].attrs.items()):
        if AttrName == 'ChemistryModel':
            for i in range(len(AttrValue[:,0])):
                print('\nChemistry Model %s: '%(i+1), AttrValue[i,0].decode('utf-8'))
                if AttrValue.shape[1] != 1:
                    print('Non-Reactive species: ', AttrValue[i,1].decode('utf-8'))
            print('')

        elif AttrName in ['Products','Reactants','ReactionModel']:
            if len(AttrValue) == 1:
                print(AttrName,": ",AttrValue[0].decode('utf-8'))
            else:
                print(AttrName, end=': ')
                i = 0
                for i in range(len(AttrValue)-1):
                    print(AttrValue[i].decode('utf-8'), end=',')
                print(AttrValue[i+1].decode('utf-8'))


        else:
            try:
                print(AttrName,': ',AttrValue.decode('utf-8'))
            except Exception:
                print(AttrName,': ',AttrValue)

def check_reaction(CURRENT_REACTION):
    existing_reaction_list = list(hdf_unified_data["Reactions"].keys())
    unique_reaction_list = []
    for item in hdf_unified_data["Reactions"].keys():
        if re.sub(r'#\d+','',item) not in unique_reaction_list:
            unique_reaction_list.append(re.sub(r'#\d+','',item))

    if CURRENT_REACTION in unique_reaction_list:
        print("It seems like the given reaction exists already in the database - "+green(CURRENT_REACTION))
        # count matches for new reaction number later
        count = 0
        #check for more than one hits of current reaction
        for REACTION in existing_reaction_list:
            if re.sub(r'#\d+','',REACTION) == CURRENT_REACTION:
                display_reaction(REACTION)
                count += 1
        user_input = get_valid_input(create_prompt('add a new reaction with different attributes',
                                                   'add a new chemistry model to existing reaction',
                                                   'to skip this reaction'),
                                                   lambda x: x == '1' or x == '2' or x == '3' or x == '4')

        if user_input == '1':
            New_Reaction_Name = CURRENT_REACTION + "#" + str(count+1)
            return New_Reaction_Name
        elif user_input == '2':
            # new chem model to existing reaction
            if count != 1:
                user_input = input(bold('\nPlease enter the number of the reaction where you want add a new chemistry model\n-->'))
                chem_model_reaction_name = CURRENT_REACTION + "#" + str(user_input)
            else:
                chem_model_reaction_name = CURRENT_REACTION + "#1"
            only_new_chem_model(chem_model_reaction_name)
        elif user_input == '3':
            return None
        elif user_input == '4':
            own_exit()
    else:
        New_Reaction_Name = CURRENT_REACTION + "#1"
        return New_Reaction_Name

# create new chemistry model in the reaction attribute 'AvailableChemistryModels,Reference'
def create_new_ref_attribute(reaction_name, CHEM_MODEL):
    # add reference to attribute of reactions group
    Reference = input(bold('\nPlease enter the Reference of current reaction %s\n-->') % reaction_name)
    AvailChemModels = hdf_unified_data["Reactions"].attrs['AvailableChemistryModels,Reference']
    AvailChemModels_unicode = np.char.decode(AvailChemModels, 'utf-8')

    if not (np.any(AvailChemModels_unicode[:, 0] == CHEM_MODEL)):
        size = len(AvailChemModels_unicode[:, 0])
        new_chem_models = np.zeros((size+1,2), dtype=datatype_h5)
        for i in range(size):
            new_chem_models[i,0] = AvailChemModels_unicode[i,0]
            new_chem_models[i,1] = AvailChemModels_unicode[i,1]
        new_chem_models[-1,0] = CHEM_MODEL
        new_chem_models[-1,1] = Reference

        del hdf_unified_data["Reactions"].attrs['AvailableChemistryModels,Reference']
        hdf_unified_data['Reactions'].attrs.create('AvailableChemistryModels,Reference', new_chem_models ,dtype=datatype_h5)

def create_reaction(CURRENT_REACTION):
    reaction_name = check_reaction(CURRENT_REACTION)
    if reaction_name != None:
        # get attributes
        try:
            Reactants, Products = CURRENT_REACTION.split('_')
        except Exception as e:
            print("Read in of products and reactants from reaction name has failed with error %s. Please ensure reaction has correct format: C2+M_C+M+C (products and reactants separated by underscore). %s will be skipped"%(e,reaction_name))
            return
        Products = Products.split('+')
        prodList = remove_from_list(Products, 'M', '+', 'A')
        Reactants = re.sub(r'#\d+','',Reactants)
        Reactants = Reactants.split('+')
        reacList = remove_from_list(Reactants, 'M', '+', 'A')


        # create reaction
        dataset = hdf_unified_data["Reactions"].create_dataset(reaction_name, data=0)
        ReactionModel = get_valid_input(bold('\nPlease enter the reaction model of current reaction %s (Options: "QK" or "TCE")\n-->') % reaction_name, lambda x: x == 'QK' or x == 'TCE')
        dataset.attrs.create('ReactionModel', [ReactionModel], dtype=datatype_h5)
        if ReactionModel == 'TCE':          # get other parameter
            for ArrString in ["Arrhenius-Powerfactor", "Arrhenius-Prefactor" , "Activation-Energy_K"]:
                ArrValue = input(bold('\nPlease enter the '+ArrString+' of current reaction %s\n-->') % reaction_name)
                dataset.attrs.create(ArrString, ArrValue, dtype='float64')
        elif ReactionModel == 'QK':         # qk only used for ionization or dissociation
            if any('Ion' in element for element in prodList):       # check if ion is in prodcuts of reaction -> ionization
                # sanity check if required values are set for species: last electronic level with degeneracy 1
                for spec in reacList:
                    if spec != 'el':  # check all reactants except electrons
                        if hdf_unified_data['Species'][spec][:][-1,0] != 1:
                            print("Error: Degeneracy of last electronic level of species %s not equal to 1! Reaction will be skipped" % spec)
                            del hdf_unified_data["Reactions"][reaction_name]
                            return
                        else:
                            print("Please double check the ionization energy of species %s: %s" % (spec,hdf_unified_data['Species'][spec][:][-1,1]))
            else:       # dissociation reaction
                # sanity check of reactants: EdissEV in species set
                for spec in reacList:
                    if spec != 'el':  # check all reactants except electrons
                        if hdf_unified_data['Species'][spec].attrs['Ediss_eV'] == None:
                            print(red("Error:")+" Ediss_eV of species %s not found in database! Reaction will be skipped" % spec)
                            del hdf_unified_data["Reactions"][reaction_name]
                            return
                        else:
                            print("Please double check the dissociation energy of species %s: %s" % (spec,hdf_unified_data['Species'][spec].attrs['Ediss_eV']))

        ##### TODO: show available chem models and be able to select one
        ChemModel = input(bold('\nPlease enter the chemistry model to which the currenet reaction %s belongs\n-->') % reaction_name)
        if 'M' in Reactants or 'M' in Products or 'A' in Reactants or 'A' in Products:
            NonReactives = input(bold('\nPlease enter the non reactive species of currenet reaction %s as comma separated string, e.g "H2,H,O"\n-->') % reaction_name)
            dataset.attrs.create('ChemistryModel', [[ChemModel,NonReactives]], dtype=datatype_h5)
        else:
            dataset.attrs.create('ChemistryModel', [[ChemModel]], dtype=datatype_h5)

        dataset.attrs.create('* Created', date.today().strftime("%B %d, %Y"), dtype=datatype_h5)
        dataset.attrs.create('Products', np.array(prodList), dtype=datatype_h5)
        dataset.attrs.create('Reactants', np.array(reacList), dtype=datatype_h5)

        create_new_ref_attribute(reaction_name, ChemModel)

        print("New reaction created!")

def delete_reaction(CURRENT_REACTION):
    existing_reaction_list = list(hdf_unified_data["Reactions"].keys())
    unique_reaction_list = []
    for item in hdf_unified_data["Reactions"].keys():
        if re.sub(r'#\d+','',item) not in unique_reaction_list:
            unique_reaction_list.append(re.sub(r'#\d+','',item))

    if CURRENT_REACTION in unique_reaction_list:
        print("It seems likt the given reaction exists already in the database - "+green(CURRENT_REACTION))
        # counter for break if only one reaction exists
        counter = 0
        #check for more than one hits of current reaction
        for REACTION in existing_reaction_list:
            if re.sub(r'#\d+','',REACTION) == CURRENT_REACTION:
                counter += 1

        chemModelCheck = [] # save chemistry models of deleted reactions for check later
        if counter == 1:
            delete_reaction_name = CURRENT_REACTION + "#1"
            display_reaction(delete_reaction_name)
            user_input = get_valid_input(create_prompt('to delete this reaction',
                                                       'to skip this reaction'),
                                                       lambda x: x == '1' or x == '2' or x == '3')
            if user_input == '1':
                try:        # catch reaction falsly without chem model to still delete
                    chemModelCheck.append(hdf_unified_data["Reactions"][delete_reaction_name].attrs['ChemistryModel'][:,0])
                except:
                    pass
                del hdf_unified_data["Reactions"][delete_reaction_name]
            else:
                return

        else:
            # display other reactions
            for REACTION in existing_reaction_list:
                if re.sub(r'#\d+','',REACTION) == CURRENT_REACTION:
                    display_reaction(REACTION)
            # get which reaction to delete
            user_input = input(bold('\nPlease enter the number(s) of the reaction(s) you want to delete as comma separated string, e.g. 1,2,3\n-->'))
            user_input = user_input.split(',')
            # create reaction names to delete
            for num in user_input:
                delete_reaction_name = CURRENT_REACTION + "#" + str(num)
                try:        # catch reaction falsly without chem model to still delete
                    chemModelCheck.append(hdf_unified_data["Reactions"][delete_reaction_name].attrs['ChemistryModel'][:,0])
                except:
                    pass
                del hdf_unified_data["Reactions"][delete_reaction_name]

            # get new list without deleted reaction
            existing_reaction_list = list(hdf_unified_data["Reactions"].keys())
            # renumber remaining reactions
            counter = 1
            for REACTION in existing_reaction_list:
                if re.sub(r'#\d+','',REACTION) == CURRENT_REACTION:
                    attrs = hdf_unified_data["Reactions"][REACTION].attrs
                    del hdf_unified_data["Reactions"][REACTION]
                    new_name = CURRENT_REACTION + f'#{counter}'
                    dataset = hdf_unified_data["Reactions"].create_dataset(new_name, data=0)
                    # Restore attributes
                    for attr_name, attr_value in attrs.items():
                        dataset.attrs[attr_name] = attr_value
                    counter += 1

        # check if reaction with chem model still exists or if it should be deleted from the reactions attribute
        existing_reaction_list = list(hdf_unified_data["Reactions"].keys())
        # loop over all reactions
        availChemModels = hdf_unified_data["Reactions"].attrs['AvailableChemistryModels,Reference']
        for MODEL in chemModelCheck:
            for REATION in existing_reaction_list:
                if MODEL in hdf_unified_data["Reactions"][REATION].attrs['ChemistryModel'][:,0]:
                    break
            # model in remaining reactions not found - will be deleted from available chemistry models
            index_match = np.where(availChemModels[:, 0] == MODEL)[0]

            NewAvailChemModels = np.delete(availChemModels, index_match, axis=0)

        if len(chemModelCheck) != 0:    # catch reaction falsly without chem model to still delete
            del hdf_unified_data["Reactions"].attrs['AvailableChemistryModels,Reference']
            hdf_unified_data['Reactions'].attrs.create('AvailableChemistryModels,Reference', NewAvailChemModels ,dtype=datatype_h5)

    else:
        print("It seems like the reaction is not found in the database. Please ensure the right syntax (e.g. C+N_CNIon1+electron) or check by hand if the reaction exists.")

def only_new_chem_model(chem_model_reaction_name):
    # get chem model (and non reactives from reaction)
    ChemModelNonReacArray = hdf_unified_data["Reactions"][chem_model_reaction_name].attrs['ChemistryModel']
    ChemModelNonReacArray_unicode = np.char.decode(ChemModelNonReacArray, 'utf-8')
    size = len(ChemModelNonReacArray_unicode[:, 0])
    NumCols = ChemModelNonReacArray_unicode.shape[1]
    New_ChemModelNonReacArray = np.zeros((size+1,NumCols), dtype=datatype_h5)

    for i in range(size):
        New_ChemModelNonReacArray[i,0] = ChemModelNonReacArray_unicode[i,0]
        if NumCols == 2:
            New_ChemModelNonReacArray[i,1] = ChemModelNonReacArray_unicode[i,1]

    # add new chem model and non reacs
    ChemistryModel = input(bold('\nPlease enter the name of the new chemistry model\n-->'))
    New_ChemModelNonReacArray[-1,0] = ChemistryModel
    if NumCols == 2:
        NewNonReactives = input(bold('\nPlease enter the non reactive species of the new chemistry model as comma separated string, e.g "H2,H,O"\n-->'))
        New_ChemModelNonReacArray[-1,1] = NewNonReactives

    del hdf_unified_data["Reactions"][chem_model_reaction_name].attrs['ChemistryModel']
    hdf_unified_data['Reactions'][chem_model_reaction_name].attrs.create('ChemistryModel', New_ChemModelNonReacArray ,dtype=datatype_h5)

    create_new_ref_attribute(chem_model_reaction_name, ChemistryModel)
