import sys

def custom_sort_reactants(species_list):
    species_order = [
        'SF6',
        'SO2',
        'CH3CH3',
        'CH2CH2',
        'CHCH',
        'CH4',
        'CH3',
        'CH2',
        'CH',
        'CO2',
        'CO2Ion1',
        'CO',
        'COIon1',
        'HCN',
        'CN',
        'CNIon1',
        'C2',
        'C2Ion1',
        'N2',
        'N2Ion1',
        'N2O',
        'NO',
        'NOIon1',
        'NH3',
        'NH2',
        'NH',
        'H2O',
        'HO',
        'D2O',
        'DO',
        'O3',
        'O2',
        'O2Ion1',
        'H3',
        'H3Ion1',
        'H2',
        'H2Ion1',
        'D2',
        'D2Ion1',
        'HD',
        'HDIon1',
        'F2',
        'HF',
        'Cl2',
        'HCl',
        'Br2',
        'HBr',
        'C',
        'CIon1',
        'CIon2',
        'CIon3',
        'CIon4',
        'CIon5',
        'CIon6',
        'N',
        'NIon1',
        'O',
        'OIon1',
        'H',
        'HIon1',
        'D',
        'DIon1',
        'Xe',
        'XeIon1',
        'Ar',
        'ArIon1',
        'Kr',
        'Ne',
        'He',
        'HeIon1',
        'HeIon2',
        'Fe',
        'FeIon1',
        'electron',
        'M',
        'A'
        ]

    for val in species_list:
        if val not in species_order:
            print( 'The reactant ', val, 'is not defined as a species, please add it to the sorting function!')
            sys.exit()

    custom_order = {key:index + 1 for index, key in enumerate(species_order)}
    sorted_list = sorted(species_list, key=lambda x: custom_order[x])

    return sorted_list
    
    
def custom_sort_products(species_list):
    species_order = [
        'SF6',
        'SO2',
        'CH3CH3',
        'CH2CH2',
        'CHCH',
        'CH4',
        'CH3',
        'CH2',
        'CH',
        'CO2',
        'CO2Ion1',
        'CO',
        'COIon1',
        'HCN',
        'CN',
        'CNIon1',
        'C2',
        'C2Ion1',
        'N2',
        'N2Ion1',
        'N2O',
        'NO',
        'NOIon1',
        'NH3',
        'NH2',
        'NH',
        'H2O',
        'HO',
        'D2O',
        'DO',
        'O3',
        'O2',
        'O2Ion1',
        'H3',
        'H3Ion1',
        'H2',
        'H2Ion1',
        'D2',
        'D2Ion1',
        'HD',
        'HDIon1',
        'F2',
        'HF',
        'Cl2',
        'HCl',
        'Br2',
        'HBr',
        'C',
        'CIon1',
        'CIon2',
        'CIon3',
        'CIon4',
        'CIon5',
        'CIon6',
        'N',
        'NIon1',
        'O',
        'OIon1',
        'H',
        'HIon1',
        'D',
        'DIon1',
        'Xe',
        'XeIon1',
        'Ar',
        'ArIon1',
        'Kr',
        'Ne',
        'He',
        'HeIon1',
        'HeIon2',
        'Fe',
        'FeIon1',
        'electron',
        'M',
        'A'
        ]
    
    for val in species_list:
        if val not in species_order:
            print( 'The product ', val, 'is not defined as a species, please add it to the sorting function!')
            sys.exit()

    custom_order = {key:index + 1 for index, key in enumerate(species_order)}
    sorted_list = sorted(species_list, key=lambda x: custom_order[x])
    
    if 'M' in sorted_list:
        sorted_list.remove('M')
        sorted_list.insert(1,'M')
        
    if 'A' in sorted_list:
        sorted_list.remove('A')
        sorted_list.insert(1,'A')

    return sorted_list


