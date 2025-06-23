import h5py
import numpy as np

###################################################################################################
# This file contains the configuration and global variables the script uses
###################################################################################################
# Read in Species Database
relative_path = '../../SpeciesDatabase.h5'
hdf_unified_data = h5py.File(relative_path, 'a')
# select ATcT Version (script build for 1.130)
version = '1.130'
ATcT_URL = f'https://atct.anl.gov/Thermochemical%20Data/version%20{version}'
# Base URL of the query
URL_base = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
electron_charge = -1.60217653E-19 # coulombs
# datatype for strings in database file
datatype_h5 = h5py.string_dtype(encoding='utf-8')

###################################################################################################
# This section contains the data types and needed attributes for the different types of species
# if a new attribute is added using this script it will be added to this file automatically
###################################################################################################

# setup for types of attributes
attribute_types = {
    '* Created'         : str,
    '* Reference'       : str,
    'PolyatomicMol'     : int,
    'InteractionID'     : int,
    'HeatOfFormation_K' : float,
    'Ediss_eV'          : float,
    'ChargeIC'          : float,
    'MassIC'            : float,
    'SymmetryFactor'    : float,
    'CharaTempRot'      : np.array,
    'CharaTempVib'      : np.array,
    'MomentOfInertia'   : np.array,
    'Tref'              : float,
    'dref'              : float,
    'omega'             : float,
    'Vib-OmegaE'        : float,
    'Vib-ChiE'          : float,
}

# setup of different empty classes
# these attributes include all attributes needed for the different types of species, which cannot be obtained from the name of the species
# attributes which can be obtained from the name of the species are directly written to the class in __init__ function of the class in edit_species.py
atom_attributes = {
    '* Created'         : None,
    '* Reference'       : None,
    'MassIC'            : None,
    'HeatOfFormation_K' : None,
    'Tref'              : None,
    'dref'              : None,
    'omega'             : None,
}

diatomic_attributes = {
    '* Created'         : None,
    '* Reference'       : None,
    'MassIC'            : None,
    'HeatOfFormation_K' : None,
    'Ediss_eV'          : None,
    'SymmetryFactor'    : None,
    'CharaTempRot'      : None,
    'CharaTempVib'      : None,
    'MomentOfInertia'   : None,
    'Tref'              : None,
    'dref'              : None,
    'omega'             : None,
    # AHO model attributes
    'Vib-OmegaE'        : None,
    'Vib-ChiE'          : None,
}

polyatomic_attributes = {
    'NewData'         : None,
    '* Created'         : None,
    '* Reference'       : None,
    'MassIC'            : None,
    'HeatOfFormation_K' : None,
    'Ediss_eV'          : None,
    'SymmetryFactor'    : None,
    'CharaTempRot'      : None,
    'CharaTempVib'      : None,
    'Tref'              : None,
    'dref'              : None,
    'omega'             : None,
    'PolyatomicMol'     : 1,
}