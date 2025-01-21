import h5py
import os
import ast
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
    'CharaTempRot'      : float,
    'CharaTempRot1'     : float,
    'CharaTempRot2'     : float,
    'CharaTempRot3'     : float,
    'CharaTempVib'      : float,
    'CharaTempVib1'     : float,
    'CharaTempVib2'     : float,
    'CharaTempVib3'     : float,
    'CharaTempVib4'     : float,
    'CharaTempVib5'     : float,
    'CharaTempVib6'     : float,
    'CharaTempVib7'     : float,
    'CharaTempVib8'     : float,
    'CharaTempVib9'     : float,
    'MomentOfInertia'   : float,
    'MomentOfInertia1'  : float,
    'MomentOfInertia2'  : float,
    'MomentOfInertia3'  : float,
    'Tref'              : float,
    'dref'              : float,
    'omega'             : float,
    'Vib-OmegaE'        : float,
    'Vib-ChiE'          : float,
}