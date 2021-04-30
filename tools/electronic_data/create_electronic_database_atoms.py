import requests
import h5py
import pandas as pd
import io
import re
from datetime import date

# Function to determine the Roman number from an integer
def int_to_Roman(num):
   val = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   syb = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   roman_num = ""
   for i in range(len(val)):
      count = int(num / val[i])
      roman_num += syb[i] * count
      num -= val[i] * count
   return roman_num

### USER INPUT BEGIN
# Output of the HDF5 database
hdf = h5py.File('Electronic-State-Database.h5', 'a')
# Species list to be included in the output database
# species_list = ['N','O','NIon1','OIon1','H','Xe','XeIon1','XeIon2','Ar','ArIon1','ArIon2']
species_list = ['OIon1']
### USER INPUT END

# Base URL of the query
URL_base = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
# Loop over species list
for current_species in species_list:
  # Get the ionization number as the last digits and 1 to comply with the NIST standard
  ion_level = int(re.sub('.*?([0-9]*)$',r'\1',current_species) or 0) + 1
  # Get the species (everything before Ion)
  species = current_species.split('Ion')[0]
  # Build the species in the NIST database query format: Xe+I
  current_species_NIST = species + '+' + int_to_Roman(ion_level)
  print('Converted '+current_species+' to '+current_species_NIST+' for download from NIST database.')
  # Build the request URL
  URL_request = URL_base + '?de=0&spectrum=' + current_species_NIST + '&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&level_out=on&g_out=on&temp=&submit=Retrieve+Data'
  response = requests.get(URL_request)
  # Check the website
  if response.status_code == 200:
    print('Download successful.')
  # Get the data from the NIST website
  data = response.content.decode('utf-8')
  # Clean the data
  data = data.replace('[','')
  data = data.replace(']','')
  data = data.replace('(','')
  data = data.replace(')','')
  data = data.replace('= ','')
  data = data.replace('=','')
  data = data.replace(' ','')
  data = data.replace('"','')
  # Read-in data as a pandas dataframe
  data = pd.read_csv(io.StringIO(data),skipinitialspace=True,delimiter=",",usecols=['g','Levelcm-1'], na_values=['---'])
  # Drop rows with an empty cell and '---' entries (which were converted to N/A during read-in)
  data.dropna(inplace=True)
  # Check the datatype: if its a float, then all non-numerical characters have already been removed
  if data['Levelcm-1'].dtype != 'float64':
    # Drop rows with a question mark ("This level/line may not be real.")
    data.drop(data[data['Levelcm-1'].str.contains(r'[?]')].index,inplace=True)
    # Drop rows with a +x ("The relative positions of the levels within such a system are accurate within experimental uncertainties, but no experimental connection between this system and the other levels of the spectrum has been made.")
    data.drop(data[data['Levelcm-1'].str.contains(r'[+x]')].index,inplace=True)
  # Convert 1/cm to K
  data['Levelcm-1'] = data['Levelcm-1'].astype(float)
  data['Levelcm-1'] = 100 * data['Levelcm-1'] * 1.986E-025 / 1.38065E-023           # 1/cm * 100cm/m * (J m) / (J/K) = K
  data.columns=['Degeneracy', 'Level(K)']
  # Write to hdf: If data set already exists, delete the old set first
  if current_species in hdf.keys():
    del hdf[current_species]
    print('Old dataset replaced.')
  dataset = hdf.create_dataset(current_species, data=data)
  # Write attributes for source and time of retrieval
  dataset.attrs['Reference'] = 'Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2020). NIST Atomic Spectra Database (ver. 5.8), [Online]. Available: https://physics.nist.gov/asd. National Institute of Standards and Technology, Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F. Retrieved on ' + date.today().strftime("%B %d, %Y") + '.'

print('Output successful.')
hdf.close()
