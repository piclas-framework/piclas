import requests
import h5py
import pandas as pd
import io
import re
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
# Base URL of the query
URL_base = 'https://physics.nist.gov/cgi-bin/ASD/energy1.pl'
# Output of the HDF5 database
hdf = h5py.File('Electronic_State_Database.h5', 'w')
# Species list to be included in the output database
species_list = ['Xe','XeIon1','XeIon2']
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
  URL_request = URL_base + '?de=0&spectrum=' + current_species_NIST + '&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&level_out=on&j_out=on&temp=&submit=Retrieve+Data'
  response = requests.get(URL_request)
  # Check the website
  if response.status_code == 200:
    print('Download successful.')
  # Get the data from the NIST website
  data = response.content.decode('utf-8')
  # Clean the data
  data = data.replace('"','')
  data = data.replace('= ','')
  data = data.replace('=','')
  data = data.replace(' ','')
  data = data.replace('[','')
  data = data.replace(']','')
  # Read-in data as a pandas dataframe
  data = pd.read_csv(io.StringIO(data),skipinitialspace=True,delimiter=",",usecols=['J','Level(cm-1)'], na_values=['---'])
  data.dropna(inplace=True)
  # Convert fraction entries to floats (5/2 -> 2.5)
  if data['J'].dtype == 'object':
    data['J'] = data.stack().map(lambda x: pd.eval(x)).unstack()
  # Convert total angeluar momentum to degeneracy
  data['J'] = 2 * data['J'] + 1
  # Convert 1/cm to K
  data['Level(cm-1)'] = 100 * data['Level(cm-1)'] * 1.986E-025 / 1.38065E-023           # 1/cm * 100cm/m * (J m) / (J/K) = K
  data.columns=['Degeneracy', 'Level(K)']
  # Write to hdf
  hdf.create_dataset(current_species, data=data)
print('Output successful.')
hdf.close()