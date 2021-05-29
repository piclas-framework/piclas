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
species_list = ['N','O','NIon1','OIon1','H','Xe','XeIon1','XeIon2','Ar','ArIon1','ArIon2']
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
  URL_request = URL_base + '?de=0&spectrum=' + current_species_NIST + '&units=0&format=2&output=0&page_size=15&multiplet_ordered=1&level_out=on&j_out=on&temp=&submit=Retrieve+Data'
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
  data = pd.read_csv(io.StringIO(data),skipinitialspace=True,delimiter=",",usecols=['J','Levelcm-1'], na_values=['---'])
  # Drop rows with an empty cell and '---' entries (which were converted to N/A during read-in)
  # and also drop all rows following the first '---' entry
  J     = data.columns.get_loc("J")
  Level = data.columns.get_loc("Levelcm-1")
  rows_with_nan = [index for index, row in data.iterrows() if row.isnull().any()]
  max_level = rows_with_nan[0]
  drop_to_end = 1 - (len(data)-max_level)
  if drop_to_end < 0:
    data = data.iloc[:drop_to_end]
  elif drop_to_end == 0:
    pass
  else:
    print("ERROR: drop_to_end must be negative!")
    exit(1)
  if data['J'].dtype == 'float64':
      data.iloc[max_level,J] = 0.0
  else:
      data.iloc[max_level,J] = "0.0"

  # Check the datatype: if its a float, then all non-numerical characters have already been removed
  if data['Levelcm-1'].dtype != 'float64':
    # Drop rows with a question mark ("This level/line may not be real.")
    data.drop(data[data['Levelcm-1'].str.contains(r'[?]')].index,inplace=True)
    # Drop rows with a +x ("The relative positions of the levels within such a system are accurate within experimental uncertainties, but no experimental connection between this system and the other levels of the spectrum has been made.")
    data.drop(data[data['Levelcm-1'].str.contains(r'[+x]')].index,inplace=True)

  # Execute fractions and convert J to g
  if data['J'].dtype != 'float64':
      for i in range(len(data['J'])):
        # Convert fraction entries to floats (5/2 -> 2.5)
        #print(type(data.iloc[i,J]))
        found = re.search(r'\d+/\d+', data.iloc[i,J])
        if found:
          numbers = data.iloc[i,J].split("/")
          data.iloc[i,J] = float(numbers[0])/float(numbers[1])
        else:
          data.iloc[i,J] = float(data.iloc[i,J])

  # Convert J to g
  data.iloc[:,J] = 2 * data.iloc[:,J] + 1.0

  # Convert type for HDF5 output
  data['J'] = data['J'].astype(float)

  # Convert 1/cm to K
  data['Levelcm-1'] = data['Levelcm-1'].astype(float)
  data['Levelcm-1'] = 100 * data['Levelcm-1'] * 1.986E-025 / 1.38065E-023           # 1/cm * 100cm/m * (J m) / (J/K) = K

  x_old=0.
  for i in range(len(data['Levelcm-1'])):
    x = data.iloc[i,Level]
    if x < x_old:
      print('Error in level %s: the energy is not increasing with the levels E2=%s < E1=%s' % (i,x,x_old))
      exit(1)
    else:
      x_old = x

  # Write to hdf: If data set already exists, delete the old set first
  data.columns=['Degeneracy', 'Level(K)']
  if current_species in hdf.keys():
    del hdf[current_species]
    print('Old dataset replaced.')
  dataset = hdf.create_dataset(current_species, data=data)
  # Write attributes for source and time of retrieval
  dataset.attrs['Reference'] = 'Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2020). NIST Atomic Spectra Database (ver. 5.8), [Online]. Available: https://physics.nist.gov/asd. National Institute of Standards and Technology, Gaithersburg, MD. DOI: https://doi.org/10.18434/T4W30F. Retrieved on ' + date.today().strftime("%B %d, %Y") + '.'

print('Output successful.')
hdf.close()
