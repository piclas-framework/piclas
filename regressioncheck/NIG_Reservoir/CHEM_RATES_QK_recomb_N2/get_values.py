import pandas as pd
import os
import re

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

directory = os.getcwd()

for filename in sorted(os.listdir(directory)):
    if filename.startswith("Database"):
        temp = get_numbers_from_filename(filename)
        df = pd.read_csv(filename)
        print(temp,*df['007-Reaction003     '].tail(1).values,*df['008-Reaction004      '].tail(1).values, sep=",")
