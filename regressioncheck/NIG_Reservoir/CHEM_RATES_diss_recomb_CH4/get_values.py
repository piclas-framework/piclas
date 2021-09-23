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
        print(temp,*df['002-Reaction001     '].tail(1).values,*df['003-Reaction002     '].tail(1).values,*df['004-Reaction003     '].tail(1).values,*df['005-Reaction004     '].tail(1).values,*df['006-Reaction005     '].tail(1).values,*df['007-Reaction006     '].tail(1).values,*df['008-Reaction007     '].tail(1).values,*df['009-Reaction008      '].tail(1).values, sep=",")
