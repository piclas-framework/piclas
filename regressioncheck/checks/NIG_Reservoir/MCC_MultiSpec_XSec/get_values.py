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
        print(temp,*df['005-CollRate001+004     '].tail(1).values,*df['008-CollRate002+004     '].tail(1).values,*df['010-CollRate003+004     '].tail(1).values, sep=",")