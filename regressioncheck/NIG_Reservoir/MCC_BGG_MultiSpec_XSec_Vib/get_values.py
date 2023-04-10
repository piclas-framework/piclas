import pandas as pd
import os
import re

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

directory = os.getcwd()

varName_01 = '009-VibRelaxRate001+004'
varName_02 = '010-VibRelaxRate002+004 '

for filename in sorted(os.listdir(directory)):
    if filename.startswith("Database"):
        temp = get_numbers_from_filename(filename)
        df = pd.read_csv(filename)
        print(temp,*df[varName_01].tail(1).values,*df[varName_02].tail(1).values,sep=",")