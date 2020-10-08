import pandas as pd
import os
import re

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

directory = os.getcwd()

varName_01 = '003-Reaction001     '
varName_02 = '004-Reaction002     '
varName_03 = '005-Reaction003      '
varName_04 = '002-VibRelaxRate001+004'

for filename in sorted(os.listdir(directory)):
    if filename.startswith("PartAnalyze"):
        temp = get_numbers_from_filename(filename)
        df = pd.read_csv(filename)
        print(temp, *df[varName_01].tail(1).values,*df[varName_02].tail(1).values,*df[varName_03].tail(1).values,*df[varName_04].tail(1).values,sep=",")