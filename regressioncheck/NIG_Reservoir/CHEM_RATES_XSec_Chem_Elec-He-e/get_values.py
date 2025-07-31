import pandas as pd
import os
import re

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

directory = os.getcwd()

varName_01 = '003-CollRate001+002'
varName_02 = '010-ElecRelaxRate001+002-19.82'
varName_03 = '011-ElecRelaxRate001+002-20.61'
varName_04 = '012-Reaction001'
varName_05 = '004-CollRate001+003'
varName_06 = '009-BackscatterCollRate001+003'

for filename in sorted(os.listdir(directory)):
    if filename.startswith("PartAnalyze"):
        temp = get_numbers_from_filename(filename)
        df = pd.read_csv(filename)
        print(temp, *df[varName_01].tail(1).values,*df[varName_02].tail(1).values, *df[varName_03].tail(1).values, *df[varName_04].tail(1).values, *df[varName_05].tail(1).values,sep=",")