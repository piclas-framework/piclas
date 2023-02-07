import pandas as pd
import os
import re

def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)

directory = os.getcwd()

# H2 + H3Ion reactions
varName_01 = '005-Reaction001'
varName_02 = '006-Reaction002'
varName_03 = '007-Reaction003'
print('===============================================')
print('Reactions: H2 + H3Ion')
for filename in sorted(os.listdir(directory)):
    if filename.startswith("PartAnalyze"):
        temp = get_numbers_from_filename(filename)
        df = pd.read_csv(filename)
        print(temp, *df[varName_01].tail(1).values,*df[varName_02].tail(1).values,*df[varName_03].tail(1).values,sep=",")
