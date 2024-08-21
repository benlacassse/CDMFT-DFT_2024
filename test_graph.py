import pandas as pd
import numpy as np

df = pd.read_csv('CCOC_Dop.csv', index_col = 'Couche')

print(df.loc["3l_C"].tolist())