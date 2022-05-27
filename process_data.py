import pandas as pd
from math import isnan

# Copy paste data from homepage, resolve all connected cells and run this script
df = pd.read_excel("raw_data.xlsx", sheet_name='Tabellenblatt1', engine="openpyxl")
column_names = df.columns
new_names = [x for x in df.iloc[0]]
df = df.drop(columns=column_names[0], axis=1)
df = df.drop([0])
df = df.rename(columns={old:new for old, new in zip(column_names, new_names)})
df = df.drop([i+1 for i,row in enumerate(df.iterrows()) if isnan(row[1][0])])
df = df.reset_index(drop=True)
df.to_excel("data.xlsx")