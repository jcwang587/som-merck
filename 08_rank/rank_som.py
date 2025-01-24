import pandas as pd


data = pd.read_csv('../data/dataset/dataset_merck_bde.csv')
data_bde = data[['zaretzki_index', 'bde', 'som']].copy()

data_bde['bde_rank'] = data_bde.groupby('zaretzki_index')['bde'].rank(ascending=True)
data_bde.to_csv('./dataset_merck_bde_rank.csv', index=False)

# Get the total number of som
som_count = data_bde[data_bde['som'] == 1].shape[0]

# Count the number of data points where som is 1 and bde_rank is 1, 2, or 3
top_1_count = data_bde[(data_bde['som'] == 1) & (data_bde['bde_rank'] == 1)].shape[0]
top_2_count = data_bde[(data_bde['som'] == 1) & (data_bde['bde_rank'] == 2)].shape[0]
top_3_count = data_bde[(data_bde['som'] == 1) & (data_bde['bde_rank'] == 3)].shape[0]
top_4_count = data_bde[(data_bde['som'] == 1) & (data_bde['bde_rank'] == 4)].shape[0]
top_5_count = data_bde[(data_bde['som'] == 1) & (data_bde['bde_rank'] == 5)].shape[0]

print(f"som count: {som_count}")
print(f"Top 1 count: {top_1_count}")
print(f"Top 2 count: {top_2_count}")
print(f"Top 3 count: {top_3_count}")
print(f"Top 4 count: {top_4_count}")
print(f"Top 5 count: {top_5_count}")

