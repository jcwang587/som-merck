import pandas as pd

import os

data = pd.read_csv('../data/dataset/dataset_merck_bde.csv')
data_bde = data[['zaretzki_index', 'molecule_title', 'bde', 'som']].copy()

data_bde['bde_rank'] = data_bde.groupby('zaretzki_index')['bde'].rank(ascending=True)
data_bde.to_csv('./dataset_merck_bde_rank.csv', index=False)

# Get the count of unique zaretzki_index
unique_zaretzki_index = data_bde['zaretzki_index'].nunique()
print(f"number of unique zaretzki_index: {unique_zaretzki_index}")

# Initialize counters
top_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}

os.remove('bde_rank_molecules.txt')

# Group by 'zaretzki_index' and iterate over each group
for _, group in data_bde.groupby('zaretzki_index'):
    for rank in range(1, len(group)+1):
        if any((group['som'] == 1) & (group['bde_rank'] == rank)):
            if rank <= 10:
                top_counts[rank] += 1
            # print the molecule title and the rank
            with open('bde_rank_molecules.txt', 'a') as f:
                f.write(f"molecule title: {group['molecule_title'].values[0]}, rank: {rank}\n")
            break 


# Print the results
print("bde")
for rank in range(1, 11):
    print(f"Top {rank} count: {top_counts[rank]}")


data_ir = data[['zaretzki_index', 'molecule_title', 'ir', 'som']].copy()
data_ir['ir_rank'] = data_ir.groupby('zaretzki_index')['ir'].rank(ascending=False)
data_ir.to_csv('./dataset_merck_ir_rank.csv', index=False)

# Get the total number of som
som_count = data_ir[data_ir['som'] == 1].shape[0]

# Initialize counters
top_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0} 

for _, group in data_ir.groupby('zaretzki_index'):
    for rank in range(1, 11):
        if any((group['som'] == 1) & (group['ir_rank'] == rank)):
            top_counts[rank] += 1
            break 
print("intrinsic reactivity")
for rank in range(1, 11):
    print(f"Top {rank} count: {top_counts[rank]}")


data_sasa = data[['zaretzki_index', 'molecule_title', 'sasa_hydrogen_maestro', 'som']].copy()
data_sasa['sasa_rank'] = data_sasa.groupby('zaretzki_index')['sasa_hydrogen_maestro'].rank(ascending=True)
data_sasa.to_csv('./dataset_merck_sasa_rank.csv', index=False)

# Get the total number of som
som_count = data_sasa[data_sasa['som'] == 1].shape[0]

# Initialize counters
top_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0} 

for _, group in data_sasa.groupby('zaretzki_index'):
    for rank in range(1, 11):
        if any((group['som'] == 1) & (group['sasa_rank'] == rank)):
            top_counts[rank] += 1
            break

print("sasa")
for rank in range(1, 11):
    print(f"Top {rank} count: {top_counts[rank]}")

