import pandas as pd


data = pd.read_csv('../data/dataset/dataset_merck_bde.csv')
data_bde = data[['zaretzki_index', 'bde', 'som']].copy()

data_bde['bde_rank'] = data_bde.groupby('zaretzki_index')['bde'].rank(ascending=True)
data_bde.to_csv('./dataset_merck_bde_rank.csv', index=False)

# Get the total number of som
som_count = data_bde[data_bde['som'] == 1].shape[0]

# Initialize counters
top_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

# Group by 'zaretzki_index' and iterate over each group
for _, group in data_bde.groupby('zaretzki_index'):
    for rank in range(1, 6):
        if any((group['som'] == 1) & (group['bde_rank'] == rank)):
            top_counts[rank] += 1
            break 

# Print the results
print("bde")
for rank in range(1, 6):
    print(f"Top {rank} count: {top_counts[rank]}")


data_ir = data[['zaretzki_index', 'ir', 'som']].copy()
data_ir['ir_rank'] = data_ir.groupby('zaretzki_index')['ir'].rank(ascending=False)
data_ir.to_csv('./dataset_merck_ir_rank.csv', index=False)

# Get the total number of som
som_count = data_ir[data_ir['som'] == 1].shape[0]

# Initialize counters
top_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0} 

for _, group in data_ir.groupby('zaretzki_index'):
    for rank in range(1, 6):
        if any((group['som'] == 1) & (group['ir_rank'] == rank)):
            top_counts[rank] += 1
            break 
print("intrinsic reactivity")
for rank in range(1, 6):
    print(f"Top {rank} count: {top_counts[rank]}")


