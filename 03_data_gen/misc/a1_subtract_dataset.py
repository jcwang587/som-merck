import pandas as pd
import os


raw_data = pd.read_csv("../data/dataset/dataset_all_sites.csv")

merck_dir = "../data/merck"
merck_files = os.listdir(merck_dir)

# remove the first 4 characters from each file name and change to lowercase
merck_files = [file[4:].replace(".mae", "") for file in merck_files]
print(merck_files)

# leave only the values in the molecule_title column that are in merck_files
raw_data = raw_data[raw_data["molecule_title"].isin(merck_files)]

# save the result to a new csv file
raw_data.to_csv("../data/dataset/dataset_merck_sub.csv", index=False)

# get the unique values in the molecule_title column
unique_molecules = raw_data["molecule_title"].unique()
print(len(unique_molecules))
print(unique_molecules)

# get the name in the unique_molecules list that is not in merck_files
for molecule in merck_files:
    if molecule not in unique_molecules:
        print(molecule)
