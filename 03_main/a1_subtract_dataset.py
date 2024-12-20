import pandas as pd
import os


raw_data = pd.read_csv("../data/dataset/dataset_all_sites.csv")

merck_dir = "../data/merck"
merck_files = os.listdir(merck_dir)
merck_files = [file.split("_")[-1].replace(".mae", "") for file in merck_files]

# leave only the values in the molecule_title column that are in merck_files
raw_data = raw_data[raw_data["molecule_title"].isin(merck_files)]

# save the result to a new csv file
raw_data.to_csv("../data/dataset/dataset_merck.csv", index=False)

# get the unique values in the molecule_title column
unique_molecules = raw_data["molecule_title"].unique()
print(len(unique_molecules))
