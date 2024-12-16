import pandas as pd

# List of file names
file_list = [
    "max_fixed10_purgefalse_thread69.csv",
    "max_fixed10_purgetrue_thread59.csv",
    "max_fixed20_purgefalse_thread35.csv",
    "max_fixed20_purgetrue_thread54.csv",
    "max_fixed30_purgefalse_thread79.csv",
    "max_fixed30_purgetrue_thread53.csv",
    "max_fixed40_purgefalse_thread78.csv",
    "max_fixed40_purgetrue_thread29.csv"
]

Purge = [False, True] * 4 
FixedValues = [10, 10, 20, 20, 30, 30, 40, 40]  

target_index = 2000
extracted_data = []

for purge, num_values, file_name in zip(Purge, FixedValues, file_list):
    data = pd.read_csv(file_name)
    
    filtered_data = data[data["Iter"] == target_index]
    
    filtered_data["Purge"] = purge
    filtered_data["FixedValues"] = num_values
    
    extracted_data.append(filtered_data[["BestFit", "Purge", "FixedValues"]])

result_df = pd.concat(extracted_data, ignore_index=True)
result_df.to_csv("max.csv", index=False)

print("Processing complete! Combined file saved as 'output_combined.csv'.")
