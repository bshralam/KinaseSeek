import pandas as pd

# Path to raw ChEMBL export
input_path = "../data/chembl_export.csv"
output_path = "../data/chembl_export_clean.csv"

# Read CSV safely, skipping any problematic rows
df = pd.read_csv(input_path, sep=";", on_bad_lines="skip")

# Keep only the relevant columns
df_filtered = df[['Name', 'Smiles']].dropna()

# Rename columns to lowercase as expected by RDKit script
df_filtered.columns = ['name', 'smiles']

# Save the cleaned CSV
df_filtered.to_csv(output_path, index=False)

print(f"Clean file saved to: {output_path}")

