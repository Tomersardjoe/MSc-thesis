import csv
import os
import pandas as pd
import shutil

def convert_gpa_to_tab(gpa_csv, output_tab):
    """
    Convert a presence/absence CSV matrix into Coinfinder's 2-column .tab format:
        Gene    Genome
    Only presence calls (non-zero) are written.
    """
    with open(gpa_csv, newline='') as inf, open(output_tab, 'w', newline='') as outf:
        reader = csv.reader(inf)  # defaults to comma delimiter
        header = next(reader)
        
        for row in reader:
            gene = row[0]
            for i in range(1, len(row)):
                try:
                    if int(row[i]):  # 1 means present
                        outf.write(f"{gene}\t{header[i]}\n")
                except ValueError:
                    # Ignore non-integer cells if any
                    continue

# === Part 1: Merge Franz_data.csv with pangenome_lifestyles.csv ===
# Load files
df1 = pd.read_csv("pangenome_lifestyles.csv", keep_default_na=False, na_values=[])
df2 = pd.read_csv("Franz_data.csv", keep_default_na=False, na_values=[])

# Normalize species names: strip, lowercase, replace spaces with underscores
df1['Species'] = df1['Species'].str.strip().str.lower().str.replace(' ', '_')
df2['species_name'] = df2['species_name'].str.strip().str.lower().str.replace(' ', '_')

# Merge df1 into df2
merged_df = df2.merge(
    df1,
    left_on='species_name',
    right_on='Species',
    how='left'
)

# Drop helper column and rows with empty Species
merged_df = merged_df.drop(columns=['species_name'])
merged_df = merged_df.dropna(subset=['Species'])
merged_df = merged_df[merged_df['Species'].str.strip() != '']

# === Part 2: Match and copy tree files ===
directory = "Trees_red"
tree_new_directory = "tree_matches_all"

# Extract species_taxids from merged_df
species_taxids = merged_df['species_taxid'].astype(str).tolist()

# List all .nwk files
nwk_files = [f for f in os.listdir(directory) if f.endswith('.nwk')]

# Create the tree matches directory if it doesn't exist
os.makedirs(tree_new_directory, exist_ok=True)

# Track which taxids had a matching file
matched_taxids = set()

# Filter for matching files and copy them
for nwk_file in nwk_files:
    file_taxid = os.path.splitext(nwk_file)[0].split('_')[-1]
    if file_taxid in species_taxids:
        matched_taxids.add(file_taxid)
        source_path = os.path.join(directory, nwk_file)
        destination_path = os.path.join(tree_new_directory, nwk_file)
        shutil.copy(source_path, destination_path)

print(f"Tree files for {len(matched_taxids)} matched taxa have been copied to '{tree_new_directory}'.")

# === Part 3: Match and copy GPA matrix files ===
gpa_directory = "GPA_matrices_red"
gpa_new_directory = "gpa_matches_all"

# Create the GPA matches directory if it doesn't exist
os.makedirs(gpa_new_directory, exist_ok=True)

# Loop over matched taxids and look for corresponding GPA files
for taxid in matched_taxids:
    gpa_filename = f"{taxid}_REDUCED.csv"
    gpa_path = os.path.join(gpa_directory, gpa_filename)

    if os.path.exists(gpa_path):
        print(f"Processing {gpa_filename}...")
        
        # Read, transpose, and save CSV
        df_gpa = pd.read_csv(gpa_path, keep_default_na=False, na_values=[])
        
        print(f"Starting transpose for {gpa_filename}...")
        df_gpa_t = df_gpa.transpose()
        print(f"Finished transpose for {gpa_filename} (shape: {df_gpa_t.shape})")

        # Make original first row, the headers
        df_gpa_t.columns = df_gpa_t.iloc[0]
        df_gpa_t = df_gpa_t.drop(df_gpa_t.index[0])
        
        # Save transposed .csv
        dest_csv = os.path.join(gpa_new_directory, gpa_filename)
        df_gpa_t.to_csv(dest_csv, index=True)

        # Vectorised .tab creation
        stacked = df_gpa_t.stack().reset_index()
        stacked.columns = ['Gene', 'Genome', 'Value']

        # Keep only presence calls (non-zero)
        stacked = stacked[stacked['Value'].astype(str) != '0']

        # Save to .tab (no header, tab-separated)
        dest_tab = os.path.join(gpa_new_directory, f"{taxid}_REDUCED.tab")
        stacked[['Gene', 'Genome']].to_csv(dest_tab, sep='\t', index=False, header=False)

print(
    f"GPA matrix files for {len(matched_taxids)} matched taxa "
    f"have been transposed, saved as CSV, and converted to .tab in '{gpa_new_directory}'."
)

# Keep only rows in merged_df that matched a file
merged_df = merged_df[merged_df['species_taxid'].astype(str).isin(matched_taxids)]

# === Save final merged dataframe ===
merged_df.to_csv("matched_all.csv", index=False)

print(f"{len(merged_df)} matching rows kept in merged_df.")
print("Final merged dataframe saved as 'matched_all.csv'.")
