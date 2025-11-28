#!/usr/bin/env python3

import csv
import os
import argparse
import pandas as pd
import shutil

def convert_gpa_to_tab(gpa_csv, output_tab, transposed=False):
    """
    Convert presence/absence CSV to Coinfinder .tab (Gene\tGenome).
    If transposed=True, assumes gpa_csv genomes are rows and genes are columns.
    """
    with open(gpa_csv, newline='') as inf, open(output_tab, 'w', newline='') as outf:
        reader = csv.reader(inf)
        header = next(reader)

        if not transposed:
            # Expect: row[0]=gene, header[1:]=genomes
            for row in reader:
                gene = row[0]
                for i in range(1, len(row)):
                    try:
                        if int(row[i]):
                            outf.write(f"{gene}\t{header[i]}\n")
                    except ValueError:
                        continue
        else:
            # Expect transposed: row[0]=genome, header[1:]=genes
            for row in reader:
                genome = row[0]
                for i in range(1, len(row)):
                    try:
                        if int(row[i]):
                            outf.write(f"{header[i]}\t{genome}\n")
                    except ValueError:
                        continue

# Parse input arguments
parser = argparse.ArgumentParser(description="Merge Dewar_data.csv with Barth_Baumdicker_Weigel_data.csv")
parser.add_argument("--dewar", required=True, help="Path to Dewar_data.csv")
parser.add_argument("--barth", required=True, help="Path to Barth_Baumdicker_Weigel_data.csv")
parser.add_argument("--gpa", required=True, help="Path to GPA_matrices_red/")
parser.add_argument("--tree", required=True, help="Path to Trees_red/")
parser.add_argument("--outdir", default=".", help="Path to save output files (default: current working directory)")
args = parser.parse_args()

# Merge Barth_Baumdicker_Weigel_data.csv with Dewar_data.csv
# Load files
df1 = pd.read_csv(args.dewar, keep_default_na=False, na_values=[])
df2 = pd.read_csv(args.barth, keep_default_na=False, na_values=[])

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

# Match and copy tree files
tree_directory = args.tree
tree_new_directory = os.path.join(args.outdir, "tree_matches_all")

# Extract species_taxids from merged_df
species_taxids = merged_df['species_taxid'].astype(str).tolist()

# List all .nwk files
nwk_files = [f for f in os.listdir(tree_directory) if f.endswith('.nwk')]

# Create the tree matches directory if it doesn't exist
os.makedirs(tree_new_directory, exist_ok=True)

# Track which taxids had a matching file
matched_taxids = set()

# Filter for matching files and copy them
for nwk_file in nwk_files:
    file_taxid = os.path.splitext(nwk_file)[0].split('_')[-1]
    if file_taxid in species_taxids:
        print(f"Matched tree file: {nwk_file} (taxid={file_taxid})")
        matched_taxids.add(file_taxid)
        source_path = os.path.join(tree_directory, nwk_file)
        destination_path = os.path.join(tree_new_directory, nwk_file)
        shutil.copy(source_path, destination_path)

print(f"Tree files for {len(matched_taxids)} matched taxa have been copied to '{tree_new_directory}'.")

# Match and copy GPA matrix files
gpa_directory = args.gpa
gpa_new_directory = os.path.join(args.outdir, "gpa_matches_all")

# Create the GPA matches directory if it doesn't exist
os.makedirs(gpa_new_directory, exist_ok=True)

# Loop over matched taxids and look for corresponding GPA files
num_taxa_by_taxid = {}
num_genes_by_taxid = {}

for taxid in matched_taxids:
    taxid_str = str(taxid)
    gpa_filename = f"{taxid_str}_REDUCED.csv"
    gpa_path = os.path.join(gpa_directory, gpa_filename)

    if os.path.exists(gpa_path):
        print(f"Processing {gpa_filename}...")
        
        # Read, transpose, extract num_taxa and num_genes, and save CSV and .tab
        df_gpa = pd.read_csv(gpa_path, keep_default_na=False, na_values=[])
        
        # Extract num_taxa and num_genes
        num_taxa = df_gpa.shape[0]
        num_genes = df_gpa.shape[1] - 1 

        num_taxa_by_taxid[taxid_str] = num_taxa
        num_genes_by_taxid[taxid_str] = num_genes
                
        print(f"Starting transpose for {gpa_filename}...")
        df_gpa_t = df_gpa.transpose()
        print(f"Finished transpose for {gpa_filename} (shape: {df_gpa_t.shape})")

        # Make original first row, the headers
        df_gpa_t.columns = df_gpa_t.iloc[0]
        df_gpa_t = df_gpa_t.drop(df_gpa_t.index[0])
        
        # Save transposed .csv
        dest_csv = os.path.join(gpa_new_directory, gpa_filename)
        df_gpa_t.to_csv(dest_csv, index=True)

        print(f"{gpa_filename}: num_taxa={num_taxa}, num_genes={num_genes}")

        # Save to .tab (no header, tab-separated)
        dest_tab = os.path.join(gpa_new_directory, f"{taxid}_REDUCED.tab")
        convert_gpa_to_tab(gpa_path, dest_tab, transposed = True)

print(
    f"GPA matrix files for {len(matched_taxids)} matched taxa "
    f"have been transposed, saved as CSV, and converted to .tab in '{gpa_new_directory}'."
)

# Keep only rows in merged_df that matched a file and drop unwanted columns
merged_df['species_taxid'] = merged_df['species_taxid'].astype(str)
merged_df['num_genes'] = merged_df['species_taxid'].map(num_genes_by_taxid)

merged_df = merged_df[merged_df['species_taxid'].isin(matched_taxids)]
merged_df = merged_df.drop(columns=['Number_of_Strain', 'genome_size', 'number_of_samples'], errors='ignore')

# Find the index of the species_taxid column
taxid_index = merged_df.columns.get_loc("species_taxid")

# Overwrite num_taxa and num_genes
merged_df['num_taxa'] = merged_df['species_taxid'].map(num_taxa_by_taxid)
merged_df['num_genes'] = merged_df['species_taxid'].map(num_genes_by_taxid)

taxid_index = merged_df.columns.get_loc("species_taxid")
cols = list(merged_df.columns)
for c in ['num_taxa','num_genes']:
    if c in cols:
        cols.remove(c)
        
# Insert after species_taxid
cols = cols[:taxid_index+1] + ['num_taxa','num_genes'] + cols[taxid_index+1:]
merged_df = merged_df[cols]

# Save final merged dataframe
merged_df.to_csv(os.path.join(args.outdir, "matched_all.csv"), index=False)

print(f"{len(merged_df)} matching rows kept in merged_df.")
print("Final merged dataframe saved as 'matched_all.csv'.")
