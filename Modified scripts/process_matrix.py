#!/usr/bin/env python3

import argparse
import os
import sys
import math
import pandas as pd
import rf_module as rf

def get_args():
    """Parse and validate command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", dest="infile", type=str)
    parser.add_argument("-r", "--roary", dest="roary", action="store_true")
    parser.add_argument("-o", "--output", dest="outfile", type=str)
    parser.add_argument("-d", "--directory", dest="outdir", type=str, default="process_matrix_outfiles")
    args = parser.parse_args()
    if None in [args.infile, args.outfile]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.infile, args.outfile, args.roary, args.outdir]

def safe_multiindex_from_tuples(tuples, original_names):
    """Create a safe MultiIndex with adjusted level names."""
    level_count = len(tuples[0])
    names = original_names[:level_count] + [''] * (level_count - len(original_names))
    return pd.MultiIndex.from_tuples(tuples, names=names)

def write_gene_lists(matrix, outdir):
    """Identify and write constant genes, core genes, and singletons to text files."""
    constant = matrix[matrix.sum(axis=1) == matrix.shape[1]].index
    core = matrix[matrix.sum(axis=1) >= math.ceil(0.05 * matrix.shape[1])].index
    singletons = matrix[matrix.sum(axis=1) == 1].index

    with open(os.path.join(outdir, "constant_genes.txt"), "w") as out:
        out.write("\n".join(map(str, constant)))
    with open(os.path.join(outdir, "core_genes.txt"), "w") as out:
        out.write("\n".join(map(str, core)))
    with open(os.path.join(outdir, "singletons.txt"), "w") as out:
        out.write("\n".join(map(str, singletons)))

    matrix = matrix.drop(index=constant.union(singletons))
    return matrix

def collapse_genes(matrix, outdir, simulated=False):
    """Collapse genes with identical presence/absence profiles into families."""
    collapsed_rows = pd.DataFrame(columns=matrix.columns)
    identical_sets = {}
    pattern = None
    indices = []

    non_unique = matrix[matrix.duplicated(keep=False)].sort_values(by=list(matrix.columns))

    for i, row in non_unique.iterrows():
        current = list(row)
        if current == pattern:
            identical_sets[group].append(i)
        else:
            pattern = current
            group = f"family_group_{len(identical_sets)+1}"
            collapsed_rows.loc[len(identical_sets)] = current
            indices.append((group, "", "family_group"))
            identical_sets[group] = [i]

    if indices:
        print(f"Found {len(indices)} duplicate gene entries.")
        collapsed_rows.index = safe_multiindex_from_tuples(indices, matrix.index.names)
    else:
        print("No duplicate genes found. Skipping collapsing step.")
        collapsed_rows.index = pd.MultiIndex.from_tuples([], names=matrix.index.names)


    with open(os.path.join(outdir, "non-unique_genes.csv"), "w") as out:
        for key, value in identical_sets.items():
            out.write(f"{key}\t{','.join(map(str, value))}\n")

    matrix = matrix.drop_duplicates(keep=False)
    return pd.concat([matrix, collapsed_rows])

def collapse_genomes(matrix, outdir):
    """Collapse genomes with identical gene family profiles into unified groups."""
    matrix = matrix.transpose()
    collapsed = pd.DataFrame(columns=matrix.columns)
    identical_sets = {}
    pattern = None
    indices = []

    non_unique = matrix[matrix.duplicated(keep=False)].sort_values(by=list(matrix.columns))

    for i, row in non_unique.iterrows():
        current = list(row)
        if current == pattern:
            identical_sets[group].append(i)
        else:
            pattern = current
            group = f"genome_group_{len(identical_sets)+1}"
            collapsed.loc[len(identical_sets)] = current
            indices.append(group)
            identical_sets[group] = [i]

    collapsed.index = indices

    with open(os.path.join(outdir, "non-unique_genomes.csv"), "w") as out:
        for key, value in identical_sets.items():
            out.write(f"{key}\t{','.join(map(str, value))}\n")

    matrix = matrix.drop_duplicates(keep=False)
    return pd.concat([matrix, collapsed]).transpose()

def convert_roary(roary_matrix):
    """Format Roary input to conform with standard metadata structure."""
    index = roary_matrix.index
    roary_matrix.index = pd.MultiIndex.from_tuples(
        [(i[0], i[1], i[2].replace(",", ";")) for i in index],
        names=index.names
    )
    roary_matrix.replace(",", ";", regex=True, inplace=True)
    return roary_matrix.drop(roary_matrix.columns[:11], axis=1)

def write_flattened_matrix(matrix, outfile, label_name="Gene"):
    """Export matrix with unpacked metadata for simulated input."""
    meta_cols = [label_name, "", " "]
    genome_cols = [c for c in matrix.columns if c not in meta_cols]
    matrix = matrix[meta_cols + genome_cols]
    matrix.index.name = None
    matrix.reset_index(drop=True, inplace=True)
    matrix.to_csv(outfile, index=False, sep=",", doublequote=False, quoting=0)

def main():
    """Run complete preprocessing pipeline."""
    infile, outfile, roary, outdir = get_args()
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, os.path.basename(outfile))

    matrix = pd.read_csv(infile, dtype=str)
    simulated = False

    if roary:
        meta = ["Gene", "Non-unique Gene name", "Annotation"]
        matrix.index = pd.MultiIndex.from_frame(matrix[meta])
        matrix = matrix.drop(columns=meta)
        matrix = convert_roary(matrix)
    elif set(matrix.columns[:3]) == {"Gene", "Non-unique Gene name", "Annotation"}:
        meta = ["Gene", "Non-unique Gene name", "Annotation"]
        matrix.index = pd.MultiIndex.from_frame(matrix[meta])
        matrix = matrix.drop(columns=meta)
    else:
        matrix.index = pd.Index(matrix.iloc[:, 0].astype(str), name="Gene")
        matrix = matrix.iloc[:, 1:]
        simulated = True

    matrix = rf.preprocess_df(matrix, 0, 0, 0)

    print("Writing singletons, core and constant genes")
    matrix = write_gene_lists(matrix, outdir)

    print("Collapsing identical genes and genomes")
    matrix = collapse_genes(matrix, outdir, simulated)
    matrix = collapse_genomes(matrix, outdir)
    
    # Stratification compatibility
    matrix = matrix[matrix.sum(axis=1) > 1] # Keep genes present in at least 2 genomes
    matrix = matrix[(matrix.shape[1] - matrix.sum(axis=1)) >= 2] # Keep genes absent in at least 2 genomes


    print("Writing collapsed matrix")
    if simulated:
        matrix["Gene"] = [i[0] if isinstance(i, tuple) else i for i in matrix.index]
        matrix[""] = [i[1] if isinstance(i, tuple) and len(i) > 1 else "" for i in matrix.index]
        matrix[" "] = [i[2] if isinstance(i, tuple) and len(i) > 2 else "" for i in matrix.index]
        write_flattened_matrix(matrix, outfile)
    else:
        matrix["Gene"] = [i[0] for i in matrix.index]
        matrix["Non-unique Gene name"] = [i[1] for i in matrix.index]
        matrix["Annotation"] = [i[2] for i in matrix.index]
        matrix.reset_index(drop=True, inplace=True)
        meta_cols = ["Gene", "Non-unique Gene name", "Annotation"]
        genome_cols = [c for c in matrix.columns if c not in meta_cols]
        matrix = matrix[meta_cols + genome_cols]
        matrix.to_csv(outfile, sep=",", index=False, doublequote=False, quoting=0)

if __name__ == "__main__":
    main()
