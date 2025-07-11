#!/usr/bin/env python3
"""Process the matrix so that it's ready for analysis."""

import argparse
import re
import sys
import math
import pandas as pd
import rf_module as rf


def get_args():
    """Get user arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", dest="infile",
                        type=str, help="input file")
    parser.add_argument("-r", "--roary", dest="roary",
                        help="toggle for roary input", action="store_true")
    parser.add_argument("-o", "--output", dest="outfile",
                        type=str, help="Output matrix file")
    args = parser.parse_args()
    if None in [args.infile, args.outfile]:
        parser.print_help(sys.stderr)
        sys.exit(0)
    return [args.infile, args.outfile, args.roary]


# MODIFICATION [Detect whether input is simulated or real and build 3-level index]
def ensure_index_format(matrix):
    """Ensure matrix has a 3-level index: Gene, Non-unique Gene name, Annotation."""
    if matrix.index.nlevels == 1:
        matrix.index = pd.MultiIndex.from_arrays(
            [matrix.index, [''] * len(matrix), [''] * len(matrix)],
            names=['Gene', 'Non-unique Gene name', 'Annotation']
        )
        return matrix, True  # Simulated input
    return matrix, False


# MODIFICATION [Safely create MultiIndex with correct level names]
def safe_multiindex_from_tuples(tuples, original_names):
    """Safely create MultiIndex with matching name count."""
    level_count = len(tuples[0])
    if len(original_names) < level_count:
        names = original_names + [''] * (level_count - len(original_names))
    else:
        names = original_names[:level_count]
    return pd.MultiIndex.from_tuples(tuples, names=names)


def write_gene_lists(matrix):
    """Write and filter constant, core, and singleton gene families."""
    constant = list(matrix[matrix.sum(axis=1) == matrix.shape[1]].index)
    threshold = math.ceil(0.05 * matrix.shape[1])
    core = list(matrix[matrix.sum(axis=1) >= threshold].index)
    singletons = list(matrix[matrix.sum(axis=1) == 1].index)

    with open("constant_genes.txt", "w", encoding="utf-8") as out:
        out.write("\n".join(map(str, constant)))
    with open("core_genes.txt", "w", encoding="utf-8") as out:
        out.write("\n".join(map(str, core)))
    with open("singletons.txt", "w", encoding="utf-8") as out:
        out.write("\n".join(map(str, singletons)))

    matrix = matrix.drop(index=constant)
    matrix = matrix.drop(index=singletons)
    return matrix


# MODIFICATION [Collapse gene families using consistent MultiIndex format]
def collapse_genes(matrix, simulated=False):
    """Collapse genes with identical presence/absence patterns."""
    collapsed_rows = pd.DataFrame(columns=matrix.columns)
    identical_sets = {}
    count = 0
    pattern = 0
    indices = []

    non_unique = matrix[matrix.duplicated(keep=False)]
    non_unique = non_unique.sort_values(by=list(non_unique.columns))

    for i in range(len(non_unique.index)):
        if list(non_unique.iloc[i]) == pattern:
            identical_sets[f"family_group_{count}"].append(non_unique.iloc[i].name)
        else:
            collapsed_rows.loc[count] = list(non_unique.iloc[i])
            pattern = list(non_unique.iloc[i])
            count += 1
            index_tuple = (f"family_group_{count}", "", "family_group")
            indices.append(index_tuple)
            identical_sets[f"family_group_{count}"] = [non_unique.iloc[i].name]

    # MODIFICATION [Use safe_multiindex_from_tuples to avoid index errors]
    collapsed_rows.index = safe_multiindex_from_tuples(indices, matrix.index.names)

    with open("non-unique_genes.csv", "w", encoding="utf-8") as out:
        for key, value in identical_sets.items():
            out.write(key + "\t" + ",".join(map(str, value)) + "\n")

    matrix = matrix.drop_duplicates(keep=False)
    matrix = pd.concat([matrix, collapsed_rows])
    return matrix


def collapse_genomes(matrix):
    """Collapse genomes with identical presence absence patterns."""
    matrix = matrix.transpose()
    collapsed_genomes = pd.DataFrame(columns=matrix.columns)
    identical_sets = {}
    count = 0
    pattern = 0
    indices = []

    non_unique = matrix[matrix.duplicated(keep=False)]
    non_unique = non_unique.sort_values(by=list(non_unique.columns))

    for i in range(len(non_unique.index)):
        if list(non_unique.iloc[i]) == pattern:
            identical_sets[f"genome_group_{count}"].append(non_unique.iloc[i].name)
        else:
            pattern = list(non_unique.iloc[i])
            collapsed_genomes.loc[count] = list(non_unique.iloc[i])
            count += 1
            indices.append(f"genome_group_{count}")
            identical_sets[f"genome_group_{count}"] = [non_unique.iloc[i].name]

    collapsed_genomes.index = indices

    with open("non-unique_genomes.csv", "w", encoding="utf-8") as out:
        for key, value in identical_sets.items():
            out.write(key + "\t" + ",".join(map(str, value)) + "\n")

    matrix = matrix.drop_duplicates(keep=False)
    matrix = pd.concat([matrix, collapsed_genomes])
    matrix = matrix.transpose()
    return matrix


def convert_roary(roary_matrix):
    """Convert Roary matrix into Panaroo format."""
    index = roary_matrix.index
    names = index.names
    new_index = [(ind[0], ind[1], ind[2].replace(",", ";")) for ind in index]
    roary_matrix.index = pd.MultiIndex.from_tuples(new_index, names=names)
    roary_matrix = roary_matrix.replace(",", ";", regex=True)
    roary_matrix.drop(roary_matrix.columns[range(11)], axis=1, inplace=True)
    return roary_matrix

# MODIFICATION [Flatten nested tuple to access index elements]
def write_flattened_matrix(matrix, outfile, label_name="Gene"):
    """Export flattened matrix for simulated data with proper metadata structure."""

    matrix = matrix.copy()

    # Reorder columns
    meta_cols = [label_name, "", " "]
    genome_cols = [col for col in matrix.columns if col not in meta_cols]
    matrix = matrix[meta_cols + genome_cols]
    print("🧬 Final export columns:", matrix.columns.tolist())

    # Export
    matrix.index.name = None
    matrix.reset_index(drop=True, inplace=True)
    matrix.to_csv(outfile, index=False, sep=",", doublequote=False, quoting=0)

# MODIFICATION [Handle real vs simulated input and apply complete pipeline]
def main():
    """Main preprocessing pipeline."""
    infile, outfile, roary = get_args()
    matrix = pd.read_csv(infile, header=0, dtype=str, index_col=None)
    print("\n📥 Initial input matrix shape:", matrix.shape)
    print("📄 Initial column headers:", matrix.columns.tolist())
    print("🔎 First few rows:\n", matrix.head())

    # Detect input data type and assign metadata columns accordingly
    if roary:
        # Roary input has extra metadata columns upfront
        metadata_cols = ["Gene", "Non-unique Gene name", "Annotation"]
        matrix.columns = metadata_cols + list(matrix.columns[len(metadata_cols):])
        matrix.index = pd.MultiIndex.from_frame(matrix[metadata_cols])
        matrix = matrix.drop(columns=metadata_cols)
        matrix = convert_roary(matrix)
        simulated = False

    elif matrix.shape[1] >= 3 and set(matrix.columns[:3]) == {"Gene", "Non-unique Gene name", "Annotation"}:
        # Panaroo-style input with exactly 3 metadata columns
        metadata_cols = ["Gene", "Non-unique Gene name", "Annotation"]
        matrix.index = pd.MultiIndex.from_frame(matrix[metadata_cols])
        matrix = matrix.drop(columns=metadata_cols)
        simulated = False

    else:
        # Simulated input with only Gene column
        matrix.index = pd.Index(matrix.iloc[:, 0].astype(str), name="Gene")
        matrix = matrix.iloc[:, 1:]
        simulated = True

    print("\n✅ Matrix index after format detection:")
    print("🔢 Index type:", type(matrix.index))
    print("🧬 Index levels:", matrix.index.nlevels)
    print("🧾 Index names:", matrix.index.names)
    print("👁️ First few index entries:", matrix.index[:5])

    matrix = rf.preprocess_df(matrix, 0, 0, 0)

    print("\n🧪 After preprocess_df:")
    print("🔢 Index type:", type(matrix.index))
    print("🧬 Index levels:", matrix.index.nlevels)
    print("🧾 Index names:", matrix.index.names)
    print("📐 Matrix shape:", matrix.shape)


    print("Writing singletons, core and constant genes")
    matrix = write_gene_lists(matrix)

    print("Collapsing identical genes and genomes")
    matrix = collapse_genes(matrix, simulated)
    matrix = collapse_genomes(matrix)

    # MODIFICATION [Avoid singleton families created by collapsing]
    matrix = matrix[matrix.sum(axis=1) > 1]

    # MODIFICATION [Conditional output formatting]
    print("Writing collapsed matrix")
    if simulated:
        # Unpack index tuples into columns before writing
        matrix = matrix.copy()
        matrix["Gene"] = [idx[0] if isinstance(idx, tuple) else idx for idx in matrix.index]
        matrix[""] = [idx[1] if isinstance(idx, tuple) and len(idx) > 1 else "" for idx in matrix.index]
        matrix[" "] = [idx[2] if isinstance(idx, tuple) and len(idx) > 2 else "" for idx in matrix.index]

        matrix.index.name = None
        matrix.reset_index(drop=True, inplace=True)
        print("\n🧾 Columns before export (simulated):", matrix.columns.tolist())
        print("👁️ First few rows:\n", matrix.head())
        print("🔢 Index type:", type(matrix.index))
        print("🧬 Index levels:", matrix.index.nlevels)
        print("🧾 Index names:", matrix.index.names)

        write_flattened_matrix(matrix, outfile, label_name="Gene")
    else:
        matrix = matrix.copy()
        matrix["Gene"] = [idx[0] for idx in matrix.index]
        matrix["Non-unique Gene name"] = [idx[1] for idx in matrix.index]
        matrix["Annotation"] = [idx[2] for idx in matrix.index]
        matrix.reset_index(drop=True, inplace=True)

        # Reoder columns
        meta_cols = ["Gene", "Non-unique Gene name", "Annotation"]
        genome_cols = [col for col in matrix.columns if col not in meta_cols]
        matrix = matrix[meta_cols + genome_cols]

        matrix.to_csv(outfile, sep=",", index=False, doublequote=False, quoting=0)

if __name__ == "__main__":
    main()
