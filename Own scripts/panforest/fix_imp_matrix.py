#!/usr/bin/env python3
import sys
import pandas as pd

if len(sys.argv) != 3:
    sys.stderr.write(f"Usage: {sys.argv[0]} input_imp.csv output_imp_fixed.csv\n")
    sys.exit(1)

inp, outp = sys.argv[1], sys.argv[2]

# Read the importance matrix
df = pd.read_csv(inp, index_col=0)

# Strip trailing "_family_group" from both rows and columns
df.index = df.index.str.replace(r'_family_group$', '', regex=True)
df.columns = df.columns.str.replace(r'_family_group$', '', regex=True)

# Check that it's a square matrix
if df.shape[0] != df.shape[1]:
    sys.stderr.write(
        f"Warning: matrix is not square after cleanup "
        f"({df.shape[0]} rows vs {df.shape[1]} cols)\n"
    )

# Write back out
df.to_csv(outp)
print(f"Fixed file written to {outp}")
