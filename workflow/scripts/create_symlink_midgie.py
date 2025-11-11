import os
import sys
import pandas as pd

# Usage: python script.py output_dir df1.csv df2.csv
outdir, df1_path, df2_path = sys.argv[1], sys.argv[2], sys.argv[3]

df1 = pd.read_csv(df1_path)
df2 = pd.read_csv(df2_path)

# Extract identifier before first space in df1["name"]
df1["ident"] = df1["name"].str.split().str[0]

# Merge on ident
merged = df1.merge(df2, on="ident", how="inner")

# Create output dir if missing
os.makedirs(outdir, exist_ok=True)

# Symlink each genome file
for _, row in merged.iterrows():
    src = row["genome_filename"]
    dst = os.path.join(outdir, os.path.basename(src))
    if not os.path.exists(dst):
        os.symlink(src, dst)
