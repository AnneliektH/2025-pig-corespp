import pandas as pd
import sys

# usage: python script.py input.csv "Escherichia coli" output.csv --add-fasta
input_file = sys.argv[1]
species_name = sys.argv[2]
output_file = sys.argv[3]
add_fasta = "--add-fasta" in sys.argv

df_combined = pd.read_csv('/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/csv_files/251201_ident.contaminated.csv')
df = pd.read_csv(input_file)

filtered = df[df["species"] == f"s__{species_name}"]
filtered = filtered[~filtered["ident"].isin(df_combined["ident"])]

if add_fasta:
    filtered["ident"] = filtered["ident"].astype(str) + ".fasta"

filtered.to_csv(output_file, index=False)
