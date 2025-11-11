import pandas as pd
import os
import argparse

# Argument parser setup
parser = argparse.ArgumentParser(description="Create symlinks for files listed in a CSV column.")
parser.add_argument("input", help="Path to the input CSV file")
parser.add_argument("source_folder", help="Path to the folder containing source files")
parser.add_argument("destination_folder", help="Path to the folder where symlinks will be created")
args = parser.parse_args()

# Read the CSV
df = pd.read_csv(args.input, sep=',')

# Ensure destination folder exists
os.makedirs(args.destination_folder, exist_ok=True)

# Create symlinks
for genome in df["ident"]:
    source_path = os.path.join(args.source_folder, genome)
    symlink_path = os.path.join(args.destination_folder, genome)

    if os.path.exists(source_path):
        os.symlink(source_path, symlink_path)
    else:
        print(f"Warning: {source_path} does not exist.")

print("Symlinks created.")
