import pandas as pd
import glob, os

BRANCHW_OUT = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/branchwater"
PANG_OUT = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/pangenome"

# Load pang_org list
pang_df = pd.read_csv("../resources/pang_df_exEcoli.tsv", sep='\t')  # Only needs 'pang_org' column

# Create pang_folder by replacing spaces with underscores
pang_df["pang_folder"] = pang_df["pang_org"].str.replace(" ", "_")

pang_folders = pang_df["pang_folder"].tolist()

# Gives species name, e.g. Lactobacillus amylovorus
def get_species(pang_folder):
    return pang_df.loc[pang_df["pang_folder"] == pang_folder, "pang_org"].item()

# Gives exact name, e.g. GCA_004552585 s__Lactobacillus_amylovorus
def get_exact_name(pang_folder):
    return pang_df.loc[pang_df["pang_folder"] == pang_folder, "exact_name"].item()

rule all:
    input:
        expand(f"{BRANCHW_OUT}/{{pang_folder}}.gtdb+mags.k21.csv", pang_folder=pang_folders),
 

rule pangenome_merge:
    input:
        sig_gtdb_k21 = f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb.k21.zip",
        sig_gtdbmag_k21 = f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb+mags.k21.zip"
    output:
        merged_gtdb_k21=f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb.k21.pang.sig",
        merged_gtdbmag_k21=f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb+mags.k21.sig",
    conda:
        "pangenomics_dev"
    shell:
        """ 
        sourmash scripts pangenome_merge {input.sig_gtdb_k21} -k 21 -o {output.merged_gtdb_k21} --scaled 1000 && \
        sourmash scripts pangenome_merge {input.sig_gtdbmag_k21} -k 21 -o {output.merged_gtdbmag_k21} --scaled 1000
        """

rule branchwater:
    input:
        sig = f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb+mags.k21.sig",
    output:
        csv = f"{BRANCHW_OUT}/{{pang_folder}}.gtdb+mags.k21.csv",
    conda: "branchwater"
    threads: 1
    shell:
        """
        branchwater-client --full --sig {input.sig} -o {output.csv} --retry 10
        """

rule branchwater_gtdb:
    input:
        sig = f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb.k21.pang.sig",
    output:
        csv = f"{BRANCHW_OUT}/{{pang_folder}}.gtdb.k21.csv",
    conda: "branchwater"
    threads: 1
    shell:
        """
        branchwater-client --full --sig {input.sig} -o {output.csv} --retry 10
        """