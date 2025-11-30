import os
import pandas as pd

# set db path 
EUK = '/group/ctbrowngrp/sourmash-db/genbank-euks-2025.01/vertebrates.k31.sig.zip'
OUT_FOLDER = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/euk_gather"
PANG_OUT = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/pangenome"


pang_df = pd.read_csv("../resources/pang_df_exEcoli.tsv", sep='\t')  # Only needs 'pang_org' column
pang_df["pang_folder"] = pang_df["pang_org"].str.replace(" ", "_")
pang_folders = pang_df["pang_folder"].tolist()

rule all:
    input:
        expand(f"{OUT_FOLDER}/{{pang_folder}}.euk.k31.csv", pang_folder=pang_folders),
        

rule gather:
    input:
        query=f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb+mags.k31.zip"
    output:
        csv = f"{OUT_FOLDER}/{{pang_folder}}.euk.k31.csv",
    conda: 
        "branchwater-skipmer"
    threads: 1
    shell: """
       sourmash gather {input.query} {EUK} -k 31 --scaled 1000 --threshold-bp 0 -o {output.csv} 
    """

