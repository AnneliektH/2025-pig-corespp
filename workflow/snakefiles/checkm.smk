## Get the funprofiler for all core spp pangenomes. 
import pandas as pd
import glob, os
from pathlib import Path


GTDB_K21 = '/group/ctbrowngrp/sourmash-db/gtdb-rs226/gtdb-rs226.k21.sig.zip'
GTDB_K31 = '/group/ctbrowngrp/sourmash-db/gtdb-rs226/gtdb-rs226.k31.sig.zip'
GTDB_TAX  = '/group/ctbrowngrp/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv'
MAG_TAX = '/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.lineages.csv'
MIDGIE_RENAME = '/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/manysketch-renamed.csv'
MGLIST = '/group/ctbrowngrp2/amhorst/2025-pigparadigm/resources/3217_metag.txt'
# KSIZE = 31 
SCALED = 1000

MAG_DIR = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/single_copy_genes"
OUTPUT_DIR = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/checkm"


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
        # expand(f"{OUTPUT_DIR}/{{pang_folder}}/check/symlink_v3.check", pang_folder=pang_folders),
        expand( f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.genepredict.fna", pang_folder=pang_folders),
        expand(f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.marker.k31.sig", pang_folder=pang_folders),

       



rule checkm:
    input:
        csv = f"{MAG_DIR}/{{pang_folder}}/{{pang_folder}}xdirectsketch.csv",
        check = f"{MAG_DIR}/{{pang_folder}}/check/symlink_v3.check",
    output:
        check = f"{OUTPUT_DIR}/{{pang_folder}}/checkm.done",
    conda: 
        "checkM"
    threads: 48
    params:
        input_folder=f"{MAG_DIR}/{{pang_folder}}/MAGs",
        output_folder=f"{OUTPUT_DIR}/{{pang_folder}}"
    shell:
        """
        export CHECKM_DATA_PATH=/home/amhorst/databases/checkm
        checkm analyze --nt -x fasta \
        /home/amhorst/databases/checkm/Prokaryote.ms \
        {params.input_folder} {params.output_folder} -t {threads} && \
        touch {output.check}
        """

rule checkm_qa:
    input:
        check = f"{OUTPUT_DIR}/{{pang_folder}}/checkm.done",
    output:
        check = f"{OUTPUT_DIR}/{{pang_folder}}/checkm_qa.done",
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/checkm_qa.out8.csv",
    conda: 
        "checkM"
    threads: 48
    params:
        input_folder=f"{OUTPUT_DIR}/{{pang_folder}}"
    shell:
        """
        export CHECKM_DATA_PATH=/home/amhorst/databases/checkm
        checkm qa --out_format 8 --file {output.csv} \
        /home/amhorst/databases/checkm/Prokaryote.ms \
        {params.input_folder} -t {threads} && \
        touch {output.check}
        """

# agg the .fna files from fetchmgs
rule aggregate_fna:
    input:
        done = f"{OUTPUT_DIR}/{{pang_folder}}/checkm.done",
    output:
        aggregated = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.genepredict.fna",
    shell:
        r"""
        find {OUTPUT_DIR}/{wildcards.pang_folder}/bins -type f -name "*.fna" \
            -exec cat {{}} + > {output.aggregated}
        """

rule cut_sequences:
    input:
        fasta = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.genepredict.fna",
        gene_ids = f"{OUTPUT_DIR}/{{pang_folder}}/checkm_qa.out8.csv"
    output:
        gene_ids = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.gene_id_marker.txt",
        fasta =  f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.markergenes.fna",
    conda: 
        "seqkit"
    shell:
        """
        cut -f2 {input.gene_ids} >  {output.gene_ids}
        seqkit grep -f {output.gene_ids} {input.fasta} > {output.fasta}
        """

# # sketch the single copy genes
rule sketch_scg:
    input:
        fna = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.markergenes.fna",
    output:
        sig = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.marker.k31.sig",
    conda: 
        "branchwater"
    threads: 1
    shell:
        """ 
        sourmash sketch dna -o {output.sig} {input.fna}
        """
 