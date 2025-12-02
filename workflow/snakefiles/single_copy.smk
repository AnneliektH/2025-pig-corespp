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

OUTPUT_DIR = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/single_copy_genes"
MAG_LOCATION = "/group/ctbrowngrp2/scratch/annie/2023-swine-sra/results/MAGs/genomes/all_genomes/"

#genomes_all = glob_wildcards(f"{OUTPUT_DIR}/{{pang_folder}}/MAGs/{{sample}}.fasta").sample


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

GENOMES = {
    pang_folder.name: list(pang_folder.glob("MAGs/*.fasta"))
    for pang_folder in Path(OUTPUT_DIR).iterdir()
    if (pang_folder / "MAGs").exists()
}

rule all:
    input:
        # expand(f"{OUTPUT_DIR}/{{pang_folder}}/check/symlink_v3.check", pang_folder=pang_folders),
        # expand(f"{OUTPUT_DIR}/{{pang_folder}}/fasta/{{pang_folder}}.gtdb.zip", pang_folder=pang_folders),
        expand(f"{OUTPUT_DIR}/{{pang_folder}}/check/prokka.done", pang_folder=pang_folders),

# Get MAGs from specific speices 
rule get_species_mags:
    output:
        csvown = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xownmags.csv",
        csvgtdb = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xgtdb.csv",
    conda:
        "branchwater-skipmer" 
    threads: 1
    params:
        species=lambda w: get_species(w.pang_folder)
    shell:
        """
        python scripts/250903_extract_lineages.py {GTDB_TAX} "{params.species}" {output.csvgtdb} && \
        python scripts/250903_extract_lineages.py {MAG_TAX} "{params.species}" {output.csvown}
        """

rule get_MAGS:
    input:
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xgtdb.csv",
    output:
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xdirectsketch.csv",
        sig = f"{OUTPUT_DIR}/{{pang_folder}}/fasta/{{pang_folder}}.gtdb.zip",
        failed_test = f"{OUTPUT_DIR}/{{pang_folder}}/check/{{pang_folder}}.failed.csv",
        fail_checksum= f"{OUTPUT_DIR}/{{pang_folder}}/check/{{pang_folder}}.checksum.failed.csv",
    conda: 
        "branchwater-skipmer" 
    threads: 10
    params:
        output_folder=f"{OUTPUT_DIR}/{{pang_folder}}/MAGs"
    shell:
        """
        python scripts/create_acc.py {input.csv} {output.csv} && \
        sourmash scripts gbsketch  --keep-fasta --genomes-only \
        {output.csv} -o {output.sig} -p dna,k=21,k=31,scaled=1000,abund \
        -f {params.output_folder} -k -c {threads} -n 5 \
        --failed {output.failed_test} -r 5 --checksum-fail {output.fail_checksum}
        for f in {params.output_folder}/*.fna.gz; do
            mv "$f" "${{f%.fna.gz}}.fasta"
        done
        """
# symlink pig-specific MAGs
rule symlink_MAGs:
    input:
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xownmags.csv",
    output:
        check = f"{OUTPUT_DIR}/{{pang_folder}}/check/symlink_v3.check",
    conda: 
        "branchwater-skipmer"
    threads: 1
    params:
        output_folder=f"{OUTPUT_DIR}/{{pang_folder}}/MAGs"
    shell:
        """ 
        python scripts/create_symlink_midgie.py {params.output_folder} {MIDGIE_RENAME} {input.csv} && \
        touch {output.check}
        """
# prodigal or prokka
rule prokka:
    input:
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xdirectsketch.csv",
        check = f"{OUTPUT_DIR}/{{pang_folder}}/check/symlink_v3.check",
    output:
        check = f"{OUTPUT_DIR}/{{pang_folder}}/check/prokka.done",
    conda: 
        "prokka"
    threads: 1
    params:
        input_folder=f"{OUTPUT_DIR}/{{pang_folder}}/MAGs/done",
        done_folder = f"{OUTPUT_DIR}/{{pang_folder}}/MAGs",
        output_folder=f"{OUTPUT_DIR}/{{pang_folder}}/prokka"
    shell:
        """
        mkdir -p {params.done_folder}
        mkdir -p {params.output_folder}
        for f in {params.input_folder}/*.fasta; do
            base=$(basename "$f" .fasta)
            prokka --kingdom Bacteria \
                --outdir {params.output_folder}/$base \
                --norrna --notrna --prefix $base --force \
                --locustag $base "$f" && \
                mv $f {params.done_folder}
        done
        touch {output.check}
        """

rule fetch_mg:
    input:
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}xdirectsketch.csv",
        check = f"{OUTPUT_DIR}/{{pang_folder}}/check/symlink_v3.check",
    output:
        check = f"{OUTPUT_DIR}/{{pang_folder}}/check/fetchmg.done",
    conda: 
        "fetch_mgs"
    threads: 12
    params:
        input_folder=f"{OUTPUT_DIR}/{{pang_folder}}/MAGs",
        done_folder = f"{OUTPUT_DIR}/{{pang_folder}}/MAGs/done",
        output_folder=f"{OUTPUT_DIR}/{{pang_folder}}/fetch_mg"
    shell:
        """
        mkdir -p {params.done_folder}
        mkdir -p {params.output_folder}
        for f in {params.input_folder}/*.fasta; do
            base=$(basename "$f" .fasta)
            /home/amhorst/.local/bin/fetchMGs extraction \
                -t {threads} "$f" "genome" {params.output_folder}/$base && \
                mv $f {params.done_folder}
        done
        touch {output.check}
        """

# agg the .fna files from fetchmgs
rule aggregate_fna:
    input:
        done = f"{OUTPUT_DIR}/{{pang_folder}}/check/fetchmg.done",
    output:
        aggregated = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.SCG.fna",
    shell:
        r"""
        find {OUTPUT_DIR}/{wildcards.pang_folder}/fetch_mg -type f -name "*.fna" \
            -exec cat {{}} + > {output.aggregated}
        """


# # sketch the single copy genes
rule sketch_scg:
    input:
        fna = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.SCG.fna",
    output:
        sig = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.SCG.k31.sig",
    conda: 
        "branchwater"
    threads: 1
    shell:
        """ 
        sourmash sketch dna -o {output.sig} {input.fna}
        """
 

 # manysearch x metagenomes
rule manysearch:
    input:
        sig = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.SCG.k31.zip",
    output:
        csv = f"{OUTPUT_DIR}/{{pang_folder}}/{{pang_folder}}.SCG.k31.manysearch.csv",
    conda: 
        "branchwater"
    threads: 40
    shell:
        """ 
        sourmash scripts manysearch {MGLIST} {input.sig} -o {output.csv} -c {threads} -t 0
        """
 

# # use colton script for funprofiler
# rule funcprofiler_MAGs:
#     input:
#         genomes= f"{OUTPUT_DIR}/{{pang_folder}}/MAGs/{{genome}}.fasta",
#     output:
#         csv = f"{OUTPUT_DIR}/{{pang_folder}}/funcprofiler/{{genome}}.prefetch.csv",
#     conda: 
#         "branchwater"
#     threads: 1
#     params:
#         output_folder=f"{OUTPUT_DIR}/{{pang_folder}}/funprofiler",
#     shell:
#         """ 
#         python fmh-funprofiler/funcprofiler.py \
#         {input.genomes} \
#         /group/ctbrowngrp2/amhorst/2025-pigparadigm/databases/KOs_sketched_scaled_1000.sig.zip \
#         11 1000 {params.output_folder} -p {output.csv}
#         """