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
OUTPUT_DIR = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/prokka"
MAG_DIR = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/single_copy_genes"

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
    for pang_folder in Path(MAG_DIR).iterdir()
    if (pang_folder / "MAGs").exists()
}

rule all:
    input:
        # expand(f"{OUTPUT_DIR}/{{pang_folder}}/check/symlink_v3.check", pang_folder=pang_folders),
        expand(f"{OUTPUT_DIR}/{{pang_folder}}.coding.k31.sig", pang_folder=pang_folders),
        expand(f"{OUTPUT_DIR}/{{pang_folder}}.coding.fa", pang_folder=pang_folders),



# FASTAS, = glob_wildcards("/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/single_copy_genes/Phocaeicola_vulgatus/MAGs/done/{sample}.fasta")

# rule all_prokka:
#     input:
#         expand("/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/prokka/Phocaeicola_vulgatus/check/{sample}.done", sample=FASTAS)


# # Get MAGs from specific speices 
# rule get_species_mags:
#     output:
#         csvown = f"{MAG_DIR}/{{pang_folder}}/{{pang_folder}}xownmags.csv",
#         csvgtdb = f"{MAG_DIR}/{{pang_folder}}/{{pang_folder}}xgtdb.csv",
#     conda:
#         "branchwater-skipmer" 
#     threads: 1
#     params:
#         species=lambda w: get_species(w.pang_folder)
#     shell:
#         """
#         python scripts/250903_extract_lineages.py {GTDB_TAX} "{params.species}" {output.csvgtdb} && \
#         python scripts/250903_extract_lineages.py {MAG_TAX} "{params.species}" {output.csvown}
#         """


# # prodigal or prokka
# rule prokka:
#     input:
#         csv = f"{MAG_DIR}/{{pang_folder}}/{{pang_folder}}xdirectsketch.csv",
#         check = f"{MAG_DIR}/{{pang_folder}}/check/symlink_v3.check",
#     output:
#         check = f"{OUTPUT_DIR}/{{pang_folder}}/check/prokka.done",
#     conda: 
#         "prokka"
#     threads: 1
#     params:
#         input_folder=f"{MAG_DIR}/{{pang_folder}}/MAGs/done",
#         done_folder = f"{MAG_DIR}/{{pang_folder}}/MAGs",
#         output_folder=f"{OUTPUT_DIR}/{{pang_folder}}"
#     shell:
#         """
#         mkdir -p {params.done_folder}
#         mkdir -p {params.output_folder}
#         for f in {params.input_folder}/*.fasta; do
#             base=$(basename "$f" .fasta)
#             prokka --kingdom Bacteria \
#                 --outdir {params.output_folder}/$base \
#                 --norrna --notrna --prefix $base --force \
#                 --locustag $base "$f" && \
#                 mv $f {params.done_folder}
#         done
#         touch {output.check}
#         """

# agg the .fna files from fetchmgs
rule aggregate_fna:
    input:
        done = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/pangenome/{pang_folder}/{pang_folder}.gtdb+mags.k21.zip",
    output:
        aggregated = f"{OUTPUT_DIR}/{{pang_folder}}.coding.fa",
    shell:
        r"""
        find {OUTPUT_DIR}/{wildcards.pang_folder}/ -type f -name "*.ffn" \
            -exec cat {{}} + > {output.aggregated}
        """

rule sketch:
    input:
        fa =f"{OUTPUT_DIR}/{{pang_folder}}.coding.fa",
    output:
        sig = f"{OUTPUT_DIR}/{{pang_folder}}.coding.k31.sig",
    conda: 
        "branchwater"
    threads: 1
    shell:
        """ 
        sourmash sketch dna -o {output.sig} {input.fa}
        """
# rule prokka_per_sample:
#     input:
#         fasta = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/single_copy_genes/Phocaeicola_vulgatus/MAGs/done/{sample}.fasta"
#     output:
#         check = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/prokka/Phocaeicola_vulgatus/check/{sample}.done"
#     conda:
#         "prokka"
#     threads: 1
#     params:
#         outdir = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/prokka/Phocaeicola_vulgatus/{sample}",
#         done_folder = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/single_copy_genes/Phocaeicola_vulgatus/MAGs/"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         prokka --kingdom Bacteria \
#             --outdir {params.outdir} \
#             --norrna --notrna --prefix {wildcards.sample} --force \
#             --locustag {wildcards.sample} {input.fasta}
#         mv {input.fasta} {params.done_folder}
#         touch {output.check}
#         """