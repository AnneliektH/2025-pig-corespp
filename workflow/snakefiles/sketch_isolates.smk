import pandas as pd
import glob, os


# GTDB_K21 = '/group/ctbrowngrp/sourmash-db/gtdb-rs226/gtdb-rs226.k21.sig.zip'
# GTDB_K31 = '/group/ctbrowngrp/sourmash-db/gtdb-rs226/gtdb-rs226.k31.sig.zip'
# GTDB_TAX  = '/group/ctbrowngrp/sourmash-db/gtdb-rs226/gtdb-rs226.lineages.csv'
# MAG_TAX = '/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/MAGs/250827_taxfor_sourmash.gtdb226.csv'
# OWN_MAG_SIG = '/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/sketches/MAGs.all_taxed.zip'
# # KSIZE = 31 
# SCALED = 1000
# OUTPUT_DIR = "/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/functional_profile_corespp"
# MAG_LOCATION = "/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/MAGs/all_genomes"


# Load pang_org list
species_all = glob_wildcards("/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}_filtered.csv").name



rule all:
    input:
        expand("/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}.pang.rn.sig", name=species_all),



rule get_MAGS:
    input:
        csv = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}_filtered.csv",
    output:
        sig = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}.gtdb.sig",
        failed_test = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}/{name}.failed.csv",
        fail_checksum= "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}/{name}.checksum.failed.csv",
    conda: 
        "branchwater-skipmer"  
    threads: 1
    params:
        output_folder="/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}"
    shell:
        """
        sourmash scripts gbsketch  --keep-fasta --genomes-only \
        {input.csv} -o {output.sig} -p dna,k=21,k=31,scaled=100,abund \
        -f {params.output_folder} -k -c {threads} -n 5 \
        --failed {output.failed_test} -r 5 --checksum-fail {output.fail_checksum}
        for f in {params.output_folder}/*.fna.gz; do
            mv "$f" "${{f%.fna.gz}}.fasta"
        done
        """

rule pangenome_merge_k31:
    input:
        sig =  "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}.gtdb.zip",
    output:
        merged="/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}.pang.sig",
    conda:
        "pangenomics_dev"
    shell:
        """ 
        sourmash scripts pangenome_merge {input.sig} -k 31 -o {output.merged} --scaled 1000
        """

rule sigrename:
    input:
        sig_gtdb = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}.pang.sig",
    output:
        rn_gtdb="/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/isolate_genomes/{name}.pang.rn.sig",
    conda:
        "pangenomics_dev"
    shell:
        """ 
        sourmash sig rename {input.sig_gtdb} {wildcards.name} -o {output.rn_gtdb} 

        """