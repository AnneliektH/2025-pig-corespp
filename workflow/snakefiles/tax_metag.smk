import os
import pandas as pd

# set db path 
MAGS = '/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/gtdb_pangenomedb/sra_mags.pangenomedb_speciestaxed.10k.zip'
GTDB = '/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/gtdb_pangenomedb/gtdb-rs226.pangenomedb_species.10k.zip'
KSIZE =  31

# set list of samples
WORT_METAG = pd.read_csv("../resources/metag-wort-hq.3217.txt", usecols=[0], header=None).squeeze().tolist()

rule all:
    input:
        expand("../results/tax_metag/check/{metag}.check", metag=WORT_METAG,),

rule tax_metag:
    input:
        csv="../results/gatherxgtdb+bins.species/gather/{metag}.csv",
        tax_db1 ="../resources/gtdb-rs226.lineages.csv",
        tax_db2 ="../resources/bin-sketches.lineages.csv"
    output:
        check = "../results/tax_metag/check/{metag}.check"
    conda: 
        "branchwater-skipmer"
    threads: 1
    shell: """
       sourmash tax metagenome -g {input.csv} \
       --taxonomy {input.tax_db1} {input.tax_db2} -o {wildcards.metag} \
       --output-dir ../results/tax_metag/summary/ && touch {output.check}
    """
