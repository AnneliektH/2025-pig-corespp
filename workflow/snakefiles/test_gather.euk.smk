import os
import pandas as pd

# set db path 
EUK = '/group/ctbrowngrp/sourmash-db/genbank-euks-2025.01/vertebrates.k51.sig.zip'

# set list of samples
METAG = ['SRR5371395', 'SRR5371371', 'SRR5371473', 'SRR5371477', 'SRR5371479', 'SRR5371493', 'SRR5371497', 'SRR5371447', 'SRR5371449', 'SRR5371453', 'SRR5371405', 'SRR5371419', 'SRR5371429', 'SRR5371431', 'SRR5371385', 'SRR5371387', 'SRR5371501', 'SRR5371503', 'SRR5371505', 'SRR5371509', 'SRR5371511']
METAG_COMP = ['SRR5371371', 'SRR8960211', 'ERR3211898']
PANG_OUT = "/group/ctbrowngrp2/amhorst/2025-pig-corespp/results/pangenome"


pang_df = pd.read_csv("../resources/pang_df_exEcoli.tsv", sep='\t')  # Only needs 'pang_org' column
pang_df["pang_folder"] = pang_df["pang_org"].str.replace(" ", "_")
pang_folders = pang_df["pang_folder"].tolist()

rule all:
    input:
        expand("../results/gather_euk.cow/plots_m100/{metag}.x.{pang_folder}.png", metag=METAG_COMP, pang_folder=pang_folders),
        

rule plot_abundhist:
    input:
        query="/group/ctbrowngrp/irber/data/wort-data/wort-sra/sigs/{metag}.sig",
        microbe_pang = f"{PANG_OUT}/{{pang_folder}}/{{pang_folder}}.gtdb+mags.k31.zip",
    output:
        png = "../results/gather_euk.cow/plots_m100/{metag}.x.{pang_folder}.png",
    conda: 
        "branchwater"
    threads: 1
    shell: """
       sourmash scripts abundhist {input.query} -I {input.microbe_pang} \
       --figure {output.png} --bins 100 --min 100
    """

rule gather:
    input:
        query="/group/ctbrowngrp/irber/data/wort-data/wort-sra/sigs/{metag}.sig"
    output:
        csv = "../results/gather_euk/{metag}.euk.k51.csv",
    conda: 
        "branchwater-skipmer"
    threads: 1
    shell: """
       sourmash gather {input.query} {EUK} -k 51 --scaled 10000 -o {output.csv} 
    """

