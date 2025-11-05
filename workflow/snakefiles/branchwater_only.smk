import pandas as pd
import glob, os

# OUTPUT_DIR = "/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/branchwater"
METAPG_OUT = "/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/metapangenome"
SIG, = glob_wildcards('/group/ctbrowngrp2/amhorst/2025-pigparadigm/results/metapangenome/sketches/{signature}.sig')
SIG = [s for s in SIG if s.endswith("pres1") or s.endswith("pres2")]


# # Load pang_org list
# pang_df = pd.read_csv("../resources/pang_list.csv")  # Only needs 'pang_org' column

# # Create pang_folder by replacing spaces with underscores
# pang_df["pang_folder"] = pang_df["pang_org"].str.replace(" ", "_")

# pang_folders = pang_df["pang_folder"].tolist()

# def get_species(pang_folder):
#     return pang_df.loc[pang_df["pang_folder"] == pang_folder, "pang_org"].item()


rule all:
    input:
        expand(f"{METAPG_OUT}/branchwater/{{signature}}.csv", signature=SIG),
        #expand(f"{METAPG_OUT}/branchwater_analysis/{{pang_folder}}.branchw_corehash.csv", pang_folder=pang_folders)


rule branchwater_corehash:
    input:
        sig = f"{METAPG_OUT}/sketches/{{signature}}.sig",
    output:
        csv = f"{METAPG_OUT}/branchwater/{{signature}}.csv",
    conda: "branchwater"
    threads: 1
    shell:
        """
        branchwater-client --full --sig {input.sig} -o {output.csv} --retry 10
        """


# # branchwater client
# rule branchwater:
#     input:
#         sig = f"{OUTPUT_DIR}/{{pang_folder}}/sourmash/{{pang_folder}}.merge.k21.sig",
#     output:
#         csv = f"{OUTPUT_DIR}/branchwater_analysis/{{pang_folder}}.branchw.csv",
#     conda: "branchwater"
#     threads: 1
#     shell:
#         """
#         branchwater-client --full --sig {input.sig} -o {output.csv} 
#         """