import subprocess

# Run Diamond blast on protein sequences against themselves
diamond_path = "/home/atoffano/these-antoine/utils"
datapath = "/home/atoffano/PFP_baselines/2024_01"
print("Creating Diamond database...")
# subprocess.run(
#     f"{diamond_path}/diamond makedb --in {datapath}/swissprot_2024_01.fasta -d {datapath}/swissprot_2024_01_proteins_set",
#     shell=True,
#     check=True,
# )

print("Running Diamond blast on protein sequences against themselves...")
subprocess.run(
    f"{diamond_path}/diamond blastp --very-sensitive --db {datapath}/swissprot_2024_01_proteins_set.dmnd --query {datapath}/swissprot_2024_01.fasta --out {datapath}/diamond_swissprot_2024_01_alignment.tsv -e 0.001",
    shell=True,
    check=True,
)
