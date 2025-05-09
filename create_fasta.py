from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

BASE_PATH = "/home/atoffano/PFP_baselines"
dates = [
    "2024_01",
]

for date in dates:
    tsv_file = os.path.join(BASE_PATH, f"{date}", f"swissprot_{date}_annotations.tsv")
    fasta_file = os.path.join(BASE_PATH, f"{date}", f"swissprot_{date}.fasta")
    if not os.path.isfile(tsv_file):
        print(f"TSV not found: {tsv_file}")
        continue
    # Read TSV and write FASTA
    records = []
    with open(tsv_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip().split("\t")
            entry_id, entry_name, go_terms, sequence = parts
            if not sequence:
                continue
            record = SeqRecord(Seq(sequence), id=entry_name, description="")
            records.append(record)
    SeqIO.write(records, fasta_file, "fasta")
    print(f"Wrote {len(records)} records to {fasta_file}")
