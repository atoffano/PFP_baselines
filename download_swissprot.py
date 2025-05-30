import os
import re
import tqdm
import requests
import tarfile
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BASE_PATH = "."
DB_VERSIONS = [
    "2024_01",
    "2023_01",
    "2022_01",
    "2021_01",
    "2020_01",
    "2019_01",
    "2018_01",
    "2017_01",
    "2016_01",
    "2015_01",
    "2014_01",
    "2013_01",
    "2012_01",
    "2011_01",
    "2010_01",
    "15.0",
    "13.0",
    "10.0",
    "7.0",
    "4.0",
    "1.0",
]


def dl_swissprot(file, url, db_version):
    # Download if not already present
    if not os.path.isfile(file):
        print(f"Downloading {file}...")
        try:
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
        except Exception as e:
            print(f"Failed to download {url}: {e}")
            return

    # Extract tar.gz
    print(f"Extracting {file}...")
    if not os.path.isdir(db_version):
        os.makedirs(db_version)
        with tarfile.open(file, "r:gz") as tar:
            tar.extractall(path=db_version)

    # Find and gunzip uniprot_sprot.dat.gz
    print(f"Finding uniprot_sprot.dat.gz in {db_version}...")
    datgz_path = None
    for root, _, files in os.walk(db_version):
        for fname in files:
            if fname == "uniprot_sprot.dat.gz":
                datgz_path = os.path.join(root, fname)
                break
        if datgz_path:
            break

    if not datgz_path:
        print(f"No uniprot_sprot.dat.gz found for {db_version}")
        return

    dat_path = datgz_path[:-3]  # Remove .gz
    print(f"Found {datgz_path}, extracting to {dat_path}...")
    if not os.path.isfile(dat_path):
        with gzip.open(datgz_path, "rb") as f_in, open(dat_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Extracted {dat_path}")

    # Delete all .gz files in outdir
    for root, _, files in os.walk(db_version):
        for fname in files:
            if fname.endswith(".gz"):
                os.remove(os.path.join(root, fname))


for db_version in tqdm.tqdm(DB_VERSIONS, desc="Downloading SwissProt releases"):
    if "_" in db_version:
        # Download recent SwissProt releases - from 2024 to 2010
        rel = f"release-{db_version}_01"
        file = f"uniprot_sprot-only{db_version}.tar.gz"
        url = f"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/{rel}/knowledgebase/{file}"
        dl_swissprot(file, url, db_version)
        print(f"Downloaded and extracted SwissProt {db_version}_01 annotations.")
    else:
        # Download older SwissProt releases - from 2009 to 2003
        file = f"uniprot_sprot-only{db_version}.tar.gz"
        url = f"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release{rel}/knowledgebase/{file}"
        dl_swissprot(file, url, db_version)
        print(f"Downloaded and extracted SwissProt {db_version} annotations.")


def parse_uniprot_dat(filepath):
    results = []
    print(f"Parsing {filepath} entries...")
    with open(filepath, "r") as f:
        entry = []
        for line in f:
            if line.strip() == "//":
                results.append(entry)
                entry = []
            else:
                entry.append(line.rstrip())
    parsed = []
    print(f"Parsed {len(results)} entries.")
    for entry in tqdm.tqdm(results, desc="Parsing entries content", unit="entry"):
        uniprot_id = None
        go_terms = []
        sequence = ""
        in_seq = False
        for line in entry:
            if line.startswith("ID "):
                # example line: ID 001R_FRG3G Reviewed; 256 AA.
                uniprot_id = line.split()[1]
            elif line.startswith("AC "):
                entryid = line.split()[1].strip(";")
            elif line.startswith("DR   GO;"):
                # example line: DR   GO; GO:0046782; P:regulation of viral transcription; IEA:InterPro.
                m = re.match(r"DR\s+GO;\s*(GO:\d+);", line)
                if m:
                    go_terms.append(m.group(1))
            elif line.startswith("SQ "):
                in_seq = True
                continue
            elif in_seq:
                if line.strip() == "":
                    continue
                sequence += "".join(line.strip().split())
        if uniprot_id:
            parsed.append((uniprot_id, entryid, "; ".join(go_terms), sequence))
    return parsed


# Parse all SwissProt releases and save annotations to TSV files
for db_version in tqdm.tqdm(DB_VERSIONS, desc="Parsing SwissProt releases"):
    year_folder = os.path.join(BASE_PATH, db_version)
    dat_file = os.path.join(year_folder, "uniprot_sprot.dat")
    output_file = os.path.join(year_folder, f"swissprot_{db_version}_annotations.tsv")
    if not os.path.isfile(dat_file):
        print(f"File not found: {dat_file}")
        continue
    entries = parse_uniprot_dat(dat_file)
    with open(output_file, "w") as out:
        out.write("EntryID\tEntry Name\tterm\tSequence\n")
        for uniprot_id, entryid, go_terms, sequence in entries:
            out.write(f"{uniprot_id}\t{entryid}\t{go_terms}\t{sequence}\n")
    print(f"Wrote {output_file}")

# Removes entries from all SwissProt releases that are not present in the latest release (2024_01)
# This is mainly to avoid leakage from proteins that could have been renamed.
ref_file = os.path.join(BASE_PATH, "2024_01/swissprot_2024_01_annotations.tsv")
with open(ref_file) as f:
    ref_ids = set(line.split("\t")[1] for i, line in enumerate(f) if i > 0)

for db_version in tqdm.tqdm(DB_VERSIONS, desc="Filtering SwissProt releases"):
    tsv_file = os.path.join(
        BASE_PATH, f"{db_version}", f"swissprot_{db_version}_annotations.tsv"
    )
    if not os.path.isfile(tsv_file):
        continue
    # Read all lines
    with open(tsv_file) as f:
        lines = f.readlines()
    header = lines[0]
    filtered_lines = [header]
    for line in lines[1:]:
        entry_id = line.split("\t")[1]
        if entry_id in ref_ids:
            filtered_lines.append(line)
    # Overwrite file with filtered entries
    with open(tsv_file, "w") as f:
        f.writelines(filtered_lines)
    print(
        f"Filtered {tsv_file}: {len(filtered_lines)-1} entries kept among {len(lines)-1} original entries."
    )

# Create a fasta file with all sequences from the latest SwissProt release (2024_01)
# This will be used to align all proteins against each other.
tsv_file = os.path.join(BASE_PATH, f"2024_01/swissprot_2024_01_annotations.tsv")
fasta_file = os.path.join(BASE_PATH, f"2024_01/swissprot_2024_01.fasta")
if not os.path.isfile(tsv_file):
    print(f"TSV not found: {tsv_file}")
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
