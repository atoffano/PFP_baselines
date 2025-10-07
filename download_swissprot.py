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
import argparse
import sys
import os

sys.path.append(os.path.abspath(".."))
from constants import *

BASE_PATH = "."


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


def parse_uniprot_dat(db_version, experimental_only=False):
    filepath = os.path.join("swissprot", db_version, "uniprot_sprot.dat")
    if not os.path.isfile(filepath):
        print(f"File not found: {filepath}")
        return
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
    exp_codes = {"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"}
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
                # Check for experimental evidence codes if requested
                if experimental_only:
                    match db_version:
                        case "7.0" | "4.0" | "1.0":
                            # Example : DR   GO; GO:0006270; P:DNA replication initiation; TAS.
                            m = re.match(r"DR\s+GO;\s+(GO:\d+);.*;\s+([A-Z]+)\.", line)
                        case _:
                            # DR   GO; GO:0005975; P:carbohydrate metabolic process; TAS:ProtInc.
                            m = re.match(r"DR\s+GO;\s*(GO:\d+);.*;\s*([A-Z]+):", line)
                    if m and m.group(2) in exp_codes:
                        go_terms.append(m.group(1))
                else:
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


def filter_release(experimental_only=False):
    # Removes entries from all SwissProt releases that are not present in the latest release (2024_01)
    # This is mainly to avoid leakage from proteins that could have been renamed from one version to another.
    ref_file = os.path.join(
        BASE_PATH, "swissprot/2024_01/swissprot_2024_01_annotations.tsv"
    )
    with open(ref_file) as f:
        ref_ids = set(line.split("\t")[1] for i, line in enumerate(f) if i > 0)

    for db_version in tqdm.tqdm(
        SWISSPROT_VERSIONS, desc="Filtering SwissProt releases"
    ):
        if experimental_only:
            tsv_file = os.path.join(
                BASE_PATH,
                f"swissprot/{db_version}",
                f"swissprot_{db_version}_exp_annotations.tsv",
            )
        else:
            tsv_file = os.path.join(
                BASE_PATH,
                f"swissprot/{db_version}",
                f"swissprot_{db_version}_annotations.tsv",
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


def main(experimental_only=True):
    print("Experimental only mode enabled.")
    # # Download SwissProt releases
    for db_version in tqdm.tqdm(
        SWISSPROT_VERSIONS, desc="Downloading SwissProt releases"
    ):
        if "_" in db_version:
            # Download recent SwissProt releases - from 2024 to 2010
            rel = f"release-{db_version}"
            file = f"uniprot_sprot-only{db_version}.tar.gz"
            url = f"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/{rel}/knowledgebase/{file}"
            dl_swissprot(file, url, db_version)
            print(f"Downloaded and extracted SwissProt {db_version} annotations.")
        else:
            # Download older SwissProt releases - from 2009 to 2003 (different scheme)
            file = f"uniprot_sprot-only{db_version}.tar.gz"
            url = f"https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release{rel}/knowledgebase/{file}"
            dl_swissprot(file, url, db_version)
            print(f"Downloaded and extracted SwissProt {db_version} annotations.")

    # Parse all SwissProt releases and save annotations to TSV files
    for db_version in tqdm.tqdm(SWISSPROT_VERSIONS, desc="Parsing SwissProt releases"):
        year_folder = os.path.join(BASE_PATH, "swissprot", db_version)

        entries = parse_uniprot_dat(file)
        output_file = os.path.join(
            year_folder, f"swissprot_{db_version}_annotations.tsv"
        )
        with open(output_file, "w") as out:
            out.write("EntryID\tEntry Name\tterm\tSequence\n")
            for uniprot_id, entryid, go_terms, sequence in entries:
                out.write(f"{uniprot_id}\t{entryid}\t{go_terms}\t{sequence}\n")
        print(f"Wrote {output_file}")

        entries = parse_uniprot_dat(db_version, experimental_only)
        output_file = os.path.join(
            year_folder, f"swissprot_{db_version}_exp_annotations.tsv"
        )
        with open(output_file, "w") as out:
            out.write("EntryID\tEntry Name\tterm\tSequence\n")
            for uniprot_id, entryid, go_terms, sequence in entries:
                out.write(f"{uniprot_id}\t{entryid}\t{go_terms}\t{sequence}\n")
        print(f"Wrote {output_file}")

    print("Filtering SwissProt releases to keep only entries present in 2024_01...")
    filter_release(experimental_only)

    # Create a fasta file with all sequences from the latest SwissProt release (2024_01)
    # This will be used to align all proteins against each other.
    tsv_file = os.path.join(
        BASE_PATH, f"swissprot/2024_01/swissprot_2024_01_annotations.tsv"
    )
    fasta_file = os.path.join(BASE_PATH, f"swissprot/2024_01/swissprot_2024_01.fasta")
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

    # Move fasta and tsv files to swissprot/db<db_version>
    os.makedirs(f"swissprot/{db_version}", exist_ok=True)
    shutil.move(fasta_file, f"swissprot/{db_version}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Your script description")
    parser.add_argument(
        "--experimental_only",
        action="store_true",
        help="Run experimental features only",
    )
    args = parser.parse_args()

    if args.experimental_only:
        print("Experimental only mode enabled.")
    main(args.experimental_only)
    print("SwissProt download and parsing completed successfully!")
