import os
import re
import tqdm

BASE_PATH = "/home/atoffano/PFP_baselines"


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


dates = [
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
    "2009_03",
    "2008_01",
    "2007_03",
    "2006_02",
    "2005_01",
    "2003_12",
]
for date in dates:
    year_folder = os.path.join(BASE_PATH, date)
    dat_file = os.path.join(year_folder, "uniprot_sprot.dat")
    output_file = os.path.join(year_folder, f"swissprot_{date}_annotations.tsv")
    if not os.path.isfile(dat_file):
        print(f"File not found: {dat_file}")
        continue
    entries = parse_uniprot_dat(dat_file)
    with open(output_file, "w") as out:
        out.write("EntryID\tEntry Name\tterm\tSequence\n")
        for uniprot_id, entryid, go_terms, sequence in entries:
            out.write(f"{uniprot_id}\t{entryid}\t{go_terms}\t{sequence}\n")
    print(f"Wrote {output_file}")
