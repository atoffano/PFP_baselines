import os

BASE_PATH = "/home/atoffano/PFP_baselines"

dates = ['2024_01', '2023_01', '2022_01', '2021_01', '2020_01', '2019_01', '2018_01', '2017_01', '2016_01', '2015_01', '2014_01', '2013_01', '2012_01', '2011_01', '2010_01', '2009_03', '2008_01', '2007_03', '2006_02', '2005_01', '2003_12']

# Ref file: last version
ref_file = os.path.join(BASE_PATH, "2024_01/swissprot_2024_01_annotations.tsv")
with open(ref_file) as f:
    ref_ids = set(line.split('\t')[1] for i, line in enumerate(f) if i > 0)

for date in dates:
    tsv_file = os.path.join(BASE_PATH, f"{date}", f"swissprot_{date}_annotations.tsv")
    if not os.path.isfile(tsv_file):
        continue
    # Read all lines
    with open(tsv_file) as f:
        lines = f.readlines()
    header = lines[0]
    filtered_lines = [header]
    for line in lines[1:]:
        entry_id = line.split('\t')[1]
        if entry_id in ref_ids:
            filtered_lines.append(line)
    # Overwrite file with filtered entries
    with open(tsv_file, "w") as f:
        f.writelines(filtered_lines)
    print(f"Filtered {tsv_file}: {len(filtered_lines)-1} entries kept among {len(lines)-1} original entries.")
