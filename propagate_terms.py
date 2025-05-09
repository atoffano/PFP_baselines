import os
import obonet
import networkx as nx
import tqdm
import pandas as pd

os.chdir("/home/atoffano/these-antoine/")
from utils.preprocessing import obsolete_terms, alt_id_terms
from utils.ia import propagate_terms, clean_ontology_edges, fetch_aspect

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
BASE_PATH = "/home/atoffano/PFP_baselines"

# Load GO ontology
obo_file_path = "/home/atoffano/these-antoine/data/ontologies/go.obo"
ontology_graph = obonet.read_obo(obo_file_path)
ontology_graph = clean_ontology_edges(ontology_graph)

# Get obsolete and alternative terms
obsolete, old_to_new = obsolete_terms(obo_file_path)
alt_to_term = alt_id_terms(obo_file_path)

# Get subontologies and mappings
roots = {"BPO": "GO:0008150", "CCO": "GO:0005575", "MFO": "GO:0003674"}
subontologies = {
    aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots
}
aspect2go = {aspect: list(subontologies[aspect].nodes) for aspect in roots}
go2aspect = {go: aspect for aspect, go_list in aspect2go.items() for go in go_list}

ontologies = {"BPO": "bp", "CCO": "cc", "MFO": "mf"}


for date in tqdm.tqdm(dates):
    print(f"Propagating terms from version: {date}...")
    tsv_file = os.path.join(BASE_PATH, f"{date}", f"swissprot_{date}_annotations.tsv")
    # Load df
    swissprot = pd.read_csv(tsv_file, sep="\t")
    swissprot = swissprot[["Entry Name", "term"]]
    # Rename Entry Name to EntryID
    swissprot = swissprot.rename(columns={"Entry Name": "EntryID"})
    swissprot["term"] = swissprot["term"].str.replace(" ", "").str.split(";")

    df_exploded = swissprot.explode("term").dropna()
    df_exploded["aspect"] = df_exploded["term"].map(go2aspect)
    df_exploded = df_exploded.dropna(subset=["aspect"])

    df_exploded = df_exploded.drop_duplicates()
    df_exploded = df_exploded[df_exploded["term"].notna()]
    print(df_exploded)
    # Split according to aspect
    for aspect in ["BPO", "CCO", "MFO"]:
        df_aspect = df_exploded[df_exploded["aspect"] == aspect].drop_duplicates()
        # Propagate annotations
        print("Propagating annotations...")
        df_prop = propagate_terms(df_aspect, {aspect: subontologies[aspect]})

        # Group terms by protein
        df_grouped = df_prop.groupby("EntryID")["term"].apply(list).reset_index()

        # Save to TSV
        output_filename = os.path.join(
            BASE_PATH,
            f"{date}",
            f"swissprot_{date}_{aspect}_annotations.tsv",
        )
        df_grouped.to_csv(output_filename, sep="\t", index=False)
        print(f"Saved {output_filename}")
