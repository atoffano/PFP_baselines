import obonet
import networkx as nx
import tqdm
import pandas as pd
import argparse

DB_VERSIONS = [
    # "2024_01",
    # "2023_01",
    # "2022_01",
    # "2021_01",
    # "2020_01",
    # "2019_01",
    # "2018_01",
    # "2017_01",
    # "2016_01",
    # "2015_01",
    # "2014_01",
    # "2013_01",
    # "2012_01",
    # "2011_01",
    # "2010_01",
    # "15.0",
    # "13.0",
    # "10.0",
    "7.0",
    "4.0",
    "1.0",
]


def propagate_terms(terms_df, subontologies):
    """
    Propagate terms in DataFrame terms_df abbording to the structure in subontologies.
    If terms were already propagated with the same graph, the returned dataframe will be equivalent to the input

    :param terms_df: pandas DataFrame of annotated terms (column names 'EntryID', 'term' 'aspect')
    :param subontologies: dict of ontology aspects (networkx DiGraphs or MultiDiGraphs)
    """

    # Look up ancestors ahead of time for efficiency
    subont_terms = {
        aspect: set(terms_df[terms_df.aspect == aspect].term.values)
        for aspect in subontologies.keys()
    }
    ancestor_lookup = {
        aspect: {
            t: nx.descendants(subont, t) for t in subont_terms[aspect] if t in subont
        }
        for aspect, subont in subontologies.items()
    }

    propagated_terms = []
    for (protein, aspect), entry_df in tqdm.tqdm(
        terms_df.groupby(["EntryID", "aspect"]),
        desc="Propagating terms",
        total=len(terms_df),
    ):
        protein_terms = set().union(
            *[list(ancestor_lookup[aspect][t]) + [t] for t in set(entry_df.term.values)]
        )

        propagated_terms += [
            {"EntryID": protein, "term": t, "aspect": aspect} for t in protein_terms
        ]

    return pd.DataFrame(propagated_terms)


def clean_ontology_edges(ontology):
    """
    Remove all ontology edges except types "is_a" and "part_of" and ensure there are no inter-ontology edges
    :param ontology: Ontology stucture (networkx DiGraph or MultiDiGraph)
    """

    # keep only "is_a" and "part_of" edges (All the "regulates" edges are in BPO)
    remove_edges = [
        (i, j, k) for i, j, k in ontology.edges if not (k == "is_a" or k == "part_of")
    ]

    ontology.remove_edges_from(remove_edges)

    # There should not be any cross-ontology edges, but we verify here
    crossont_edges = [
        (i, j, k)
        for i, j, k in ontology.edges
        if ontology.nodes[i]["namespace"] != ontology.nodes[j]["namespace"]
    ]
    if len(crossont_edges) > 0:
        ontology.remove_edges_from(crossont_edges)

    return ontology


def fetch_aspect(ontology, root: str):
    """
    Return a subgraph of an ontology starting at node <root>

    :param ontology: Ontology stucture (networkx DiGraph or MultiDiGraph)
    :param root: node name (GO term) to start subgraph
    """

    namespace = ontology.nodes[root]["namespace"]
    aspect_nodes = [
        n for n, v in ontology.nodes(data=True) if v["namespace"] == namespace
    ]
    subont_ = ontology.subgraph(aspect_nodes)
    return subont_


def main(experimental_only=False):
    """
    Main function to propagate GO terms from SwissProt annotations.
    It reads the GO ontology, processes SwissProt releases, and saves propagated annotations.
    """
    # Load GO ontology
    obo_file_path = "./data/go.obo"
    ontology_graph = obonet.read_obo(obo_file_path)
    ontology_graph = clean_ontology_edges(ontology_graph)

    # Get subontologies and mappings
    roots = {"BPO": "GO:0008150", "CCO": "GO:0005575", "MFO": "GO:0003674"}
    subontologies = {
        aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots
    }
    aspect2go = {aspect: list(subontologies[aspect].nodes) for aspect in roots}
    go2aspect = {go: aspect for aspect, go_list in aspect2go.items() for go in go_list}

    for db_version in tqdm.tqdm(DB_VERSIONS, desc="Processing SwissProt releases"):
        print(f"Propagating terms from version: {db_version}...")
        if experimental_only:
            tsv_file = f"./data/swissprot/{db_version}/swissprot_{db_version}_exp_annotations.tsv"
        else:
            tsv_file = (
                f"./data/swissprot/{db_version}/swissprot_{db_version}_annotations.tsv"
            )
        # Load df
        swissprot = pd.read_csv(tsv_file, sep="\t")
        swissprot = swissprot[["Entry Name", "term"]]
        swissprot = swissprot.rename(columns={"Entry Name": "EntryID"})
        swissprot["term"] = swissprot["term"].str.replace(" ", "").str.split(";")

        df_exploded = swissprot.explode("term").dropna()
        df_exploded["aspect"] = df_exploded["term"].map(go2aspect)
        df_exploded = df_exploded.dropna(subset=["aspect"])

        df_exploded = df_exploded.drop_duplicates()
        df_exploded = df_exploded[df_exploded["term"].notna()]
        print(df_exploded)
        # Split by aspect and propagate annotations
        for aspect in ["BPO", "CCO", "MFO"]:
            df_aspect = df_exploded[df_exploded["aspect"] == aspect].drop_duplicates()
            # Propagate annotations
            print("Propagating annotations...")
            df_prop = propagate_terms(df_aspect, {aspect: subontologies[aspect]})

            # Group by protein and join terms
            df_grouped = df_prop.groupby("EntryID")["term"].apply(list).reset_index()
            df_grouped["term"] = df_grouped["term"].apply(
                lambda x: "; ".join(map(str, x))
            )
            if experimental_only:
                output_filename = f"./data/swissprot/{db_version}/swissprot_{db_version}_{aspect}_exp_annotations.tsv"
            else:
                output_filename = f"./data/swissprot/{db_version}/swissprot_{db_version}_{aspect}_annotations.tsv"
            df_grouped.to_csv(output_filename, sep="\t", index=False)
            print(f"Saved {output_filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Propagate GO terms from SwissProt annotations."
    )
    parser.add_argument(
        "--experimental_only",
        action="store_true",
        help="Use only experimental annotations.",
    )
    args = parser.parse_args()
    main(experimental_only=args.experimental_only)
