import sys
import argparse
import obonet
import numpy as np
import pandas as pd
import networkx as nx
from collections import Counter
from scipy.sparse import dok_matrix
import tqdm


def obsolete_terms(ontology):
    """Returns a set of obsolete terms without replacements and a dict of replaceable obsolete terms with their replacements as values."""
    graph_with_obs = (
        obonet.read_obo(ontology, ignore_obsolete=False)
        if type(ontology) == str
        else ontology
    )
    print(f"Number of terms: {len(graph_with_obs)}")
    old_to_new = dict()
    obsolete = set()
    for node, data in graph_with_obs.nodes(data=True):
        for replaced_by in data.get("replaced_by", []):
            old_to_new[node] = replaced_by
        if data.get("is_obsolete", False) and node not in old_to_new.keys():
            obsolete.add(node)

    return obsolete, old_to_new


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


def term_counts(terms_df, term_indices):
    """
    Count the number of instances of each term

    :param terms_df: pandas DataFrame of (propagated) annotated terms (column names 'EntryID', 'term', 'aspect')
    :param term_indices:
    """

    num_proteins = len(terms_df.groupby("EntryID"))
    S = dok_matrix((num_proteins + 1, len(term_indices)), dtype=np.int32)
    S[-1, :] = 1  # dummy protein

    for i, (protein, protdf) in enumerate(terms_df.groupby("EntryID")):
        row_count = {term_indices[t]: c for t, c in Counter(protdf["term"]).items()}
        for col, count in row_count.items():
            S[i, col] = count

    return S


def calc_ia(term, count_matrix, ontology, terms_index):

    parents = nx.descendants_at_distance(ontology, term, 1)

    # count of proteins with term
    prots_with_term = count_matrix[:, terms_index[term]].sum()

    # count of proteins with all parents
    num_parents = len(parents)
    prots_with_parents = (
        count_matrix[:, [terms_index[p] for p in parents]].sum(1) == num_parents
    ).sum()

    # avoid floating point errors by returning exactly zero
    if prots_with_term == prots_with_parents:
        return 0

    return -np.log2(prots_with_term / prots_with_parents)


def parse_inputs(argv):
    parser = argparse.ArgumentParser(
        description="Compute Information Accretion of GO annotations. Note: If annotations in input file have been propagated to ontology roots, the input onotology graph should be the same as the one used to propagate terms"
    )

    parser.add_argument("--annot", "-a", required=True, help="Path to annotation file")

    parser.add_argument(
        "--dataset",
        "-o",
        default=None,
        help="Dataset name used to load test proteins. If empty (default), the script will not filter out test proteins",
    )
    parser.add_argument(
        "--ontology",
        "-go",
        default=None,
        help="Path to OBO ontology graph file, if local. If empty (default) current OBO structure at run-time will be downloaded from http://purl.obolibrary.org/obo/go/go-basic.obo",
    )

    parser.add_argument(
        "--prop",
        "-p",
        action="store_true",
        help="Flag to propagate terms in annotation file according to the ontology graph",
    )

    parser.add_argument(
        "--aspect",
        "-asp",
        default=None,
        help="Compute IA for terms in this aspect only. If empty (default), IA will be computed for all terms",
    )

    return parser.parse_args(argv)


if __name__ == "__main__":

    args = parse_inputs(sys.argv[1:])

    # IA should be computed using the same ontology version that was used for term propagation.
    # Otherwise, this may result in negative IA values.
    annotation_df = pd.read_csv(args.annot, sep="\t")
    annotation_df = annotation_df[["EntryID", "term"]]

    # Data loading
    if args.aspect:
        test_df = pd.read_csv(
            f"./data/{args.dataset}/{args.dataset}_{args.aspect}_test_annotations.tsv",
            sep="\t",
            header=None,
            names=["EntryID", "term"],
        )

    if args.dataset:
        # Remove test proteins from IC computation
        annotation_df = annotation_df[
            ~annotation_df["EntryID"].isin(test_df["EntryID"])
        ]

    # Ontology Loading
    if args.ontology:
        ontology_path = args.ontology
    else:
        ontology_path = "http://purl.obolibrary.org/obo/go/go.obo"
    ontology_graph = clean_ontology_edges(
        obonet.read_obo(
            ontology_path,
            ignore_obsolete=False,
        )
    )
    roots = {"BPO": "GO:0008150", "CCO": "GO:0005575", "MFO": "GO:0003674"}
    subontologies = {
        aspect: fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots
    }
    aspect = {
        "BPO": list(subontologies["BPO"].nodes),
        "CCO": list(subontologies["CCO"].nodes),
        "MFO": list(subontologies["MFO"].nodes),
    }
    obsolete, old_to_new = obsolete_terms(ontology_graph)

    # Reverse aspect dictionary
    aspect = {go: aspect for aspect, go_list in aspect.items() for go in go_list}
    annotation_df = annotation_df.dropna(subset=["term"])

    annotation_df["term"] = annotation_df["term"].apply(lambda x: x.split("; "))
    annotation_df = annotation_df[annotation_df["term"].notna()]
    annotation_df = annotation_df[["EntryID", "term"]].explode("term")
    # annotation_df["term"] = annotation_df["term"].map(lambda x: old_to_new.get(x, x))
    # annotation_df = annotation_df[~annotation_df["term"].isin(obsolete)]
    # Remove
    annotation_df["aspect"] = annotation_df["term"].map(aspect)
    annotation_df = annotation_df.dropna(subset=["aspect"])

    print(annotation_df.head())

    # Remove aspects not matching args.aspect if specified
    if args.aspect:
        subontologies = {args.aspect: subontologies[args.aspect]}
        print(f"Computing IA for aspect {args.aspect}")

    if args.prop:
        print("Propagating Terms")
        annotation_df = propagate_terms(annotation_df, subontologies)

    # Count term instances
    print("Counting Terms")
    aspect_counts = dict()
    aspect_terms = dict()
    term_idx = dict()
    for aspect, subont in subontologies.items():
        aspect_terms[aspect] = sorted(subont.nodes)  # ensure same order
        term_idx[aspect] = {t: i for i, t in enumerate(aspect_terms[aspect])}
        aspect_counts[aspect] = term_counts(
            annotation_df[annotation_df.aspect == aspect], term_idx[aspect]
        )

        assert aspect_counts[aspect].sum() == len(
            annotation_df[annotation_df.aspect == aspect]
        ) + len(aspect_terms[aspect])

    # since we are indexing by column to compute IA,
    # let's convert to Compressed Sparse Column format
    sp_matrix = {aspect: dok.tocsc() for aspect, dok in aspect_counts.items()}

    # Compute IA
    print("Computing Information Accretion")
    aspect_ia = {
        aspect: {t: 0 for t in aspect_terms[aspect]} for aspect in aspect_terms.keys()
    }
    for aspect, subontology in subontologies.items():
        for term in aspect_ia[aspect].keys():
            aspect_ia[aspect][term] = calc_ia(
                term, sp_matrix[aspect], subontology, term_idx[aspect]
            )

    ia_df = pd.concat(
        [
            pd.DataFrame.from_dict(
                {
                    "term": aspect_ia[aspect].keys(),
                    "ic": aspect_ia[aspect].values(),
                    "aspect": aspect,
                }
            )
            for aspect in subontologies.keys()
        ]
    )

    negative_ia_terms = ia_df[ia_df["ic"] < 0]
    if not negative_ia_terms.empty:
        print("The following terms have negative IA values:")
        print(negative_ia_terms)
        print(
            "This usually happens when there is a mismatch between ontology versions used for propagation and IA computation."
        )
    # all counts should be non-negative
    assert ia_df["ic"].min() >= 0

    # Save to file
    if args.aspect:
        args.outfile = f"./data/{args.dataset}/IC_{args.dataset}_{args.aspect}.tsv"
    else:
        args.outfile = f"./data/{args.dataset}/IC_{args.dataset}.tsv"
    print(f"Saving to file {args.outfile}")
    ia_df[["term", "ic"]].to_csv(args.outfile, header=None, sep="\t", index=False)

# Example usage:
# python ia.py --annot ./data/swissprot/2024_01/swissprot_2024_01_BPO_exp_annotations.tsv --dataset H30 --ontology ./data/go.obo --aspect BPO
