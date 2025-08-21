import pandas as pd
from collections import defaultdict
import tqdm


def score(E):
    """
    Compute the Diamond Score for a query sequence against GO annotation terms.

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results filtered by e-value of 0.001,
                    with columns 'subject_id', 'bit_score' and 'subject_annotations'.
    """
    go_terms = defaultdict(float)
    for _, row in E.iterrows():
        subject_annotations = row["subject_annotations"]
        bit_score = row["bit_score"]
        for annotation in subject_annotations:
            go_terms[annotation] += bit_score
    return go_terms


def alignment_score(E):
    """
    Compute the Diamond Score for a query sequence against GO annotation terms.

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results filtered by e-value of 0.001,
                    with columns 'subject_id', 'bit_score' and 'subject_annotations'.
    """
    go_terms = score(E)
    # Normalize scores by the sum of bit scores
    for key in go_terms:
        go_terms[key] = go_terms[key] / E["bit_score"].sum()

    return go_terms


def alignment_knn(E, k=5):
    """
    Transfer annotations from the k most similar proteins based on bit score.

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results,
                  with columns 'subject_id', 'bit_score' and 'subject_annotations'.
    k (int): Number of nearest neighbors to consider
    """

    # Sort by bit_score (descending) and get top k
    top_k = E.sort_values(by="bit_score", ascending=False).head(k)
    go_terms = score(top_k)

    # Normalize scores
    total_score = top_k["bit_score"].sum()
    if total_score > 0:  # Avoid division by zero
        for key in go_terms:
            go_terms[key] = go_terms[key] / total_score

    return go_terms


def best_percent_identity(E):
    """
    Get the best percent identity from the alignments. Transfer its annotations to current protein

    Parameters:
    E (dataframe): Dataframe of similar sequences according to Diamond blast results,
                  with columns 'subject_id', 'perc_identity'.
    """
    if E.empty:
        return {}

    # Get the row with the maximum percent identity
    best_row = E.loc[E["perc_identity"].idxmax()]
    annotations = best_row["subject_annotations"]
    go_terms = {term: 1.0 for term in annotations}  # Assign a score of 1.0 to each term

    return go_terms


def naive_baseline(input_dir, train, val):
    # Compute the frequency of each GO term across all proteins in the training set
    go_term_counts = train["term"].value_counts()
    go_term_scores = go_term_counts / train["EntryID"].nunique()

    # Assign the same scores to all query proteins
    query_proteins = val["EntryID"].unique()
    predictions = pd.DataFrame(
        [
            (protein, term, score)
            for protein in query_proteins
            for term, score in go_term_scores.items()
        ],
        columns=["target_ID", "term_ID", "score"],
    )
    predictions.to_csv(
        f"{input_dir}/predictions/NaiveBaseline/predictions.tsv",
        sep="\t",
        index=False,
    )


def transfer_annotations(
    logger, pairwise_alignment, train, test, k_values, one_vs_all=False
):
    # Annotations for each sequence in the known protein set (used to transfer annotations)
    # Ts: A dictionary of annotations like {'protein_name': ['go_term1', 'go_term2', ...], ...}
    Ts = train.groupby("EntryID")["term"].apply(list).to_dict()

    pairwise_alignment["subject_annotations"] = pairwise_alignment["subject_id"].map(Ts)
    pairwise_alignment = pairwise_alignment[
        pairwise_alignment["subject_annotations"].map(
            lambda x: len(x) if isinstance(x, list) else 0
        )
        > 0
    ]  # Drop rows without annotations

    grouped = pairwise_alignment.groupby(
        "query_id"
    )  # Group by test protein for faster processing

    ascore_pred, idscore_pred = [], []
    blastknn_preds_dict = {k: [] for k in k_values}

    unaligned_proteins = 0
    unaligned_protein_ids = []  # Store unannotated protein IDs
    for protein in tqdm.tqdm(
        test["EntryID"].unique(), desc="Computing Diamond-based predictions"
    ):
        try:
            group = grouped.get_group(protein)
        except KeyError:  # No alignments for this protein
            unaligned_proteins += 1
            unaligned_protein_ids.append(protein)
            continue
        if not one_vs_all:
            try:
                assert (
                    not group["subject_id"].isin(test["EntryID"].unique()).any()
                ), "Annotation leakage has been found beetween protein sets !"
            except AssertionError as e:
                logger.warning(f"Warning for protein {protein}: {e}")
                logger.warning(
                    f"Leakage in:\n{group[group['subject_id'].isin(test['EntryID'].unique())]}"
                )
                exit(1)

        # CAFA3 baseline: Best percent identity
        idscore = best_percent_identity(group)
        if idscore:
            idscore_pred.extend(
                {"target_ID": protein, "term_ID": term_id, "score": score}
                for term_id, score in idscore.items()
            )

        # Alignment Score, DiamondKNN (based off bitscore)
        ascore_preds = alignment_score(group)  # Compute from all alignments
        for k in k_values:  # Compute from k closest alignments
            knn_preds = alignment_knn(group, k=k)

            blastknn_preds_dict[k].extend(
                [
                    {"target_ID": protein, "term_ID": term_id, "score": score}
                    for term_id, score in knn_preds.items()
                ]
            )

        ascore_pred.extend(
            [
                {"target_ID": protein, "term_ID": term_id, "score": score}
                for term_id, score in ascore_preds.items()
            ]
        )
    logger.info(
        f"Number of unaligned proteins: {unaligned_proteins} out of {len(test['EntryID'].unique())} ({unaligned_proteins / test['EntryID'].nunique() * 100} %); No annotations have been transfered for alignment-based methods."
    )
    return unaligned_protein_ids, ascore_pred, blastknn_preds_dict, idscore_pred
