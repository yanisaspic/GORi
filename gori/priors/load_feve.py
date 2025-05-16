"""Function called to load the results of a fEVE analysis as a prior for a GORi analysis.

    2025/04/30 @yanisaspic"""

import pandas as pd
from typing import Any, Callable
from gori.utils import _get_uniprot_id


def _get_genes_wrapper() -> dict[str, Callable]:
    """Get the function to parse the features of a fEVE analysis.

    Returns
        A function that takes a feature name and returns its corresponding cluster name.
    """
    genes_wrapper = {
        "up": lambda x: x > 0,
        "down": lambda x: x < 0,
        "any": lambda x: x != 0,
    }
    return genes_wrapper


def load_feve(path: str, direction: str = "any") -> dict[str, dict[str, Any]]:
    """Load a custom knowledge base from a fEVE analysis.

    ``path`` is the path to the .xlsx file containing the results of a fEVE analysis.
    ``direction`` is the direction of the genes' expression: one of 'up', 'down' or 'any'.

    Returns
        A dict with three keys:
            - "annotations": a dict associating a UniProt ID to its related clusters
            - "hierarchy": a dict associating a cluster to its parent clusters in the hierarchy
            - "translations": a placeholder for the cluster names
    """
    meta = pd.read_excel(path, sheet_name="meta", index_col=0)
    features = pd.read_excel(path, sheet_name="features", index_col=0)
    genes_wrapper = _get_genes_wrapper()
    f = genes_wrapper[direction]

    annotations = {
        _get_uniprot_id(gene): set(features.columns[f(features.loc[gene])])
        for gene in features.index
    }
    annotations = {gene: clusters for gene, clusters in annotations.items() if clusters}
    hierarchy = meta.iloc[1:].groupby(meta.iloc[1:].index).parent.apply(set).to_dict()
    translations = {clu: clu for clu in meta.index}

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }
