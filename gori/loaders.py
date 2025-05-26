"""Functions called to load knowledge bases prior to the GORi analysis.

    2025/05/20 @yanisaspic"""

import os
import json
import pandas as pd
from typing import Any, Callable
from gori.params import get_parameters
from gori.src.utils import _get_uniprot_id


def load_priors(
    priors: set[str], path: str = "./priors", params: dict[str, Any] = get_parameters()
) -> dict[str, Any]:
    """Load a subset of knowledge bases readily available in the GORi package.

    ``priors`` is a set of strings containing the names of the knowledge bases to load, including:
        `BIOP` (biological processes),
        `CELC` (cellular components),
        `CTYP` (cell types),
        `DISE` (diseases),
        `GENG` (gene groups),
        `MOLF` (molecular functions),
        `PATH` (pathways) and
        `PHEN` (phenotypes).
    ``path`` is the path to the JSON files containing the knowledge bases.
    ``params`` is a dict of parameters.

    Returns
        A dict associating knowledge base names (keys) to their contents (values).
    """
    load_wrapper = params["wrappers"]["load_wrapper"]
    if any(p in priors for p in ["BIOP", "CELC", "MOLF"]):
        go = load_wrapper["BIOP"](
            path
        )  # outputs the same GOAnnotation as CELC and MOLF

    knowledge_bases = {}
    for p in priors:
        if p not in load_wrapper.keys():
            raise ValueError(
                f"Invalid prior: {p}. Valid priors are: {set(load_wrapper)}"
            )
        elif p in ["BIOP", "CELC", "MOLF"]:
            knowledge_bases[p] = go
        else:
            knowledge_bases[p] = load_wrapper[p](path)
    return knowledge_bases


def load_feve(path: str, direction: str = "any") -> dict[str, dict[str, Any]]:
    """Load a custom knowledge base from a fEVE analysis.

    ``path`` is the path to the .xlsx file containing the results of a fEVE analysis.
    ``direction`` is the direction of the genes' expression: one of 'up', 'down' or 'any'.

    Returns
        A dict with three keys:
            ``annotations``: a dict associating a UniProt ID to its related clusters
            ``hierarchy``: a dict associating a cluster to its parent clusters in the hierarchy
            ``translations``: a placeholder for the cluster names
    """
    meta = pd.read_excel(path, sheet_name="meta", index_col=0)
    features = pd.read_excel(path, sheet_name="features", index_col=0)

    genes_wrapper = {
        "up": lambda x: x > 0,
        "down": lambda x: x < 0,
        "any": lambda x: x != 0,
    }
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


def load_local(prior: str, path: str) -> dict[str, dict[str, Any]]:
    """Load a local knowledge base.

    ``prior`` is the label of the knowledge base to add.
    ``path`` is the path to the JSON files containing the annotations, the hierarchy
    and the translations of the knowledge base.

    Returns
        A dict with three keys:
            ``annotations``: a dict associating a Uniprot ID to its related concepts
            ``ontology``: a dict associating a concept to its parent concepts in the hierarchy
            ``translations``: a dict associating a concept to its human-readable label
    """
    with open(f"{path}/{prior}_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        _get_uniprot_id(gene): set(cids) for gene, cids in annotations.items()
    }

    with open(f"{path}/{prior}_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {cid: set(parents) for cid, parents in hierarchy.items()}

    if os.path.isfile(f"{path}/{prior}_t.json"):
        with open(f"{path}/{prior}_t.json", "r") as file:
            translations = json.load(file)
    else:
        concepts = {p for parents in hierarchy.values() for p in parents}
        concepts = concepts.union(hierarchy.keys())
        translations = {c: c for c in concepts}

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }
