"""Functions called to load the knowledge base resources for the GORi analysis.

    2025/06/10 @yanisaspic"""

import os
import json
import pandas as pd
from typing import Any
from gori.params import get_parameters
from gori.src.utils import _get_uniprot_id


def load_resources(
    resources: set[str],
    path: str = "./resources",
    params: dict[str, Any] = get_parameters(),
) -> dict[str, Any]:
    """Load a subset of knowledge bases readily available in the GORi package.

    ``resources`` is a set of strings containing the names of the knowledge bases to load, including:
        `GO_BP` (biological processes),
        `GO_CC` (cellular components),
        `CellMarker2` (cell types with CellMarker 2.0 annotations),
        `CellTaxonomy` (cell types with CellTaxonomy annotations),
        `MeSH` (diseases),
        `HGNC` (gene groups),
        `GO_MF` (molecular functions),
        `Reactome` (pathways) and
        `HPO` (phenotypes).
    ``path`` is the path to the JSON files containing the knowledge bases.
    ``params`` is a dict of parameters.

    Returns
        A dict associating knowledge base names (keys) to their contents (values).
    """
    load_wrapper = params["wrappers"]["load_wrapper"]
    if any(p in resources for p in ["GO_BP", "GO_CC", "GO_MF"]):
        go = load_wrapper["GO_BP"](
            path
        )  # outputs the same GOAnnotation as GO_CC and GO_MF

    knowledge_bases = {}
    for p in resources:
        if p not in load_wrapper.keys():
            raise ValueError(
                f"Invalid resource: {p}. Valid resources are: {set(load_wrapper)}"
            )
        elif p in ["GO_BP", "GO_CC", "GO_MF"]:
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


def load_local(resource: str, path: str) -> dict[str, dict[str, Any]]:
    """Load a local knowledge base.

    ``resource`` is the label of the knowledge base to add.
    ``path`` is the path to the JSON files containing the annotations, the hierarchy
    and the translations of the knowledge base.

    Returns
        A dict with three keys:
            ``annotations``: a dict associating a Uniprot ID to its related concepts
            ``ontology``: a dict associating a concept to its parent concepts in the hierarchy
            ``translations``: a dict associating a concept to its human-readable label
    """
    with open(f"{path}/{resource}_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        _get_uniprot_id(gene): set(cids) for gene, cids in annotations.items()
    }

    with open(f"{path}/{resource}_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {cid: set(parents) for cid, parents in hierarchy.items()}

    if os.path.isfile(f"{path}/{resource}_t.json"):
        with open(f"{path}/{resource}_t.json", "r") as file:
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
