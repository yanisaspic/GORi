"""Function called to extract relevant annotations in a knowledge base.

    2025/05/12 @yanisaspic"""

import pandas as pd
from typing import Any
from pypath.utils.go import GOAnnotation
from pypath.utils.mapping import label as gene_symbol
from gori.src.utils import (
    _get_generic_ancestors,
    _get_go_ancestors,
    _get_ctyp_ancestors,
    _get_uniprot_id,
)


def _get_generic_annotations(uid: str, prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get the explicit and implicit annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A set of ids.
    """
    exp_annotations = prior["annotations"].get((uid), set())
    imp_annotations = _get_generic_ancestors(exp_annotations, prior["hierarchy"])
    return exp_annotations.union(imp_annotations)


def _get_ctyp_annotations(uid: str, prior: dict[str, Any]) -> set[str]:
    """Get the explicit and implicit CTYP annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CTYP ids.
    """
    exp_annotations = prior["annotations"].get((uid), set())
    imp_annotations = _get_ctyp_ancestors(exp_annotations, prior)
    return exp_annotations.union(imp_annotations)


def _get_dise_annotations(uid: str, prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get the explicit and implicit DISE annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A set of DISE ids.
    """
    tmp = _get_generic_annotations(uid, prior)
    annotations = {a for a in tmp if a in prior["hierarchy"].keys()}
    return annotations


def _get_geng_annotations(uid: str, prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get the explicit and implicit GENG annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A set of GENG ids.
    """
    exp_annotations = prior["annotations"].get((uid), set())
    exp_annotations = {a.mainclass for a in exp_annotations}
    imp_annotations = _get_generic_ancestors(exp_annotations, prior["hierarchy"])
    return exp_annotations.union(imp_annotations)


def _get_go_annotations(uid: str, prior: GOAnnotation, aspect: str) -> set[str]:
    """Get the explicit and implicit GO annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a GOAnnotation object.
    ``aspect`` is one of "p" (BIOP), "c" (CELC) or "f" (MOLF).

    Returns
        A set of GO IDs.
    """
    exp_annotations = prior.get_annot(uid, aspect)
    imp_annotations = _get_go_ancestors(exp_annotations, prior)
    return exp_annotations.union(imp_annotations)


def _get_phen_annotations(uid: str, prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get the explicit and implicit PHEN annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A set of PHEN ids.
    """
    exp_annotations = prior["annotations"].get((uid), set())
    exp_annotations = {a.hpo_id for a in exp_annotations}
    exp_annotations = {a for a in exp_annotations if a in prior["hierarchy"].keys()}
    imp_annotations = _get_generic_ancestors(exp_annotations, prior["hierarchy"])
    return exp_annotations.union(imp_annotations)


def _get_prior_annotations(uid: str, prior: str, data: dict[str, Any]) -> set[str]:
    """A wrapper to get the annotations of a gene from a prior.

    ``uid`` is a UniProtID.
    ``prior`` is the name of a knowledge base to use.
    ``data`` is a dict associating priors (keys) to their contents (values).

    Returns
        A set of ids.
    """
    if prior == "BIOP":
        annotations = _get_go_annotations(uid, data[prior], "p")
    elif prior == "CELC":
        annotations = _get_go_annotations(uid, data[prior], "c")
    elif prior == "CTYP":
        annotations = _get_ctyp_annotations(uid, data[prior])
    elif prior == "DISE":
        annotations = _get_dise_annotations(uid, data[prior])
    elif prior == "GENG":
        annotations = _get_geng_annotations(uid, data[prior])
    elif prior == "PHEN":
        annotations = _get_phen_annotations(uid, data[prior])
    elif prior == "MOLF":
        annotations = _get_go_annotations(uid, data[prior], "f")
    else:
        annotations = _get_generic_annotations(uid, data[prior])
    out = {f"{prior}:{a}" for a in annotations}
    return out


def _get_annotations_counter(
    annotations: dict[str, dict[str, set[str]]]
) -> pd.DataFrame:
    """Get the annotations of a geneset across multiple priors.

    ``annotations`` is a dict associating prior (keys) to their relevant annotations (values).

    Returns
        A pandas DataFrame associating each prior to its number of annotations.
    """
    annotations_counter = [
        (gene, prior, len(annotations0))
        for prior, tmp in annotations.items()
        for gene, annotations0 in tmp.items()
    ]
    df = pd.DataFrame(annotations_counter, columns=["gene", "prior", "n_annotations"])
    df = df.pivot(index="gene", columns="prior", values="n_annotations")
    return df


def _get_annotations(
    geneset: set[str], data: dict[str, Any], use_gene_symbol: bool
) -> dict[str, dict[str, set[str]]]:
    """Get the annotations of a geneset across multiple priors.

    ``geneset`` is a set of gene symbols or UniProtIDs.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``use_gene_symbol`` is a boolean indicating if the gene symbols should be used.

    Returns
        A dict associating priors (keys) to their gene-specific annotations.
    """
    f = lambda gene: gene_symbol(gene) if use_gene_symbol else gene
    annotations = {
        prior: {
            f(gene): _get_prior_annotations(_get_uniprot_id(gene), prior, data)
            for gene in geneset
        }
        for prior in data.keys()
    }
    return annotations
