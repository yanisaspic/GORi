"""Functions called to extract relevant annotations in a knowledge base.

    2025/05/12 @yanisaspic"""

import pandas as pd
from typing import Any
from pypath.utils.go import GOAnnotation
from pypath.utils.mapping import label as gene_symbol
from gori.src.utils import _get_uniprot_id, _get_prior_ancestors


def _get_generic_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of ids.
    """
    exp_annotations = data[prior]["annotations"].get((uid), set())
    imp_annotations = _get_prior_ancestors(exp_annotations, prior, data, params)
    return exp_annotations.union(imp_annotations)


def _get_prior_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """A wrapper to get the annotations of a gene from a prior.

    ``uid`` is a UniProtID.
    ``prior`` is the name of a knowledge base to use.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of ids.
    """
    annotations_wrapper = params["wrappers"]["annotations_wrapper"]
    if prior in annotations_wrapper.keys():
        annotations = annotations_wrapper[prior](uid, prior, data, params)
    else:
        annotations = _get_generic_annotations(uid, prior, data, params)
    out = {f"{prior}:{a}" for a in annotations}
    return out


def _get_annotations(
    geneset: set[str],
    data: dict[str, Any],
    params: dict[str, Any],
) -> dict[str, dict[str, set[str]]]:
    """Get the annotations of a geneset across multiple priors.

    ``geneset`` is a set of gene symbols or UniProtIDs.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A dict associating priors (keys) to their gene-specific annotations.
    """
    f = lambda gene: gene_symbol(gene) if params["use_gene_symbol"] else gene
    annotations = {
        prior: {
            f(gene): _get_prior_annotations(_get_uniprot_id(gene), prior, data, params)
            for gene in geneset
        }
        for prior in data.keys()
    }
    return annotations


def _get_annotations_counter(
    annotations: dict[str, dict[str, set[str]]]
) -> pd.DataFrame:
    """Get the annotations of a geneset across multiple priors.

    ``annotations`` is a dict associating priors (keys) to their relevant annotations (values).

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
