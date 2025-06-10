"""Functions called to extract relevant annotations in a knowledge base.

    2025/05/12 @yanisaspic"""

import pandas as pd
from typing import Any
from gori.src.utils import _get_uniprot_id, _get_resource_ancestors, _get_gene_symbol


def _get_generic_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of ids.
    """
    exp_annotations = data[resource]["annotations"].get((uid), set())
    imp_annotations = _get_resource_ancestors(exp_annotations, resource, data, params)
    return exp_annotations.union(imp_annotations)


def _get_resource_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """A wrapper to get the annotations of a gene from a resource.

    ``uid`` is a UniProtID.
    ``resource`` is the name of a knowledge base to use.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of ids.
    """
    annotations_wrapper = params["wrappers"]["annotations_wrapper"]
    if resource in annotations_wrapper.keys():
        annotations = annotations_wrapper[resource](uid, resource, data, params)
    else:
        annotations = _get_generic_annotations(uid, resource, data, params)
    out = {f"{resource}:{a}" for a in annotations}
    return out


def _get_annotations(
    geneset: set[str],
    data: dict[str, Any],
    params: dict[str, Any],
) -> dict[str, dict[str, set[str]]]:
    """Get the annotations of a geneset across multiple resources.

    ``geneset`` is a set of gene symbols or UniProtIDs.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A dict associating resources (keys) to their gene-specific annotations.
    """
    f = lambda gene: _get_gene_symbol(gene) if params["use_gene_symbol"] else gene
    annotations = {
        resource: {
            f(gene): _get_resource_annotations(
                _get_uniprot_id(gene), resource, data, params
            )
            for gene in geneset
        }
        for resource in data.keys()
    }
    return annotations


def _get_annotations_counter(
    annotations: dict[str, dict[str, set[str]]]
) -> pd.DataFrame:
    """Get the annotations of a geneset across multiple resources.

    ``annotations`` is a dict associating resources (keys) to their relevant annotations (values).

    Returns
        A pandas DataFrame associating each resource to its number of annotations.
    """
    annotations_counter = [
        (gene, resource, len(annotations0))
        for resource, tmp in annotations.items()
        for gene, annotations0 in tmp.items()
    ]
    df = pd.DataFrame(
        annotations_counter, columns=["gene", "resource", "n_annotations"]
    )
    df = df.pivot(index="gene", columns="resource", values="n_annotations")
    return df
