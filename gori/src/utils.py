"""Miscalleneous functions used at multiple points during a GORi analysis.

    2025/05/22 @yanisaspic"""

import time
import datetime
import pandas as pd
from typing import Any, Optional
from pypath.utils.mapping import label, id_from_label
from mlxtend.preprocessing import TransactionEncoder


def _get_gene_symbol(gene: str) -> str:
    """Get a gene symbol.

    ``gene`` is a gene symbol (or a UniProtID).

    Returns
        A gene symbol.
    """
    if gene.startswith(
        "MI"
    ):  # prevents errors raised from pypath querying mir-base wrongfully.
        return gene
    return label(gene)


def _get_timestamp() -> str:
    """Get a formatted timestamp.

    Returns
        A string with the date and the time.
    """
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime(
        "%Y/%m/%d_%H:%M:%S"
    )
    return timestamp


def _get_uniprot_id(gene: str) -> str:
    """Get a UniProtID.

    ``gene`` is a gene symbol (or a UniProtID).

    Returns
        A UniProtID.
    """
    uids = id_from_label(gene)  # use id_from_label to get a stable output.
    tmp = sorted(list(uids))
    uid = tmp[0]
    if len(uids) > 1:
        print(
            f"Multiple UniProtIDs found for gene {gene}: {uids}. Selecting {uid} (alphanumeric order)."
        )
    return uid


def _get_transaction_matrix(
    annotations: dict[str, dict[str, set[str]]]
) -> pd.DataFrame:
    """Get a transaction matrix from a set of annotations.

    ``annotations`` is a dict associating resources (keys) to their gene-specific annotations (values).

    Returns
        A binary pd.DataFrame where rows are genes, columns are annotations, and cell values indicate
        if a gene i is annotated by a term j.
    """
    tmp = []
    for resource in annotations.keys():
        te = TransactionEncoder()
        te_ary = te.fit(annotations[resource].values()).transform(
            annotations[resource].values()
        )
        tm = pd.DataFrame(
            te_ary, columns=te.columns_, index=annotations[resource].keys()
        )
        tmp.append(tm)
    transaction_matrix = pd.concat(tmp, axis=1)
    return transaction_matrix


def _get_generic_ancestors(
    terms: set[str], resource: dict[str, dict[str, Any]], ancestors: set[str] = set()
) -> set[str]:
    """Get the ancestors of multiple terms.
    An ancestor is defined as an element ascending any of the input terms.

    ``terms`` is a set of annotation terms.
    ``resource`` is a dict with three keys: "annotations", "hierarchy" and "translations".
    ``ancestors`` is a set of terms extended recursively.

    Returns
        A set of ancestor terms.
    """
    H = resource["hierarchy"]
    tmp = [H.get(e, set()) for e in terms]
    parents = {parent for subset in tmp for parent in subset}
    terms = parents.difference(ancestors)
    ancestors = ancestors.union(parents)
    if len(terms) > 0:
        ancestors = _get_generic_ancestors(terms, resource, ancestors)
    return ancestors


def _get_generic_descendants(
    terms: set[str], resource: dict[str, dict[str, Any]]
) -> set[str]:
    """Get the descendants of multiple terms.
    A descendant is defined as an element descending any of the input terms.

    ``terms`` is a set of annotation terms.
    ``resource`` is a dict with three keys: "annotations", "hierarchy" and "translations".
    ``ancestors`` is a set of terms extended recursively.

    Returns
        A set of descendant terms.
    """
    H = resource["hierarchy"]
    ancestry = {e: _get_generic_ancestors({e}, resource) for e in H.keys()}
    is_descendant = lambda e: len(ancestry.get(e, set()).intersection(terms)) > 0
    descendants = {d for d, ancestors in ancestry.items() if is_descendant(d)}
    return descendants


def _get_resource_ancestors(
    terms: set[str], resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the ancestor annotations of a term from a resource.

    ``terms`` is a set of annotation terms.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of ancestor terms.
    """
    ancestors_wrapper = params["wrappers"]["ancestors_wrapper"]
    if resource in ancestors_wrapper.keys():
        return ancestors_wrapper[resource](terms, data[resource])
    return _get_generic_ancestors(terms, data[resource])


def _get_resource_descendants(
    terms: set[str], resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """A wrapper to get the descendant annotations of a term from a resource.

    ``terms`` is a set of annotation terms.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of descendant terms.
    """
    descendants_wrapper = params["wrappers"]["descendants_wrapper"]
    if resource in descendants_wrapper.keys():
        return descendants_wrapper[resource](terms, data[resource])
    return _get_generic_descendants(terms, data[resource])


def _prune_hierarchy(
    resource: dict[str, dict[str, Any]],
    roots: Optional[set[str]] = None,
    leaves: Optional[set[str]] = None,
) -> dict[str, set[str]]:
    """Prune a resource's hierarchy, i.e. get a new hierarchy with specific roots and leaves.
    If both roots and leaves are indicated, terms associated to both of them are returned.

    ``resource`` is a dict with three keys: "annotations", "hierarchy" and "translations".
    ``roots`` is an optional set of terms to use as roots for the new hierarchy.
    ``leaves`` is an optional set of terms to use as leaves for the new hierarchy.

    Returns
        A dict associating terms (keys) to their parents (values).
    """
    H = resource["hierarchy"]
    ancestors = (
        set(H.keys())
        if leaves is None
        else _get_generic_ancestors(leaves, resource) | leaves
    )
    descendants = (
        set(H.keys())
        if roots is None
        else _get_generic_descendants(roots, resource) | roots
    )
    terms = ancestors.intersection(descendants)
    out = {e: H.get(e, set()) for e in terms}
    return out


def _get_resource_translation(
    term: str,
    resource: str,
    data: dict[str, Any],
    params: dict[str, Any],
    has_prefix: bool = False,
) -> str:
    """Get the human-readable label of an annotation from a resource.

    ``term`` is an annotation term.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.
    ``has_prefix`` is a boolean indicating if the resource-specific prefix should be trimmed (e.g. CellMarker2).

    Returns
        A human-readable label.
    """
    translate_wrapper = params["wrappers"]["translate_wrapper"]
    _get_generic_translation = lambda t: data[resource]["translations"][t]

    if has_prefix:
        term = term.split(":", 1)[1]
    if resource in translate_wrapper.keys():
        return translate_wrapper[resource](term, data[resource])
    return _get_generic_translation(term)


def _get_generic_inverse_translation(
    label: str, resource: dict[str, dict[str, Any]]
) -> str:
    """Generic function to get the annotation term corresponding to a human-readable label.

    ``label`` is a human-readable label.
    ``resource`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        An annotation term
    """
    tmp = {v: k for k, v in resource["translations"].items()}
    return tmp[label]


def _get_resource_inverse_translation(
    label: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> str:
    """Get the annotation term corresponding to a human-readable label from a resource.

    ``label`` is a human-readable label.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        An annotation term
    """
    inverse_translate_wrapper = params["wrappers"]["inverse_translate_wrapper"]
    if resource in inverse_translate_wrapper.keys():
        return inverse_translate_wrapper[resource](label, data[resource])
    return _get_generic_inverse_translation(label, data[resource])


def _get_resource_url(
    term: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> Optional[str]:
    """Get a url to an annotation definition from a resource.

    ``term`` is an annotation term.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A url.
    """
    headers_wrapper = params["wrappers"]["headers_wrapper"]
    urls_wrapper = params["wrappers"]["urls_wrapper"]
    _get_generic_url = lambda t: headers_wrapper[resource] + t

    tmp = term.split(":", 1)[1]
    if resource in headers_wrapper.keys():
        if resource in urls_wrapper:
            return urls_wrapper[resource](tmp, headers_wrapper[resource], data)
        return _get_generic_url(tmp)
    return None


def _get_generic_terms(resource: dict[str, dict[str, Any]]) -> set[str]:
    """Get every term in a knowledge base.

    ``resource`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A set of terms.
    """
    return set(resource["translations"].keys())


def _get_resource_terms(
    resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get every term in a knowledge base.

    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of terms.
    """
    terms_wrapper = params["wrappers"]["terms_wrapper"]
    if resource in terms_wrapper.keys():
        return terms_wrapper[resource](data[resource])
    return _get_generic_terms(data[resource])


def _get_resource_boundaries(
    resource: str, data: dict[str, Any], params: dict[str, Any]
) -> dict[str, set[str]]:
    """Get the boundaries of a resource.

    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of boundary terms.
    """
    terms = _get_resource_terms(resource, data, params)
    roots_wrapper = params["wrappers"]["roots_wrapper"]
    if resource in roots_wrapper.keys():
        roots = roots_wrapper[resource]()
    else:
        roots = terms.difference(
            _get_resource_descendants(terms, resource, data, params)
        )

    tmp = terms.difference(_get_resource_ancestors(terms, resource, data, params))
    leaves = {
        t
        for t in tmp
        if _get_resource_ancestors({t}, resource, data, params).intersection(roots)
        != set()
    }  # filter out terms that do not share an is_a path to a root
    return {"roots": roots, "leaves": leaves}
