"""Miscalleneous functions used at multiple points during a GORi analysis.

    2025/05/22 @yanisaspic"""

import time
import datetime
import pandas as pd
import networkx as nx
from pypath.utils.go import GOAnnotation
from typing import Any, Callable, Optional
from pypath.utils.mapping import id_from_label
from mlxtend.preprocessing import TransactionEncoder


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

    ``annotations`` is a dict associating priors (keys) to their gene-specific annotations (values).

    Returns
        A binary pd.DataFrame where rows are genes, columns are annotations, and cell values indicate
        if a gene i is annotated by a term j.
    """
    tmp = []
    for prior in annotations.keys():
        te = TransactionEncoder()
        te_ary = te.fit(annotations[prior].values()).transform(
            annotations[prior].values()
        )
        tm = pd.DataFrame(te_ary, columns=te.columns_, index=annotations[prior].keys())
        tmp.append(tm)
    transaction_matrix = pd.concat(tmp, axis=1)
    return transaction_matrix


def _get_generic_ancestors(
    terms: set[str], prior: dict[str, dict[str, Any]], ancestors: set[str] = set()
) -> set[str]:
    """Get the ancestors of multiple terms.
    An ancestor is defined as an element ascending any of the input terms.

    ``terms`` is a set of annotation terms.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".
    ``ancestors`` is a set of terms extended recursively.

    Returns
        A set of ancestor terms.
    """
    H = prior["hierarchy"]
    tmp = [H.get(e, set()) for e in terms]
    parents = {parent for subset in tmp for parent in subset}
    terms = parents.difference(ancestors)
    ancestors = ancestors.union(parents)
    if len(terms) > 0:
        ancestors = _get_generic_ancestors(terms, prior, ancestors)
    return ancestors


def _get_generic_descendants(
    terms: set[str], prior: dict[str, dict[str, Any]]
) -> set[str]:
    """Get the descendants of multiple terms.
    A descendant is defined as an element descending any of the input terms.

    ``terms`` is a set of annotation terms.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".
    ``ancestors`` is a set of terms extended recursively.

    Returns
        A set of descendant terms.
    """
    H = prior["hierarchy"]
    ancestry = {e: _get_generic_ancestors({e}, prior) for e in H.keys()}
    is_descendant = lambda e: len(ancestry.get(e, set()).intersection(terms)) > 0
    descendants = {d for d, ancestors in ancestry.items() if is_descendant(d)}
    return descendants


def _get_prior_ancestors(
    terms: set[str], prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the ancestor annotations of a term from a prior.

    ``terms`` is a set of annotation terms.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of ancestor terms.
    """
    ancestors_wrapper = params["wrappers"]["ancestors_wrapper"]
    if prior in ancestors_wrapper.keys():
        return ancestors_wrapper[prior](terms, data[prior])
    return _get_generic_ancestors(terms, data[prior])


def _get_prior_descendants(
    terms: set[str], prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """A wrapper to get the descendant annotations of a term from a prior.

    ``terms`` is a set of annotation terms.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of descendant terms.
    """
    descendants_wrapper = params["wrappers"]["descendants_wrapper"]
    if prior in descendants_wrapper.keys():
        return descendants_wrapper[prior](terms, data[prior])
    return _get_generic_descendants(terms, data[prior])


def _prune_hierarchy(
    prior: dict[str, dict[str, Any]],
    roots: Optional[set[str]] = None,
    leaves: Optional[set[str]] = None,
) -> dict[str, set[str]]:
    """Prune a prior's hierarchy, i.e. get a new hierarchy with specific roots and leaves.
    If both roots and leaves are indicated, terms associated to both of them are returned.

    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".
    ``roots`` is an optional set of terms to use as roots for the new hierarchy.
    ``leaves`` is an optional set of terms to use as leaves for the new hierarchy.

    Returns
        A dict associating terms (keys) to their parents (values).
    """
    H = prior["hierarchy"]
    ancestors = (
        set(H.keys())
        if leaves is None
        else _get_generic_ancestors(leaves, prior) | leaves
    )
    descendants = (
        set(H.keys())
        if roots is None
        else _get_generic_descendants(roots, prior) | roots
    )
    terms = ancestors.intersection(descendants)
    out = {e: H.get(e, set()) for e in terms}
    return out


def _get_prior_translation(
    term: str,
    prior: str,
    data: dict[str, Any],
    params: dict[str, Any],
    has_prefix: bool = False,
) -> str:
    """Get the human-readable label of an annotation from a prior.

    ``term`` is an annotation term.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.
    ``has_prefix`` is a boolean indicating if the prior-specific prefix should be trimmed (e.g. CTYP).

    Returns
        A human-readable label.
    """
    translate_wrapper = params["wrappers"]["translate_wrapper"]
    _get_generic_translation = lambda t: data[prior]["translations"][t]

    if has_prefix:
        term = term.split(":", 1)[1]
    if prior in translate_wrapper.keys():
        return translate_wrapper[prior](term, data[prior])
    return _get_generic_translation(term)


def _get_generic_inverse_translation(
    label: str, prior: dict[str, dict[str, Any]]
) -> str:
    """Generic function to get the annotation term corresponding to a human-readable label.

    ``label`` is a human-readable label.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        An annotation term
    """
    tmp = {v: k for k, v in prior["translations"].items()}
    return tmp[label]


def _get_prior_inverse_translation(
    label: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> str:
    """Get the annotation term corresponding to a human-readable label from a prior.

    ``label`` is a human-readable label.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        An annotation term
    """
    inverse_translate_wrapper = params["wrappers"]["inverse_translate_wrapper"]
    if prior in inverse_translate_wrapper.keys():
        return inverse_translate_wrapper[prior](label, data[prior])
    return _get_generic_inverse_translation(label, data[prior])


def _get_prior_url(
    term: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> Optional[str]:
    """Get a url to an annotation definition from a prior.

    ``term`` is an annotation term.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A url.
    """
    headers_wrapper = params["wrappers"]["headers_wrapper"]
    urls_wrapper = params["wrappers"]["urls_wrapper"]
    _get_generic_url = lambda t: headers_wrapper[prior] + t

    tmp = term.split(":", 1)[1]
    if prior in headers_wrapper.keys():
        if prior in urls_wrapper:
            return urls_wrapper[prior](tmp, headers_wrapper[prior], data)
        return _get_generic_url(tmp)
    return None
