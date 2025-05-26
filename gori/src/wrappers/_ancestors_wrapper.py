"""Functions used by ancestors_wrapper().

    2025/05/22 @yanisaspic"""

import networkx as nx
from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_ancestors(terms: set[str], prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get every CTYP ancestor terms.
    A term is ancestor if at least one of the input CTYP ontology term descends from it.

    ``terms`` is a set of CTYP ontology terms.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CTYP ids.
    """
    tmp = [nx.ancestors(prior["ontology"], t) for t in terms]
    ancestors = {ancestor for ancestors0 in tmp for ancestor in ancestors0}
    return ancestors


def _get_go_ancestors(terms: set[str], prior: GOAnnotation) -> set[str]:
    """Get every GO ancestor terms.
    A term is ancestor if at least one of the input GO ontology term descends from it.

    ``terms`` is a set of GO ontology terms.
    ``prior`` is a GOAnnotation object.

    Returns
        A set of GO ids.
    """
    _get_ancestors0 = lambda t: prior.ontology.get_all_ancestors(
        t, relations="is_a"
    ).difference({t})
    tmp = [_get_ancestors0(t) for t in terms]
    ancestors = {ancestor for ancestors0 in tmp for ancestor in ancestors0}
    return ancestors
