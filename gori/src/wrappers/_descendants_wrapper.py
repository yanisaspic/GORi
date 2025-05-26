"""Functions used by descendants_wrapper().

    2025/05/22 @yanisaspic"""

import networkx as nx
from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_descendants(
    terms: set[str], prior: dict[str, dict[str, Any]]
) -> set[str]:
    """Get every CTYP descendant terms.
    A term is descendant if at least one of the input CTYP ontology term is its descendant.

    ``terms`` is a set of CTYP ontology terms.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CTYP ids.
    """
    tmp = [nx.descendants(prior["ontology"], t) for t in terms]
    descendants = {descendant for descendants0 in tmp for descendant in descendants0}
    return descendants


def _get_go_descendants(terms: set[str], prior: GOAnnotation) -> set[str]:
    """Get every GO descendant terms.
    A term is descendant if at least one of the input GO ontology term is its ancestor.

    ``terms`` is a set of GO ontology terms.
    ``prior`` is a GOAnnotation object.

    Returns
        A set of GO ids.
    """
    _get_descendants0 = lambda t: prior.ontology.get_all_descendants(
        t, relations="is_a"
    ).difference({t})
    tmp = [_get_descendants0(t) for t in terms]
    descendants = {descendant for descendants0 in tmp for descendant in descendants0}
    return descendants
