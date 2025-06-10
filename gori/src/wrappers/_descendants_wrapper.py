"""Functions used by descendants_wrapper().

    2025/05/22 @yanisaspic"""

import networkx as nx
from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_descendants(
    terms: set[str], resource: dict[str, dict[str, Any]]
) -> set[str]:
    """Get every CellMarker2 descendant terms.
    A term is descendant if at least one of the input CellMarker2 ontology term is its descendant.

    ``terms`` is a set of CellMarker2 ontology terms.
    ``resource`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CellMarker2 ids.
    """
    tmp = [nx.descendants(resource["ontology"], t) for t in terms]
    descendants = {descendant for descendants0 in tmp for descendant in descendants0}
    return descendants


def _get_go_descendants(terms: set[str], resource: GOAnnotation) -> set[str]:
    """Get every GO descendant terms.
    A term is descendant if at least one of the input GO ontology term is its ancestor.

    ``terms`` is a set of GO ontology terms.
    ``resource`` is a GOAnnotation object.

    Returns
        A set of GO ids.
    """
    _get_descendants0 = lambda t: resource.ontology.get_all_descendants(
        t, relations="is_a"
    ).difference({t})
    tmp = [_get_descendants0(t) for t in terms]
    descendants = {descendant for descendants0 in tmp for descendant in descendants0}
    return descendants
