"""Functions used by inverse_translate_wrapper().

    2025/05/23 @yanisaspic"""

from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_inverse_translation(label: str, resource: dict[str, Any]) -> str:
    """Get the annotation term of a CellMarker2 human-readable label.

    ``label`` is a CellMarker2 human-readable label.
    ``resource`` is a dict with two keys: "annotations" and "ontology".

    Returns
        An annotation term.
    """
    graph = resource["ontology"]
    tmp = {data["name"]: term for term, data in graph.nodes(data=True)}
    return tmp[label]


def _get_go_inverse_translation(label: str, resource: GOAnnotation) -> str:
    """Get the annotation term of a GO human-readable label.

    ``label`` is a GO human-readable label.
    ``resource`` is a dict with two keys: "annotations" and "ontology".

    Returns
        An annotation term.
    """
    return resource.get_term(label)


def _get_geng_inverse_translation(
    label: str, resource: dict[str, dict[str, Any]]
) -> str:
    """Get the annotation term of a HGNC human-readable label.

    ``label`` is a HGNC human-readable label.
    ``resource`` is a dict with two keys: "annotations" and "ontology".

    Returns
        An annotation term.
    """
    return label
