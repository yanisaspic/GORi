"""Functions used by inverse_translate_wrapper().

    2025/05/23 @yanisaspic"""

from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_inverse_translation(label: str, prior: dict[str, Any]) -> str:
    """Get the annotation term of a CTYP human-readable label.

    ``label`` is a CTYP human-readable label.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        An annotation term.
    """
    graph = prior["ontology"]
    tmp = {data["name"]: term for term, data in graph.nodes(data=True)}
    return tmp[label]


def _get_go_inverse_translation(label: str, prior: GOAnnotation) -> str:
    """Get the annotation term of a GO human-readable label.

    ``label`` is a GO human-readable label.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        An annotation term.
    """
    return prior.get_term(label)


def _get_geng_inverse_translation(label: str, prior: dict[str, dict[str, Any]]) -> str:
    """Get the annotation term of a GENG human-readable label.

    ``label`` is a GENG human-readable label.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        An annotation term.
    """
    return label
