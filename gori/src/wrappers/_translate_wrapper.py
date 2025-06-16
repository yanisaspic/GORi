"""Functions used by translate_wrapper().

    2025/05/22 @yanisaspic"""

from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_translation(term: str, resource: dict[str, Any]) -> str:
    """Get the human-readable label of a Cell Ontology annotation.

    ``term`` is a Cell Ontology annotation.
    ``resource`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A human-readable label.
    """
    return resource["ontology"].nodes[term]["name"]


def _get_go_translation(term: str, resource: GOAnnotation) -> str:
    """Get the human-readable label of a GO annotation.

    ``term`` is a GO annotation.
    ``resource`` is a GOAnnotation object.

    Returns
        A human-readable label.
    """
    return resource.get_name(term)


def _get_geng_translation(term: str, resource: dict[str, dict[str, Any]]) -> str:
    """Get the human-readable label of a HGNC annotation.

    ``term`` is a HGNC annotation.
    ``resource`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A human-readable label.
    """
    return term
