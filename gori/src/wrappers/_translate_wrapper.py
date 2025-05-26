"""Functions used by translate_wrapper().

    2025/05/22 @yanisaspic"""

from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_translation(term: str, prior: dict[str, Any]) -> str:
    """Get the human-readable label of a CTYP annotation.

    ``term`` is a CTYP annotation.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A human-readable label.
    """
    return prior["ontology"].nodes[term]["name"]


def _get_go_translation(term: str, prior: GOAnnotation) -> str:
    """Get the human-readable label of a GO annotation.

    ``term`` is a GO annotation.
    ``prior`` is a GOAnnotation object.

    Returns
        A human-readable label.
    """
    return prior.get_name(term)


def _get_geng_translation(term: str, prior: dict[str, dict[str, Any]]) -> str:
    """Get the human-readable label of a GENG annotation.

    ``term`` is a GENG annotation.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A human-readable label.
    """
    return term
