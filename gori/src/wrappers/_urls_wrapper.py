"""Functions used by url_wrapper().

    2025/05/22 @yanisaspic"""

from typing import Any, Optional


def _get_ctyp_url(term: str, header: str, data: Optional[dict[str, Any]]) -> str:
    """Get the descriptive url to a Cell Ontology term.

    ``term`` is an annotation term.
    ``header`` is a url header.
    ``data`` is a dict associating resources (keys) to their contents (values).

    Returns
        A url.
    """
    term = term.replace(":", "_")
    return header + term


def _get_geng_url(term: str, header: str, data: dict[str, Any]) -> str:
    """Get the descriptive url to an HGNC term.

    ``term`` is an annotation term.
    ``header`` is a url header.
    ``data`` is a dict associating resources (keys) to their contents (values).

    Returns
        A url.
    """
    translations = {v: k for k, v in data["HGNC"]["translations"].items()}
    term = translations[term]
    return header + term
