"""Functions used by terms_wrapper().

    2025/05/23 @yanisaspic"""

from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_terms(prior: dict[str, Any]) -> set[str]:
    """Get every term in the CTYP knowledge base.

    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CL ids.
    """
    return set(prior["ontology"].nodes)


def _get_go_terms(prior: GOAnnotation, aspect: str) -> set[str]:
    """Get every term in a GO knowledge base.

    ``prior`` is a GO annotation.
    ``aspect`` is one of `P` (BIOP), `C` (CELC) or `F` (MOLF).

    Returns
        A set of GO ids.
    """
    return prior.ontology.all_from_aspect(aspect)


def _get_biop_terms(prior: GOAnnotation) -> set[str]:
    """Get every term in the BIOP knowledge base.

    ``prior`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return _get_go_terms(prior, "P")


def _get_celc_terms(prior: GOAnnotation) -> set[str]:
    """Get every term in the CELC knowledge base.

    ``prior`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return _get_go_terms(prior, "C")


def _get_geng_terms(prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get every term in the MOLF knowledge base.

    ``prior`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return set(prior["translations"].values())


def _get_molf_terms(prior: GOAnnotation) -> set[str]:
    """Get every term in the MOLF knowledge base.

    ``prior`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return _get_go_terms(prior, "F")
