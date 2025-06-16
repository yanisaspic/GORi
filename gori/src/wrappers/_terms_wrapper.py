"""Functions used by terms_wrapper().

    2025/05/23 @yanisaspic"""

from typing import Any
from pypath.utils.go import GOAnnotation


def _get_ctyp_terms(resource: dict[str, Any]) -> set[str]:
    """Get every term in the Cell Ontology knowledge base.

    ``resource`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CL ids.
    """
    return set(resource["ontology"].nodes)


def _get_go_terms(resource: GOAnnotation, aspect: str) -> set[str]:
    """Get every term in a GO knowledge base.

    ``resource`` is a GO annotation.
    ``aspect`` is one of `P` (GO_BP), `C` (GO_CC) or `F` (GO_MF).

    Returns
        A set of GO ids.
    """
    return resource.ontology.all_from_aspect(aspect)


def _get_biop_terms(resource: GOAnnotation) -> set[str]:
    """Get every term in the GO_BP knowledge base.

    ``resource`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return _get_go_terms(resource, "P")


def _get_celc_terms(resource: GOAnnotation) -> set[str]:
    """Get every term in the GO_CC knowledge base.

    ``resource`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return _get_go_terms(resource, "C")


def _get_geng_terms(resource: dict[str, dict[str, Any]]) -> set[str]:
    """Get every term in the HGNC knowledge base.

    ``resource`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return set(resource["translations"].values())


def _get_molf_terms(resource: GOAnnotation) -> set[str]:
    """Get every term in the GO_MF knowledge base.

    ``resource`` is a GO annotation.

    Returns
        A set of GO ids.
    """
    return _get_go_terms(resource, "F")
