"""Functions used by annotations_wrapper().

    2025/05/22 @yanisaspic"""

from typing import Any
from gori.src.utils import _get_resource_ancestors


def _get_ctyp_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit CellMarker2 or CellTaxonomy annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of Cell Ontology ids.
    """
    exp_annotations = data[resource]["annotations"].get((uid), set())
    imp_annotations = _get_resource_ancestors(exp_annotations, resource, data, params)
    return exp_annotations.union(imp_annotations)


def _get_dise_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit MeSH annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of MeSH ids.
    """
    exp_annotations = data[resource]["annotations"].get((uid), set())
    imp_annotations = _get_resource_ancestors(exp_annotations, resource, data, params)
    tmp = exp_annotations.union(imp_annotations)
    annotations = {a for a in tmp if a in data[resource]["hierarchy"].keys()}
    return annotations


def _get_geng_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit HGNC annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of HGNC ids.
    """
    exp_annotations = data[resource]["annotations"].get((uid), set())
    exp_annotations = {a.mainclass for a in exp_annotations}
    imp_annotations = _get_resource_ancestors(exp_annotations, resource, data, params)
    return exp_annotations.union(imp_annotations)


def _get_go_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any], aspect: str
) -> set[str]:
    """Get the explicit and implicit GO annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.
    ``aspect`` is one of `p` (GO_BP), `c` (GO_CC) or `f` (GO_MF).

    Returns
        A set of GO IDs.
    """
    exp_annotations = data[resource].get_annot(uid, aspect)
    imp_annotations = _get_resource_ancestors(exp_annotations, resource, data, params)
    return exp_annotations.union(imp_annotations)


def _get_biop_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit GO_BP annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GO IDs.
    """
    return _get_go_annotations(uid, resource, data, params, "p")


def _get_celc_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit GO_CC annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GO IDs.
    """
    return _get_go_annotations(uid, resource, data, params, "c")


def _get_molf_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit GO_MF annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GO IDs.
    """
    return _get_go_annotations(uid, resource, data, params, "f")


def _get_phen_annotations(
    uid: str, resource: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit HPO annotations of a gene.

    ``uid`` is a UniProtID.
    ``resource`` is a resource label.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of HPO ids.
    """
    exp_annotations = data[resource]["annotations"].get((uid), set())
    exp_annotations = {a.hpo_id for a in exp_annotations}
    exp_annotations = {
        a for a in exp_annotations if a in data[resource]["hierarchy"].keys()
    }
    imp_annotations = _get_resource_ancestors(exp_annotations, resource, data, params)
    return exp_annotations.union(imp_annotations)
