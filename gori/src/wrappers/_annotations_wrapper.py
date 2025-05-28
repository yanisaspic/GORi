"""Functions used by annotations_wrapper().

    2025/05/22 @yanisaspic"""

from typing import Any
from gori.src.utils import _get_prior_ancestors


def _get_ctyp_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit CTYP annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of CTYP ids.
    """
    exp_annotations = data[prior]["annotations"].get((uid), set())
    imp_annotations = _get_prior_ancestors(exp_annotations, prior, data, params)
    return exp_annotations.union(imp_annotations)


def _get_dise_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit DISE annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of DISE ids.
    """
    exp_annotations = data[prior]["annotations"].get((uid), set())
    imp_annotations = _get_prior_ancestors(exp_annotations, prior, data, params)
    tmp = exp_annotations.union(imp_annotations)
    annotations = {a for a in tmp if a in data[prior]["hierarchy"].keys()}
    return annotations


def _get_geng_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit GENG annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GENG ids.
    """
    exp_annotations = data[prior]["annotations"].get((uid), set())
    exp_annotations = {a.mainclass for a in exp_annotations}
    imp_annotations = _get_prior_ancestors(exp_annotations, prior, data, params)
    return exp_annotations.union(imp_annotations)


def _get_go_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any], aspect: str
) -> set[str]:
    """Get the explicit and implicit GO annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.
    ``aspect`` is one of `p` (BIOP), `c` (CELC) or `f` (MOLF).

    Returns
        A set of GO IDs.
    """
    exp_annotations = data[prior].get_annot(uid, aspect)
    imp_annotations = _get_prior_ancestors(exp_annotations, prior, data, params)
    return exp_annotations.union(imp_annotations)


def _get_biop_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit BIOP annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GO IDs.
    """
    return _get_go_annotations(uid, prior, data, params, "p")


def _get_celc_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit CELC annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GO IDs.
    """
    return _get_go_annotations(uid, prior, data, params, "c")


def _get_molf_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit MOLF annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of GO IDs.
    """
    return _get_go_annotations(uid, prior, data, params, "f")


def _get_phen_annotations(
    uid: str, prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get the explicit and implicit PHEN annotations of a gene.

    ``uid`` is a UniProtID.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of PHEN ids.
    """
    exp_annotations = data[prior]["annotations"].get((uid), set())
    exp_annotations = {a.hpo_id for a in exp_annotations}
    exp_annotations = {
        a for a in exp_annotations if a in data[prior]["hierarchy"].keys()
    }
    imp_annotations = _get_prior_ancestors(exp_annotations, prior, data, params)
    return exp_annotations.union(imp_annotations)
