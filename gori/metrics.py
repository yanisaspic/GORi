"""Functions called to measure the information content of a term.

    2025/05/13 @yanisaspic"""

import pandas as pd
from math import log2
from typing import Any
from gori.params import get_parameters
from gori.src.utils import (
    _get_prior_ancestors,
    _get_prior_descendants,
    _get_prior_translation,
)


def _get_generic_terms(prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get every term in a knowledge base.

    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A set of terms.
    """
    return set(prior["translations"].keys())


def _get_prior_terms(
    prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get every term in a knowledge base.

    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of terms.
    """
    terms_wrapper = params["wrappers"]["terms_wrapper"]
    if prior in terms_wrapper.keys():
        return terms_wrapper[prior](data[prior])
    return _get_generic_terms(data[prior])


def _get_prior_limits(
    prior: str, data: dict[str, Any], params: dict[str, Any]
) -> dict[str, set[str]]:
    """Get every root and leaf in a knowledge base.

    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of terms.
    """
    terms = _get_prior_terms(prior, data, params)
    roots = terms.difference(_get_prior_descendants(terms, prior, data, params))
    leaves = terms.difference(_get_prior_ancestors(terms, prior, data, params))
    return {"roots": roots, "leaves": leaves}


def _get_prior_iic0(
    term: str,
    limits: dict[str, set[str]],
    prior: str,
    data: dict[str, Any],
    params: dict[str, Any],
) -> float:
    """Get the Intrinsic Information Content (IIC) of a term.

    IIC(t) = -log[ ((leaves(t) / subsumers(t)) + 1) / (max_leaves + 1) ], where:
        `IIC(t)` is the intrinsic information content of the term t,
        `leaves(t)` is the number of leaf descendant terms of t,
        `subsumers(t)` is the number of ancestor terms of t, and t,
        `max_leaves` is the number of leaf terms in the GO ontology.

    Sánchez, David, et al.
    "Ontology-based information content computation."
    Knowledge-based systems 24.2 (2011): 297-303.
    https://doi.org/10.1016/j.knosys.2010.10.001

    ``term`` is an annotation term.
    ``n_limits`` is a dict with the ``roots`` and ``leaves`` count of the prior.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A float ranging from 0 to 1.
    """
    subsumers_t = _get_prior_ancestors({term}, prior, data, params) | {term}
    leaves_t = _get_prior_descendants({term}, prior, data, params).intersection(
        limits["leaves"]
    )
    tmp = ((len(leaves_t) / len(subsumers_t)) + 1) / (len(limits["leaves"]) + 1)
    iic = -log2(tmp)
    return iic


def _get_prior_iics(
    prior: str, data: dict[str, Any], params: dict[str, Any]
) -> pd.DataFrame:
    """Get the Intrinsic Information Content (IIC) of every term in a knowledge base,
    cf. _get_prior_iic0().

    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.Series associating terms to their iics.
    """
    terms = _get_prior_terms(prior, data, params)
    limits = _get_prior_limits(prior, data, params)
    tmp = {t: _get_prior_iic0(t, limits, prior, data, params) for t in terms}
    out = pd.DataFrame.from_dict(tmp, orient="index", columns=["iic"])
    return out


def _get_iics(
    data: dict[str, Any], params: dict[str, Any] = get_parameters()
) -> pd.DataFrame:
    """Get a table associating each annotation term to its raw and normalized IICs, its prior and
    its human-readable label.

    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame where rows are terms, and four columns: `label`, `iic` and `prior`.
    """
    tmp = []
    for p in data.keys():
        p_out = _get_prior_iics(p, data, params)
        p_out["n_iic"] = p_out.iic / p_out.iic.max()
        p_out["label"] = [
            _get_prior_translation(t, p, data, params, has_prefix=False)
            for t in p_out.index
        ]
        p_out["prior"] = p
        tmp.append(p_out)
    iic = pd.concat(tmp)
    return iic


def get_prior_eics(
    corpus: dict[str, set[str]],
    prior: str,
    data: dict[str, Any],
    params: dict[str, Any] = get_parameters(),
) -> dict[str, float]:
    """Get the Extrinsic Information Content (EIC) of every term in a corpus. It ranges from 0 to +Inf.

    EIC(t) = -log[ p(t) ], where:
        `EIC(t)` is the extrinsic information content of the term t,
        `p(t)` is the frequency of a term t in the corpus.

    Resnik, Philip.
    "Using information content to evaluate semantic similarity in a taxonomy."
    arXiv preprint cmp-lg/9511007 (1995).
    See: https://doi.org/10.48550/arXiv.cmp-lg/9511007

    ``corpus`` is a dict associating some group (keys) to their annotation terms (values).
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A dict associating terms (keys) to their EIC (values).
    """
