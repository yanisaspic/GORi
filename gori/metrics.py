"""Functions called to measure the information content of a term.

    2025/05/13 @yanisaspic"""

import numpy as np
import pandas as pd
from math import log2
from typing import Any
from gori.params import get_parameters
from gori.src.utils import (
    _get_prior_ancestors,
    _get_prior_boundaries,
    _get_prior_descendants,
    _get_prior_inverse_translation,
    _get_prior_terms,
    _get_transaction_matrix,
    _get_prior_translation,
)


def _get_prior_iic0(
    term: str,
    boundaries: dict[str, set[str]],
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
    ``boundaries`` is a dict with two keys: `roots` and `leaves`.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A float ranging from 0 to 1.
    """
    subsumers_t = _get_prior_ancestors({term}, prior, data, params) | {term}
    leaves_t = _get_prior_descendants({term}, prior, data, params).intersection(
        boundaries["leaves"]
    )
    tmp = ((len(leaves_t) / len(subsumers_t)) + 1) / (len(boundaries["leaves"]) + 1)
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
    boundaries = _get_prior_boundaries(prior, data, params)
    tmp = {t: _get_prior_iic0(t, boundaries, prior, data, params) for t in terms}
    out = pd.DataFrame.from_dict(tmp, orient="index", columns=["iic"])
    return out


def get_iics(
    data: dict[str, Any], params: dict[str, Any] = get_parameters()
) -> pd.DataFrame:
    """Get a table associating each annotation term to its raw and normalized IICs, its prior and
    its human-readable label.

    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame where rows are terms, with four columns: `iic`, `label`, `n_iic` and `prior`.
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
    return pd.concat(tmp)


def _setup_eics(
    associations: pd.DataFrame,
    data: dict[str, Any],
    params: dict[str, Any] = get_parameters(),
) -> dict[str, dict[str, set[str]]]:
    """Setup the results of a GORI enrichment analysis to compute the Extrinsic Information Content (EIC) of terms.

    ``associations`` is a pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
        `n_genes`, `genes`, `url_a` and `url_c`.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A dict associating priors (keys) to their group-specific annotation terms (values).
    """
    f = np.vectorize(lambda C, P: _get_prior_inverse_translation(C, P, data, params))
    associations["inverse_consequents"] = f(
        associations.consequents, associations.prior_c
    )
    corpus = {
        prior: group.groupby("antecedents")["inverse_consequents"].apply(set).to_dict()
        for prior, group in associations.groupby("prior_c")
    }
    return corpus


def _get_prior_eics(
    prior: str,
    corpus: dict[str, dict[str, set[str]]],
    data: dict[str, Any],
    params: dict[str, Any],
) -> pd.DataFrame:
    """Get the Extrinsic Information Content (EIC) of every term from a prior.

    ``prior`` is a prior label.
    ``corpus`` is a dict associating some priors (keys) to their group-specific annotation terms (values).
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame where rows are terms, with four columns: `eic`, `n_eic`, `label` and `prior`.
    """
    tmp = {prior: {}}  # type: dict[str, dict[str, set[str]]]
    for g, terms in corpus[prior].items():
        lineage = _get_prior_ancestors(terms, prior, data, params) | terms
        tmp[prior][g] = lineage

    tm = _get_transaction_matrix(tmp)
    tm = tm.fillna(0)
    tm = tm.mean(axis=0)
    tm = tm.apply(lambda x: -log2(x))
    out = pd.DataFrame(tm, columns=["eic"])

    out["prior"] = prior
    out["n_eic"] = out.eic / out.eic.max()
    out["label"] = [_get_prior_translation(t, prior, data, params) for t in out.index]
    return out


def get_eics(
    corpus: dict[str, dict[str, set[str]]],
    data: dict[str, Any],
    params: dict[str, Any] = get_parameters(),
) -> pd.DataFrame:
    """Get the Extrinsic Information Content (EIC) of every term in a corpus. It ranges from 0 to +Inf.

    EIC(t) = -log[ p(t) ], where:
        `EIC(t)` is the extrinsic information content of the term t,
        `p(t)` is the frequency of a term t in the corpus.

    Resnik, Philip.
    "Using information content to evaluate semantic similarity in a taxonomy."
    arXiv preprint cmp-lg/9511007 (1995).
    See: https://doi.org/10.48550/arXiv.cmp-lg/9511007

    ``prior`` is a prior label.
    ``corpus`` is a dict associating some priors (keys) to their group-specific annotation terms (values).
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame where rows are terms, with four columns: `eic`, `n_eic`, `label` and `prior`.
    """
    tmp = []  # type: list[pd.DataFrame]
    for p in corpus.keys():
        tmp.append(_get_prior_eics(p, corpus, data, params))
    eics = pd.concat(tmp)

    corpus_terms = set()  # type: set[str]
    for prior, _tmp in corpus.items():
        for group, terms in _tmp.items():
            corpus_terms = corpus_terms | terms

    eics = eics.loc[list(corpus_terms), :]
    return eics
