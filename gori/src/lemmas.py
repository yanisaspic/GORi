"""Functions called to identify relevant annotation lemmas.

    2025/05/23 @yanisaspic"""

import re
import nltk
import numpy as np
import pandas as pd
from typing import Any
from nltk.corpus import stopwords
from nltk.stem import WordNetLemmatizer
from gori.src.utils import (
    _get_transaction_matrix,
    _get_prior_boundaries,
    _get_prior_translation,
)


def _get_lemmas0(label: str) -> set[str]:
    """Get a subset of lemmas from a human-readable label.

    ``label`` is a human-readable label.

    Returns
        A set of lemmas.
    """
    lemmatizer = WordNetLemmatizer()
    a = re.sub(r"[^a-zA-Z0-9]", " ", label)
    b = re.sub(r"[^a-zA-Z0-9-]", " ", label)
    tokens = set(a.split()) | set(b.split())
    tmp = {t.lower() for t in tokens if not t.isnumeric()}
    lemmas = {lemmatizer.lemmatize(t) for t in tmp}
    stop_words = set(stopwords.words("english"))
    out = {l for l in lemmas if l not in stop_words}
    return out


def _get_lemmas_scores(associations: pd.DataFrame) -> pd.DataFrame:
    """Get a table reporting the group-specific scores of each lemma.

    ``associations`` is a pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
        `n_genes`, `genes`, `url_a` and `url_c`.

    Returns
        A pd.DataFrame where rows are groups, columns are lemmas, and cell values are scores.
    """
    pairs = associations.groupby("antecedents")["consequents"].apply(set).to_dict()
    weights = {
        a: dict(group[["consequents", "lift"]].values)
        for a, group in associations.groupby("antecedents")
    }

    _scores = []
    for group, terms in pairs.items():
        lemmas = {t: _get_lemmas0(t) for t in terms}
        lm = _get_transaction_matrix({"tmp": lemmas})
        lm = lm.replace({True: 1, False: np.nan})
        lm = lm.add(lm.index.map(weights[group]), axis=0)  # add to increment the scores
        lm = lm.mean(axis=0)
        lm.name = group
        _scores.append(lm)

    scores = pd.concat(_scores, axis=1)
    return scores


def _get_top_lemmas(
    lemmas_scores: pd.DataFrame, data: dict[str, Any], params: dict[str, Any]
) -> pd.DataFrame:
    """Get a table ranking the best lemmas for each group.

    ``associations`` is a pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
        `n_genes`, `genes`, `url_a` and `url_c`.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame where rows are ranks, columns are groups, and cell values are lemmas.
    """
    root_lemmas = set()  # type: set[str]
    for p in data.keys():
        tmp = _get_prior_boundaries(p, data, params)["roots"]
        roots = {
            _get_prior_translation(r, p, data, params, has_prefix=False) for r in tmp
        }
        for r in roots:
            root_lemmas = root_lemmas | _get_lemmas0(r)

    scores = lemmas_scores.drop(root_lemmas, axis=0, errors="ignore")
    scores.index.name = "lemma"

    _tmp = {}   # type: dict[str, pd.Series]
    for group in scores.columns:
        scores = scores.sort_values(by=[group, "lemma"], ascending=[False, True])
        _tmp[group] = scores[group].loc[scores[group] > 0].index.tolist()

    top_lemmas = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in _tmp.items()]))
    return top_lemmas
