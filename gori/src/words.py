"""Functions called to identify relevant annotation words.

    2025/05/23 @yanisaspic"""

import re
import numpy as np
import pandas as pd
from typing import Any
from gori.src.utils import (
    _get_transaction_matrix,
    _get_resource_boundaries,
    _get_resource_translation,
)


def _get_words0(label: str, params: dict[str, Any]) -> set[str]:
    """Get a subset of words from a human-readable label.

    ``label`` is a human-readable label.
    ``params`` is a dict of parameters.

    Returns
        A set of words.
    """
    a = re.sub(r"[^a-zA-Z0-9]", " ", label)
    b = re.sub(r"[^a-zA-Z0-9-]", " ", label)
    tokens = set(a.split()) | set(b.split())
    tmp = {t.lower() for t in tokens if not t.isnumeric()}
    out = tmp.difference(params["stopwords"])
    return out


def _get_words_scores(
    associations: pd.DataFrame, data: dict[str, Any], params: dict[str, Any]
) -> pd.DataFrame:
    """Get a table reporting the group-specific scores of each word.

    ``associations`` is a pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
        `n_genes`, `genes`, `url_a` and `url_c`.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame where rows are groups, columns are words, and cell values are scores.
    """
    pairs = associations.groupby("antecedents")["consequents"].apply(set).to_dict()
    weights = {
        a: dict(group[["consequents", "lift"]].values)
        for a, group in associations.groupby("antecedents")
    }

    _scores = []
    for group, terms in pairs.items():
        words = {t: _get_words0(t, params) for t in terms}
        wm = _get_transaction_matrix({"tmp": words})
        wm = wm.astype(bool).replace({True: 1, False: np.nan})
        wm = wm.add(wm.index.map(weights[group]), axis=0)  # add to increment the scores
        wm = wm.mean(axis=0)
        wm.name = group
        _scores.append(wm)

    scores = pd.concat(_scores, axis=1)
    scores.index.name = "word"
    return scores
