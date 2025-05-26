"""Functions called to identify relevant annotation lemmas.

    2025/05/23 @yanisaspic"""

import re
import nltk
import numpy as np
import pandas as pd
from typing import Any, Optional
from nltk.stem import WordNetLemmatizer
from gori.src.utils import _get_transaction_matrix


def _get_lemmas0(label: str) -> set[str]:
    """Get a subset of lemmas from a human-readable label.

    ``label`` is a human-readable label.

    Returns
        A set of lemmas.
    """
    lemmatizer = WordNetLemmatizer()
    stop_words = set(nltk.corpus.stopwords.words("english"))
    label = re.sub("\W+", " ", label)
    tokens = nltk.word_tokenize(label)
    tmp = {t.lower() for t in tokens if t.isalpha()}
    _tmp = {lemmatizer.lemmatize(t) for t in tmp}
    lemmas = {t for t in _tmp if t not in stop_words}
    return lemmas


def _get_lemmas_iics(iics: pd.DataFrame) -> dict[str, set[str]]:
    """Get the Intrinsic Information Content (iic) of every lemma from an IIC table.

    ``iics`` is a pd.DataFrame with five columns: `term`, `label`, `iic`, `n_iic`, and `prior`.

    Returns
        A pd.Series associating lemmas (indexes) to their information content (values).
    """
    lemmas = {i: _get_lemmas0(iics.label[i]) for i in iics.index}
    lm = _get_transaction_matrix({"tmp": lemmas})
    lm = lm.replace({True: 1, False: np.nan})
    lm = lm.mul(lm.index.map(iics.n_iic), axis=0)
    lemmas_iics = lm.mean(axis=0) * lm.min(axis=0)  # multiply by min to penalize root lemmas
    lemmas_iics.name = "iic"
    return lemmas_iics


def _get_lemmas_lifts(associations: pd.DataFrame) -> dict[str, set[str]]:
    """Get the group-specific lifts of every lemma from an associations table.

    ``associations`` is a pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
        `n_genes`, `genes`, `url_a` and `url_c`

    Returns
        A pd.DataFrame where rows are groups, columns are lemmas, and cell values are lifts.
    """
    pairs = associations.groupby("antecedents")["consequents"].apply(set).to_dict()
    weights = {
        a: dict(group[["consequents", "lift"]].values)
        for a, group in associations.groupby("antecedents")
    }

    lifts = []
    for group, terms in pairs.items():
        lemmas = {t: _get_lemmas0(t) for t in terms}
        lm = _get_transaction_matrix({"tmp": lemmas})
        lm = lm.replace({True: 1, False: np.nan})
        lm = lm.mul(lm.index.map(weights[group]), axis=0)
        tmp = lm.mean(axis=0)
        tmp.name = group
        lifts.append(tmp)

    lemmas_lifts = pd.concat(lifts, axis=1)
    return lemmas_lifts


import numpy as np


def _get_lemmas_importances(associations: pd.DataFrame, iics_path: str = "./misc") -> pd.DataFrame:
    """Get the group-specific importances of every lemma.
    The most important lemmas are intrinsically informative and strongly associated to groups.

    ``associations`` is a pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
    `n_genes`, `genes`, `url_a` and `url_c`
    ``iics_path`` is a path where the Intrinsic Information Content (IIC) table is stored.

    Returns
        A pd.DataFrame where rows are groups, columns are lemmas, and cell values are importances.
    """
    lm_ic = pd.read_csv(f"{iics_path}/LIIC.csv", index_col=0)
    lm_ic = lm_ic.to_dict()
    lm_lift = _get_lemmas_lifts(associations)
    lm_lift = np.square(lm_lift)
    lm_imp = lm_lift.mul(lm_lift.index.map(lm_ic["iic"]), axis=0)
    return lm_imp
    return lm_lift


def _get_top_lemmas(lemmas_importances: pd.DataFrame, n: int=10) -> pd.DataFrame:
    """Get the top n lemmas of every group.

    ``lemmas_importances`` is a pd.DataFrame where rows are groups, columns are lemmas, and cell values are importances.

    Returns
        A pd.DataFrame where rows are ranks, columns are groups, and cell values are lemmas.
    """
    top_lemmas = pd.DataFrame({
        cluster: lemmas_importances[cluster].sort_values(ascending=False).head(n).index
        for cluster in lemmas_importances.columns
    })
    return top_lemmas