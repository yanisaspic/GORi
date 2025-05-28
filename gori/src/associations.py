"""Functions called to obtain every association between two knowledge bases.

    2025/05/13 @yanisaspic"""

import numpy as np
import pandas as pd
import numpy.typing as npt
from typing import Any
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection
from mlxtend.frequent_patterns import fpgrowth, association_rules
from gori.src.utils import (
    _get_prior_ancestors,
    _get_prior_descendants,
    _get_prior_translation,
    _get_prior_url,
    _get_transaction_matrix,
)


def _get_all_associations(
    transaction_matrix: pd.DataFrame, n_genes_threshold: int
) -> pd.DataFrame:
    """Get all associations from a transaction matrix.

    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.
    ``n_genes_threshold`` is the minimum number of genes to consider an association.

    Returns
        A pd.DataFrame with three columns: `antecedents`, `consequents` and `lift`.
    """
    minimum_support = n_genes_threshold / len(transaction_matrix.index)
    frequent_itemsets = fpgrowth(
        transaction_matrix, min_support=minimum_support, max_len=2, use_colnames=True
    )
    associations = association_rules(frequent_itemsets, metric="lift", min_threshold=1)
    associations = associations[["antecedents", "consequents", "lift"]]
    for col in "antecedents", "consequents":
        associations[col] = associations[col].apply(
            lambda x: next(iter(x))
        )  # remove the frozenset
    return associations


def _get_target_associations(
    associations: pd.DataFrame, antecedent_prior: str, consequent_prior: str
) -> pd.DataFrame:
    """Get associations between two specific knowledge bases.

    ``associations`` is a pd.DataFrame with three columns: `antecedents`, `consequents` and `lift`.
    ``antecedent_prior`` and ``consequent_prior`` are prior labels.

    Returns
        A pd.DataFrame with five columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c`.
    """
    get_prior = lambda x: x.split(":", 1)[0]
    for col in ["antecedents", "consequents"]:
        associations[f"prior_{col[0]}"] = associations[col].apply(get_prior)
    associations = associations.loc[
        (associations.prior_a == antecedent_prior)
        & (associations.prior_c == consequent_prior)
    ]
    return associations


def _get_contingency_matrix(
    A: str, C: str, transaction_matrix: pd.DataFrame
) -> npt.NDArray[np.int_]:
    """Get a contingency matrix of two annotations.

    ``a`` and ``b`` are two annotations.
    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.

    Returns
        A contingency matrix as a numpy array.
    """
    X = transaction_matrix[A].to_numpy()
    Y = transaction_matrix[C].to_numpy()
    a = np.sum(np.logical_and(X, Y))
    b = np.sum(X) - a
    c = np.sum(Y) - a
    d = len(transaction_matrix.index) - a - b - c
    contingency_matrix = np.array([[a, b], [c, d]])
    return contingency_matrix


def _get_fisher_exact(A: str, C: str, transaction_matrix: pd.DataFrame) -> float:
    """Get the p-value of a Fisher's exact test.

    ``A`` and ``C`` are two annotations.
    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.

    Returns
        A p-value.
    """
    contingency_matrix = _get_contingency_matrix(A, C, transaction_matrix)
    pval = fisher_exact(contingency_matrix, alternative="greater")[1]
    return pval


def _get_significant_associations(
    associations: pd.DataFrame,
    transaction_matrix: pd.DataFrame,
    params: dict[str, Any],
    fdr: bool = False,
) -> pd.DataFrame:
    """Get significant associations.

    ``associations`` is a pd.DataFrame with five or more columns: `antecedents`, `consequents`, `lift`, `prior_a` and `prior_c`.
    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.
    ``params`` is a dict of parameters.
    ``fdr`` is a boolean indicating if the p-values should be corrected for multiple testing.

    Returns
        A pd.DataFrame with six or more columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c`,
        `pval` and `fdr`.
    """
    f = np.vectorize(lambda A, C: _get_fisher_exact(A, C, transaction_matrix))
    associations["pval"] = f(associations.antecedents, associations.consequents)
    if fdr:
        associations["fdr"] = fdrcorrection(associations.pval)[1]
        return associations.loc[associations.fdr < params["pvalue_threshold"]]
    return associations.loc[associations.pval < params["pvalue_threshold"]]


def _get_lift_threshold(
    associations: pd.DataFrame, transaction_matrix: pd.DataFrame, params: dict[str, Any]
) -> float:
    """Get a lift threshold using a heuristic approach. Associations with a higher lift will be kept.

    ``associations`` is a pd.DataFrame with five columns: `antecedents`, `consequents`, `lift`, `prior_a` and `prior_c`.
    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.
    ``params`` is a dict of parameters.

    Returns
        A float.
    """
    associations = associations.reset_index(drop=True)
    f = np.vectorize(lambda A, C: _get_fisher_exact(A, C, transaction_matrix))
    strong_index = 0  # index of the weakest `strong` association
    weak_index = (
        len(associations.index) - 1
    )  # index of the strongest `weak` association

    while weak_index - strong_index > 1:  # binary search
        frontier_index = (weak_index + strong_index) // 2
        tmp = associations.loc[associations.lift == associations.lift[frontier_index]]
        pvals = f(tmp.antecedents, tmp.consequents)
        if np.median(pvals) < params["pvalue_threshold"]:
            strong_index = tmp.index[-1]
        else:
            weak_index = tmp.index[0]

    lift_threshold = associations.lift.iloc[weak_index]
    return lift_threshold


def _get_strong_associations(
    associations: pd.DataFrame, transaction_matrix: pd.DataFrame, params: dict[str, Any]
) -> pd.DataFrame:
    """Get strong associations.

    ``associations`` is a pd.DataFrame with five columns: `antecedents`, `consequents`, `lift`, `prior_a` and `prior_c`.
    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame with five or more columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c` and `pval`
    """
    if params["use_heuristic"]:
        lift_threshold = _get_lift_threshold(associations, transaction_matrix, params)
        return associations.loc[associations.lift > lift_threshold]
    associations = _get_significant_associations(
        associations, transaction_matrix, params, fdr=False
    )
    return associations


def _get_unique_terms(
    terms: set[str], prior: str, data: dict[str, Any], params: dict[str, Any]
) -> set[str]:
    """Get unique terms, i.e. terms that are not redundant with each other.

    ``terms`` is a set of annotation terms.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A set of annotation terms.
    """
    tmp = {t.split(":", 1)[1] for t in terms}
    ancestors = _get_prior_ancestors(tmp, prior, data, params)
    unique_terms = tmp.difference(ancestors)
    out = {f"{prior}:{t}" for t in unique_terms}
    return out


def _get_unique_indexes(
    a: str, associations: pd.DataFrame, data: dict[str, Any], params: dict[str, Any]
) -> pd.DataFrame:
    """Get unique indexes.

    ``a`` is an annotation term, and antecedent.
    ``associations`` is a pd.DataFrame with six columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c` and `pval`.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame with six columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c` and
        `pval`.
    """
    tmp = associations.loc[associations.antecedents == a]
    terms = set(tmp.consequents)
    unique_terms = _get_unique_terms(terms, associations.prior_c.iloc[0], data, params)

    # remove t_{n+x} and t_{n-x} from consequents:
    ancestors = _get_prior_ancestors({a}, associations.prior_a.iloc[0], data, params)
    descendants = _get_prior_descendants(
        {a}, associations.prior_a.iloc[0], data, params
    )
    throwaway = ancestors | descendants

    unique_terms = unique_terms.difference(throwaway)
    out = tmp.index[tmp.consequents.isin(unique_terms)].to_numpy()
    return out


def _get_unique_associations(
    associations: pd.DataFrame, data: dict[str, Any], params: dict[str, Any]
) -> pd.DataFrame:
    """Get unique associations.

    ``associations`` is a pd.DataFrame with six columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c` and `pval`.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame with six columns: `antecedents`, `consequents`, `lift`, `prior_a`, `prior_c` and
        `pval`.
    """
    f = np.vectorize(
        lambda A: _get_unique_indexes(A, associations, data, params),
        otypes=[np.ndarray],
    )
    unique_indexes = f(associations.antecedents.unique())
    unique_indexes = np.hstack(unique_indexes)
    associations = associations.loc[unique_indexes]
    return associations


def _get_association_genes(
    a: str, c: str, transaction_matrix: pd.DataFrame
) -> list[str]:
    """Get the genes associated with two annotations.

    ``a`` and ``c`` are two annotations.
    ``transaction_matrix`` is a binary pd.DataFrame where rows are genes and columns are annotations.

    Returns
        A list of genes.
    """
    tmp = transaction_matrix[a] & transaction_matrix[c]
    genes = transaction_matrix.index[tmp].to_list()
    return genes


def _enrich_associations(
    associations: pd.DataFrame,
    data: dict[str, Any],
    transaction_matrix: pd.DataFrame,
    params: dict[str, Any],
) -> pd.DataFrame:
    """Enrich associations with their genes, urls and human-readable translations.

    ``associations`` is a pd.DataFrame with seven columns: `antecedents`, `consequents`, `lift`,
    `prior_a`, `prior_c`, `pval` and `fdr`.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame with eleven columns: `antecedents`, `consequents`, `lift`,
        `prior_a`, `prior_c`, `pval`, `fdr`, `n_genes`, `genes`, `url_a` and `url_c`.
    """
    f = np.vectorize(
        lambda A, C: _get_association_genes(A, C, transaction_matrix),
        otypes=[np.ndarray],
    )
    genes = f(associations.antecedents, associations.consequents)
    associations["n_genes"] = [len(g) for g in genes]
    associations["genes"] = [", ".join(g) for g in genes]

    g = np.vectorize(lambda X, P: _get_prior_url(X, P, data, params))
    h = np.vectorize(
        lambda X, P: _get_prior_translation(X, P, data, params, has_prefix=True)
    )
    for col in ["antecedents", "consequents"]:
        p = f"prior_{col[0]}"
        associations[f"url_{col[0]}"] = g(associations[col], associations[p])
        associations[col] = h(associations[col], associations[p])

    return associations


def _get_associations(
    antecedent_prior: str,
    consequent_prior: str,
    annotations: dict[str, dict[str, set[str]]],
    data: dict[str, Any],
    params: dict[str, Any],
) -> dict[str, Any]:
    """Get relevant associations between two knowledge bases.

    ``antecedent_prior`` and ``consequent_prior`` are prior labels.
    ``annotations`` is a dict associating priors (keys) to their gene-specific annotations (values).
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.

    Returns
        A pd.DataFrame with nine columns: `antecedents`, `consequents`, `lift`, `pval`, `fdr`,
        `n_genes`, `genes`, `url_a` and `url_c`
    """
    associations_counter = {"prior": consequent_prior}  # type: dict[str, Any]
    transaction_matrix = _get_transaction_matrix(annotations)

    while True:
        associations = _get_all_associations(
            transaction_matrix, params["n_genes_threshold"]
        )
        if len(associations.index) == 0:
            print(
                f"No associations found with {params['n_genes_threshold']} genes or more."
            )
            break
        # Associations are filtered by direction
        associations = _get_target_associations(
            associations, antecedent_prior, consequent_prior
        )
        associations_counter["target"] = len(associations.index)
        if len(associations.index) == 0:
            print(
                f"No associations found between {antecedent_prior} and {consequent_prior}."
            )
            break

        # Associations are filtered by strength
        associations = associations.sort_values(
            by=["lift", "antecedents", "consequents"], ascending=[False, True, True]
        )  # sort by association strength
        associations = _get_strong_associations(
            associations, transaction_matrix, params
        )
        associations_counter["strong"] = len(associations.index)
        if len(associations.index) == 0:
            print(
                f"No strong associations found using {params['use_heuristic']} approach."
            )
            break

        # Associations are filtered by redundancy
        associations = _get_unique_associations(associations, data, params)
        associations_counter["unique"] = len(associations.index)
        # Associations are filtered by statistical significance
        associations = _get_significant_associations(
            associations, transaction_matrix, params, fdr=True
        )
        associations_counter["corrected"] = len(associations.index)
        if len(associations.index) == 0:
            print(f"No associations found with corrected p-values.")
            break

        associations = _enrich_associations(
            associations, data, transaction_matrix, params
        )
        break

    print(associations_counter)
    out = {"associations": associations, "associations_counter": associations_counter}
    return out
