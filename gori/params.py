"""Function called to initialize a set of parameters for the GORi setup and analysis.

    2025/05/22 @yanisaspic"""

from typing import Any
from gori.wrappers import (
    ancestors_wrapper,
    annotations_wrapper,
    descendants_wrapper,
    download_wrapper,
    headers_wrapper,
    inverse_translate_wrapper,
    load_wrapper,
    resources_wrapper,
    roots_wrapper,
    setup_wrapper,
    terms_wrapper,
    translate_wrapper,
    urls_wrapper,
)


def _get_wrappers() -> dict[str, dict[str, Any]]:
    """Get the wrapper functions used by GORi to handle specific priors.

    Returns
        A dict associating wrapper labels (keys) to their prior-specific functions (values).
    """
    return {
        "ancestors_wrapper": ancestors_wrapper(),
        "annotations_wrapper": annotations_wrapper(),
        "descendants_wrapper": descendants_wrapper(),
        "download_wrapper": download_wrapper(),
        "headers_wrapper": headers_wrapper(),
        "inverse_translate_wrapper": inverse_translate_wrapper(),
        "load_wrapper": load_wrapper(),
        "resources_wrapper": resources_wrapper(),
        "roots_wrapper": roots_wrapper(),
        "setup_wrapper": setup_wrapper(),
        "terms_wrapper": terms_wrapper(),
        "translate_wrapper": translate_wrapper(),
        "urls_wrapper": urls_wrapper(),
    }


def get_parameters() -> dict[str, Any]:
    """Get parameters for a GORi enrichment analysis.

    The parameters are stored in a dict:
        `n_genes_threshold` and `pvalues_threshold` are numerics used to filter out weak associations.
        `heuristic` is a boolean indicating if the heuristic approach should be used by GORi.
        `use_gene_symbol` is a boolean indicating if the gene symbols should be used.
        `sheets_path` is a path where the results of the analysis will be stored.
        `wrappers` is a dict of prior-speicifc functions.

    Returns
        A dict with 6 keys: `n_genes_threshold`, `pvalue_threshold`, `use_heuristic`,
        `use_gene_symbol`, `sheets_path` and `wrappers`.
    """
    parameters = {
        "n_genes_threshold": 5,
        "pvalue_threshold": 0.05,
        "use_heuristic": True,
        "use_gene_symbol": True,
        "sheets_path": "./GORi.xlsx",
        "wrappers": _get_wrappers(),
    }
    return parameters
