"""Function called to get the default parameters for a GORi analysis.

    2025/05/16 @yanisaspic"""

from typing import Any

def get_parameters() -> dict[str, Any]:
    """Get parameters for a GORi enrichment analysis.

    The parameters are stored in a dict:
    `n_genes_threshold` and `pvalues_threshold` are numeric. They are used to filter out weak associations.
    `heuristic` is a boolean indicating if the heuristic approach should be used by GORi.
    `use_gene_symbol` is a boolean indicating if the gene symbols should be used.
    `sheets_path` is a path where the results of the analysis will be stored.

    Returns
        A dict with five keys: `n_genes_threshold`, `pvalue_threshold`, `use_heuristic`,
        `use_gene_symbol` and `sheets_path`.
    """
    parameters = {
        "n_genes_threshold": 5,
        "pvalue_threshold": 0.05,
        "use_heuristic": True,
        "use_gene_symbol": True,
        "sheets_path": "./GORi.xlsx",
    }
    return parameters