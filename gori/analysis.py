"""Main functions called to run a GORi analysis.

    2025/05/20 @yanisaspic"""

import pandas as pd
from typing import Any, Optional
from gori.loaders import load_priors, load_feve
from gori.src.analysis._get_associations import _get_associations
from gori.src.analysis._get_annotations import (
    _get_annotations,
    _get_annotations_counter,
)


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


def gori(
    geneset: set[str],
    antecedent_prior: str,
    consequent_priors: set[str],
    data: dict[str, Any],
    params: dict[str, Any] = get_parameters(),
    sheets: bool = True,
) -> Optional[dict[str, pd.DataFrame]]:
    """Conduct a GORi enrichment analysis.

    ``geneset`` is a set of gene symbols or UniProtIDs.
    ``antecedent_prior`` is a prior label.
    ``consequent_priors`` is a set of prior labels.
    ``data`` is a dict associating priors (keys) to their contents (values).
    ``params`` is a dict of parameters.
    ``sheets`` is a boolean indicating if the results should be saved in a spreadsheet.

    Returns
        A dict with three keys: `annotations_counter`, `associations` and `associations_counter`.
    """
    results = {}
    outs = {
        "associations": [],
        "associations_counter": [],
    }  # type: dict[str, list[Any]]
    annotations = _get_annotations(geneset, data, params["use_gene_symbol"])
    results["annotations_counter"] = _get_annotations_counter(annotations)

    for c in consequent_priors:
        print(f"{antecedent_prior} -> {c}")

        tmp = results["annotations_counter"][[antecedent_prior, c]]
        if tmp.sum().sum() == 0:
            print(f"No annotations found using {antecedent_prior} and {c}.")
            continue

        c_annotations = {
            k: v for k, v in annotations.items() if k in [antecedent_prior, c]
        }
        c_out = _get_associations(antecedent_prior, c, c_annotations, data, params)
        outs["associations"].append(c_out["associations"])
        outs["associations_counter"].append(c_out["associations_counter"])

    if len(outs["associations"]) == 0:
        return None

    results["associations"] = pd.concat(outs["associations"])
    associations_counter = pd.DataFrame(outs["associations_counter"])
    associations_counter = associations_counter.fillna(0).sort_values(by="prior")
    results["associations_counter"] = associations_counter

    if sheets:
        with pd.ExcelWriter(params["sheets_path"]) as writer:
            for sheet in results.keys():
                if sheet == "annotations_counter":
                    results[sheet].to_excel(writer, sheet_name=sheet)
                else:
                    results[sheet].to_excel(writer, sheet_name=sheet, index=False)
    return results


def gorilon(
    path: str,
    direction: str = "any",
    priors: set[str] = {"BIOP", "CTYP", "GENG", "PATH"},
    params: dict[str, Any] = get_parameters(),
    sheets: bool = True,
) -> Optional[dict[str, pd.DataFrame]]:
    """Conduct a GORi analysis on the results of a fEVE analysis.

    ``path`` is the path to the .xlsx file containing the results of a fEVE analysis.
    ``direction`` is the direction of the genes' expression: one of 'up', 'down' or 'any'.
    ``priors`` is a set of prior labels.
    ``params`` is a dict of parameters.
    ``sheets`` is a boolean indicating if the results should be saved in a spreadsheet.

    Returns
        A dict with three keys: `annotations_counter`, `associations` and `associations_counter`.
    """
    data = load_priors(priors)
    feve = load_feve(path, direction)
    data["FEVE"] = feve
    geneset = data["FEVE"]["annotations"].keys()
    results = gori(geneset, "FEVE", priors, data, params, sheets)
    return results
