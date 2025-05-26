"""Main functions called to run a GORi analysis.

    2025/05/20 @yanisaspic"""

import pandas as pd
from typing import Any, Optional
from gori.params import get_parameters
from gori.loaders import load_priors, load_feve
from gori.src.associations import _get_associations
from gori.src.lemmas import _get_lemmas_scores, _get_top_lemmas
from gori.src.annotations import _get_annotations, _get_annotations_counter


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
    annotations = _get_annotations(geneset, data, params)
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
    results["lemmas_scores"] = _get_lemmas_scores(results["associations"])
    results["top_lemmas"] = _get_top_lemmas(results["lemmas_scores"], data, params)

    if sheets:
        with pd.ExcelWriter(params["sheets_path"]) as writer:
            for sheet in results.keys():
                if sheet in ["annotations_counter", "lemmas_scores"]:
                    results[sheet].to_excel(writer, sheet_name=sheet)
                else:
                    results[sheet].to_excel(writer, sheet_name=sheet, index=False)
    return results


def gorilon(
    path: str,
    direction: str = "any",
    priors: set[str] = {"BIOP", "CELC", "CTYP", "GENG", "MOLF", "PATH"},
    priors_path: str = "./priors",
    params: dict[str, Any] = get_parameters(),
    sheets: bool = True,
) -> Optional[dict[str, pd.DataFrame]]:
    """Conduct a GORi analysis on the results of a fEVE analysis.

    ``path`` is the path to the .xlsx file containing the results of a fEVE analysis.
    ``direction`` is the direction of the genes' expression: one of 'up', 'down' or 'any'.
    ``priors`` is a set of prior labels.
    ``priors_path`` is the path to the JSON files containing the knowledge bases.
    ``params`` is a dict of parameters.
    ``sheets`` is a boolean indicating if the results should be saved in a spreadsheet.

    Returns
        A dict with three keys: `annotations_counter`, `associations` and `associations_counter`.
    """
    data = load_priors(priors, priors_path, params)
    data["FEVE"] = load_feve(path, direction)
    geneset = data["FEVE"]["annotations"].keys()
    results = gori(geneset, "FEVE", priors, data, params, sheets)
    return results
