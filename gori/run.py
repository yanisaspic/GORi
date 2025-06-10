"""Main functions called to run a GORi analysis.

    2025/05/20 @yanisaspic"""

import pandas as pd
from typing import Any, Optional
from gori.params import get_parameters
from gori.loaders import load_resources, load_feve
from gori.src.words import _get_words_scores
from gori.src.notebook import _write_notebook
from gori.src.associations import _get_associations
from gori.src.annotations import _get_annotations, _get_annotations_counter


def gori(
    geneset: set[str],
    antecedent_resource: str,
    consequent_resources: set[str],
    data: dict[str, Any],
    params: dict[str, Any] = get_parameters(),
    save: bool = True,
) -> Optional[dict[str, pd.DataFrame]]:
    """Conduct a GORi enrichment analysis.

    ``geneset`` is a set of gene symbols or UniProtIDs.
    ``antecedent_resource`` is a resource label.
    ``consequent_resources`` is a set of resource labels.
    ``data`` is a dict associating resources (keys) to their contents (values).
    ``params`` is a dict of parameters.
    ``save`` is a boolean indicating if the results should be saved in a spreadsheet and a notebook.

    Returns
        A dict with four keys: `annotations_counter`, `associations`, `associations_counter` and `words`.
    """
    results = {}
    outs = {
        "associations": [],
        "associations_counter": [],
    }  # type: dict[str, list[Any]]
    annotations = _get_annotations(geneset, data, params)
    results["annotations_counter"] = _get_annotations_counter(annotations)

    for c in consequent_resources:
        print(f"{antecedent_resource} -> {c}")

        tmp = results["annotations_counter"][[antecedent_resource, c]]
        if tmp.sum().sum() == 0:
            print(f"No annotations found using {antecedent_resource} and {c}.")
            continue

        c_annotations = {
            k: v for k, v in annotations.items() if k in [antecedent_resource, c]
        }
        c_out = _get_associations(antecedent_resource, c, c_annotations, data, params)
        outs["associations"].append(c_out["associations"])
        outs["associations_counter"].append(c_out["associations_counter"])

    if len(outs["associations"]) == 0:
        return None

    results["associations"] = pd.concat(outs["associations"])
    associations_counter = pd.DataFrame(outs["associations_counter"])
    associations_counter = associations_counter.fillna(0).sort_values(by="resource")
    results["associations_counter"] = associations_counter
    results["words"] = _get_words_scores(results["associations"], data, params)

    if save:
        with pd.ExcelWriter(params["sheets_path"]) as writer:
            for sheet in results.keys():
                if sheet in ["annotations_counter", "words"]:
                    results[sheet].to_excel(writer, sheet_name=sheet)
                else:
                    results[sheet].to_excel(writer, sheet_name=sheet, index=False)
        _write_notebook(params)
    return results


def gorilon(
    path: str,
    direction: str = "any",
    resources: set[str] = {
        "GO_BP",
        "GO_CC",
        "CellMarker2",
        "CellTaxonomy",
        "GO_MF",
        "Reactome",
    },
    resources_path: str = "./resources",
    params: dict[str, Any] = get_parameters(),
    save: bool = True,
) -> Optional[dict[str, pd.DataFrame]]:
    """Conduct a GORi analysis on the results of a fEVE analysis.

    ``path`` is the path to the .xlsx file containing the results of a fEVE analysis.
    ``direction`` is the direction of the genes' expression: one of 'up', 'down' or 'any'.
    ``resources`` is a set of resource labels.
    ``resources_path`` is the path to the JSON files containing the knowledge bases.
    ``params`` is a dict of parameters.
    ``save`` is a boolean indicating if the results should be saved in a spreadsheet and a notebook.

    Returns
        A dict with three keys: `annotations_counter`, `associations`, `associations_counter` and `words`.
    """
    data = load_resources(resources, resources_path, params)
    data["fEVE"] = load_feve(path, direction)
    geneset = data["fEVE"]["annotations"].keys()
    results = gori(geneset, "fEVE", resources, data, params, save)
    return results
