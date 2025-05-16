"""Function called to run a GORi analysis.

    2025/05/16 @yanisaspic"""

import pandas as pd
from typing import Any, Optional
from gori.mining.get_annotations import get_annotations, get_annotations_counter
from gori.priors.load_feve import load_feve
from gori.get_parameters import get_parameters
from gori.priors.load_priors import load_priors
from gori.mining.get_associations import get_associations


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
    outs = {"associations": [], "associations_counter": []}
    annotations = get_annotations(geneset, data, params["use_gene_symbol"])
    results["annotations_counter"] = get_annotations_counter(annotations)

    for c in consequent_priors:
        print(f"{antecedent_prior} -> {c}")

        tmp = results["annotations_counter"][[antecedent_prior, c]]
        if tmp.sum().sum() == 0:
            print(f"No annotations found using {antecedent_prior} and {c}.")
            continue

        c_annotations = {
            k: v for k, v in annotations.items() if k in [antecedent_prior, c]
        }
        c_out = get_associations(antecedent_prior, c, c_annotations, data, params)
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