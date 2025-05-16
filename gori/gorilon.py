"""Function called to conduct a GORi analysis on the results of a fEVE analysis.

    2025/05/16 @yanisaspic"""

import pandas as pd
from gori.gori import gori
from typing import Any, Optional
from gori.priors.load_feve import load_feve
from gori.priors.load_priors import load_priors
from gori.get_parameters import get_parameters

def gorilon(
    path: str,
    direction: str = "any",
    priors: set[str] = {"BIOP", "CTYP", "GENG", "PATH"},
    params: dict[str, Any] = get_parameters(),
    sheets: bool = True
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