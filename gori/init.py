"""Functions called to download and setup knowledge bases prior to the GORi analysis.

    2025/05/23 @yanisaspic"""

import os
import nltk
import pandas as pd
from typing import Any, Callable
from gori.params import get_parameters
from gori.src.utils import _get_timestamp
from gori.metrics import _get_iics
from gori.etc.lemmas import _get_lemmas_iics


def download_priors(
    priors: set[str] = {"CTYP", "GENG", "PATH"},
    path: str = "./.priors",
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Download knowledge bases exploitable for a GORi annotation analysis.

    ``priors`` is a set of valid knowledge categories (e.g. CTYP).
    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    download_wrapper = params["wrappers"]["download_wrapper"]
    log = open("./downloads.log", "w")
    log.write(f"{_get_timestamp()}: Running downloads_priors.py\n")
    log.close()

    if not os.path.isdir(path):
        os.mkdir(path)

    for p in priors:
        log = open("./downloads.log", "a")
        log.write(f"\t {_get_timestamp()}: Downloading {p}\n")
        log.close()
        download_wrapper[p](path, params)

    log = open("./downloads.log", "a")
    log.write(f"{_get_timestamp()}: Done")
    log.close()


def setup_priors(
    priors: set[str] = {"CTYP", "GENG", "PATH"},
    dl_path: str = "./.priors",
    su_path: str = "./priors",
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Set-up priors exploitable for a GORi annotation analysis.

    ``priors`` is a set of valid knowledge categories (e.g. CTYP).
    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    ``params`` is a dict of parameters.

    """
    setup_wrapper = params["wrappers"]["setup_wrapper"]
    if not os.path.isdir(su_path):
        os.mkdir(su_path)
    for p in priors:
        setup_wrapper[p](dl_path, su_path)


def _fill_iics(
    iics: pd.DataFrame, data: dict[str, Any], params: dict[str, Any]
) -> pd.DataFrame:
    """Add new prior data to the Intrinstic Information Content (iic) table.

        ``iics`` is a pd.DataFrame with five columns: `term`, `label`, `iic`, `n_iic`, and `prior`.
        ``data`` is a dict associating priors (keys) to their contents (values).
        ``params`` is a dict of parameters.

    Returns
        A a pd.DataFrame with five columns: `term`, `label`, `iic`, `n_iic`, and `prior`.
    """
    missing = set(data.keys()).difference(set(iics.prior))
    tmp = {p: data[p] for p in missing}
    supplementary_table = _get_iics(data, params)
    table = pd.concat([iics, supplementary_table])
    return table


def setup_iics(
    data: dict[str, Any],
    iics_path: str = "./misc",
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Set-up the Intrinsic Information Content (IIC) exploitable for a GORi annotation analysis.

    ``data`` is a dict associating priors (keys) to their contents (values).
    ``iics_path`` is a path where the Intrinsic Information Content (IIC) table is stored.
    ``params`` is a dict of parameters.

    """
    if not os.path.isdir(iics_path):
        os.mkdir(iics_path)
    path = f"{iics_path}/IIC.csv"

    if not os.path.exists(path):
        iics = _get_iics(data, params)
        iics.to_csv(path, index=True)

    else:
        iics = pd.read_csv(path, index_col=0)
        tmp = {p: v for p, v in data.items() if p not in iics.prior.unique()}
        if len(tmp) > 0:
            suppl_iics = _get_iics(tmp, params)
            iics = pd.concat([iics, suppl_iics])
            iics.to_csv(path, index=True)


def setup_lemmas(
    data: dict[str, Any],
    iics_path: str = "./misc",
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Set-up the lemma matrix exploitable for a GORi annotation analysis.

    ``data`` is a dict associating priors (keys) to their contents (values).
    ``iics_path`` is a path where the Intrinsic Information Content (IIC) table is stored.
    ``params`` is a dict of parameters.

    """
    nltk.download("stopwords")
    nltk.download("punkt_tab")
    nltk.download("wordnet")
    iics = pd.read_csv(f"{iics_path}/IIC.csv", index_col=0)
    lemmas_iics = _get_lemmas_iics(iics)
    lemmas_iics.to_csv(f"{iics_path}/LIIC.csv", index=True)
