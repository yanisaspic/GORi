"""Functions called to download and setup knowledge bases prior to the GORi analysis.

    2025/05/23 @yanisaspic"""

import os
import shutil
from typing import Any
from gori.params import get_parameters
from gori.src.utils import _get_timestamp


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
    remove_dl: bool = True,
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Set-up priors exploitable for a GORi annotation analysis.

    ``priors`` is a set of valid knowledge categories (e.g. CTYP).
    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    ``remove_dl`` is a boolean indicating if downloaded files should be removed after the setup.
    ``params`` is a dict of parameters.

    """
    setup_wrapper = params["wrappers"]["setup_wrapper"]
    if not os.path.isdir(su_path):
        os.mkdir(su_path)
    for p in priors:
        setup_wrapper[p](dl_path, su_path)
    if remove_dl:
        shutil.rmtree(dl_path)