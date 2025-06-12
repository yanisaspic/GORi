"""Functions called to download and setup the knowledge base resources for the GORi analysis.

    2025/06/10 @yanisaspic"""

import os
import shutil
from typing import Any
from gori.params import get_parameters
from gori.src.utils import _get_timestamp


def download_resources(
    resources: set[str] = {"CellMarker2", "CellTaxonomy", "HGNC", "Reactome"},
    path: str = "./.resources",
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Download knowledge bases exploitable for a GORi annotation analysis.

    ``resources`` is a set of valid knowledge categories (e.g. CellMarker2).
    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    download_wrapper = params["wrappers"]["download_wrapper"]
    log = open("./downloads.log", "w")
    log.write(f"Downloading resources ({_get_timestamp()})\n")
    log.close()

    if not os.path.isdir(path):
        os.mkdir(path)

    for p in resources:
        log = open("./downloads.log", "a")
        log.write(f"\t +++ {p} resource ({_get_timestamp()})\n")
        log.close()

        download_wrapper[p](path, params)

        log = open("./downloads.log", "a")
        log.write(f"\t --- DONE ({_get_timestamp()})\n")
        log.close()

    log = open("./downloads.log", "a")
    log.write(f"DONE ({_get_timestamp()})")
    log.close()


def setup_resources(
    resources: set[str] = {"CellMarker2", "CellTaxonomy", "HGNC", "Reactome"},
    dl_path: str = "./.resources",
    su_path: str = "./resources",
    remove_dl: bool = True,
    params: dict[str, Any] = get_parameters(),
) -> None:
    """Set-up resources exploitable for a GORi annotation analysis.

    ``resources`` is a set of valid knowledge categories (e.g. CellMarker2).
    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    ``remove_dl`` is a boolean indicating if downloaded files should be removed after the setup.
    ``params`` is a dict of parameters.

    """
    setup_wrapper = params["wrappers"]["setup_wrapper"]
    if not os.path.isdir(su_path):
        os.mkdir(su_path)
    for p in resources:
        setup_wrapper[p](dl_path, su_path)
    if remove_dl:
        shutil.rmtree(dl_path)
