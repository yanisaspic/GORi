"""Functions called to download and setup knowledge bases prior to the GORi analysis.

    2025/05/20 @yanisaspic"""

import os
from typing import Callable
from gori.src.setup.download_priors import (
    _download_cell_types,
    _download_diseases,
    _download_gene_groups,
    _download_pathways,
    _get_timestamp,
)
from gori.src.setup.setup_priors import (
    _setup_cell_types,
    _setup_diseases,
    _setup_gene_groups,
    _setup_pathways,
)


def get_download_wrapper() -> dict[str, Callable]:
    """A wrapper to download curated knowledge bases (i.e. priors).

    Returns
        A dict associating prior labels (keys) to their download functions (values).
    """
    return {
        "CTYP": _download_cell_types,
        "DISE": _download_diseases,
        "GENG": _download_gene_groups,
        "PATH": _download_pathways,
    }


def download_priors(
    priors: set[str] = {"CTYP", "GENG", "PATH"},
    path: str = "./.priors",
    download_wrapper=get_download_wrapper(),
) -> None:
    """Download knowledge bases exploitable for a GORi annotation analysis.

    ``priors`` is a set of valid knowledge categories (e.g. CTYP).
    ``path`` is a path to store the downloaded files.
    ``download_wrapper`` is a dict associating prior labels (keys) to their download functions (values).
    """
    log = open("./downloads.log", "w")
    log.write(f"{_get_timestamp()}: Running downloads_priors.py\n")
    log.close()

    if not os.path.isdir(path):
        os.mkdir(path)

    for p in priors:
        log = open("./downloads.log", "a")
        log.write(f"\t {_get_timestamp()}: Downloading {p}\n")
        log.close()
        download_wrapper[p](path)

    log = open("./downloads.log", "a")
    log.write(f"{_get_timestamp()}: Done")
    log.close()


def get_setup_wrapper() -> dict[str, Callable]:
    """A wrapper to set-up curated knowledge bases (i.e. priors).

    Returns
        A dict associating prior labels (keys) to their set-up functions (values).
    """
    return {
        "CTYP": _setup_cell_types,
        "DISE": _setup_diseases,
        "GENG": _setup_gene_groups,
        "PATH": _setup_pathways,
    }


def setup_priors(
    priors: set[str] = {"CTYP", "GENG", "PATH"},
    dl_path: str = "./.priors",
    su_path: str = "./priors",
    setup_wrapper=get_setup_wrapper(),
) -> None:
    """Set-up priors exploitable for a GORi annotation analysis.

    ``priors`` is a set of valid knowledge categories (e.g. CTYP).
    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    ``setup_wrapper`` is a dict associating prior labels (keys) to their set-up functions (values).
    """
    if not os.path.isdir(su_path):
        os.mkdir(su_path)
    for p in priors:
        setup_wrapper[p](dl_path, su_path)
