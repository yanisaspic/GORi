"""Functions used by download_wrapper().

    2025/05/20 @yanisaspic"""

import os
import gzip
import pandas as pd
from typing import Any, Callable
from shutil import copyfileobj
from urllib.request import urlretrieve
from gori.src.utils import _get_timestamp


def _download_cell_types(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to cell types (CTYP).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/cell_types"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["resources_wrapper"]["CTYP"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            urlretrieve(url, f"{_path}/{file}")
            log.write(f"\t\t {_get_timestamp()}: Downloaded {file} from {url}\n")

    annotations = pd.read_csv(f"{_path}/raw_CellTaxonomy_annotations.txt", sep="\t")
    annotations = annotations.loc[annotations.Species == "Homo sapiens"]
    annotations = annotations[["Specific_Cell_Ontology_ID", "Uniprot"]]
    annotations = annotations.dropna()
    annotations.to_csv(f"{_path}/CellTaxonomy_annotations.csv")
    os.remove(f"{_path}/raw_CellTaxonomy_annotations.txt")


def _download_diseases(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to diseases (DISE).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/diseases"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["resources_wrapper"]["DISE"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            urlretrieve(url, f"{_path}/{file}")
            log.write(f"\t\t {_get_timestamp()}: Downloaded {file} from {url}\n")

    with gzip.open(f"{_path}/raw_CTD_annotations.csv.gz", "rb") as f_in:
        with open(f"{_path}/raw_CTD_annotations.csv", "wb") as f_out:
            copyfileobj(f_in, f_out)
    os.remove(f"{_path}/raw_CTD_annotations.csv.gz")

    with open(f"{_path}/raw_CTD_annotations.csv", "r") as f_in:
        with open(f"{_path}/CTD_annotations.csv", "w") as f_out:
            for line in f_in.readlines():
                if "MESH:" in line:
                    f_out.write(line)
    os.remove(f"{_path}/raw_CTD_annotations.csv")


def _download_gene_groups(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to gene groups (GENG).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/gene_groups"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["resources_wrapper"]["GENG"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            urlretrieve(url, f"{_path}/{file}")
            log.write(f"\t\t {_get_timestamp()}: Downloaded {file} from {url}\n")


def _download_pathways(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to pathways (PATH).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/pathways"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["resources_wrapper"]["PATH"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            urlretrieve(url, f"{_path}/{file}")
            log.write(f"\t\t {_get_timestamp()}: Downloaded {file} from {url}\n")
