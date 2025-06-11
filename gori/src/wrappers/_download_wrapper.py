"""Functions used by download_wrapper().

    2025/05/20 @yanisaspic"""

import os
import gzip
import requests
import pandas as pd
from typing import Any
from shutil import copyfileobj
from gori.src.utils import _get_timestamp


def _download_file(url: str, path: str) -> None:
    """Download a file. This function is adapted for large files (>100MB), unlike urlretrieve.

    ``url`` is a url leading to a file to download.
    ``path`` is a path to store the downloaded file.
    """
    with requests.get(url, stream=True) as r:
        r.raise_for_status()  # ensures we notice bad responses
        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:  # filters out keep-alive new chunks
                    f.write(chunk)


def _download_cellmarker2_cell_types(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to cell types and CellMarker 2.0.

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.
    """
    _path = f"{path}/cell_types"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["wrappers"]["resources_wrapper"]["CellMarker2"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            log.write(f"\t\t Downloading {file} from {url} ({_get_timestamp()})\n")
            response = requests.get(url)
            with open(f"{_path}/{file}", "wb") as output:
                output.write(response.content)  
            log.write(f"\t\t\t DONE ({_get_timestamp()})")

    annotations = pd.read_excel(f"{_path}/raw_CellMarker2_annotations.xlsx")
    annotations = annotations[["cellontology_id", "UNIPROTID"]]
    annotations = annotations.dropna()
    annotations.cellontology_id = [
        i.replace("_", ":") for i in annotations.cellontology_id
    ]
    annotations.to_csv(f"{_path}/CellMarker2_annotations.csv")
    os.remove(f"{_path}/raw_CellMarker2_annotations.xlsx")


def _download_celltaxonomy_cell_types(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to cell types and CellTaxonomy.

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.
    """
    _path = f"{path}/cell_types"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["wrappers"]["resources_wrapper"]["CellTaxonomy"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            log.write(f"\t\t Downloading {file} from {url} ({_get_timestamp()})\n")
            _download_file(url, f"{_path}/{file}")
            log.write(f"\t\t\t DONE ({_get_timestamp()})")

    annotations = pd.read_csv(f"{_path}/raw_CellTaxonomy_annotations.txt", sep="\t")
    annotations = annotations.loc[annotations.Species == "Homo sapiens"]
    annotations = annotations[["Specific_Cell_Ontology_ID", "Uniprot"]]
    annotations = annotations.dropna()
    annotations.to_csv(f"{_path}/CellTaxonomy_annotations.csv")
    os.remove(f"{_path}/raw_CellTaxonomy_annotations.txt")


def _download_diseases(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to diseases (MeSH).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/diseases"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["wrappers"]["resources_wrapper"]["MeSH"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            log.write(f"\t\t Downloading {file} from {url} ({_get_timestamp()})\n")
            _download_file(url, f"{_path}/{file}")
            log.write(f"\t\t\t DONE ({_get_timestamp()})")

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
    """Download knowledge bases associated to gene groups (HGNC).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/gene_groups"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["wrappers"]["resources_wrapper"]["HGNC"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            log.write(f"\t\t Downloading {file} from {url} ({_get_timestamp()})\n")
            _download_file(url, f"{_path}/{file}")
            log.write(f"\t\t\t DONE ({_get_timestamp()})")


def _download_pathways(path: str, params: dict[str, Any]) -> None:
    """Download knowledge bases associated to pathways (Reactome).

    ``path`` is a path to store the downloaded files.
    ``params`` is a dict of parameters.

    """
    _path = f"{path}/pathways"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = params["wrappers"]["resources_wrapper"]["Reactome"]
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            log.write(f"\t\t Downloading {file} from {url} ({_get_timestamp()})\n")
            _download_file(url, f"{_path}/{file}")
            log.write(f"\t\t\t DONE ({_get_timestamp()})")