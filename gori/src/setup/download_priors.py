"""Functions used by gori.init.download_priors().

    2025/05/20 @yanisaspic"""

import os
import gzip
import time
import datetime
import pandas as pd
from typing import Callable
from shutil import copyfileobj
from urllib.request import urlretrieve


def _get_timestamp() -> str:
    """Get a formatted timestamp.

    Returns
        A string with the date and the time."""
    timestamp = datetime.datetime.fromtimestamp(time.time()).strftime(
        "%Y/%m/%d_%H:%M:%S"
    )
    return timestamp


def _download_cell_types(path: str) -> None:
    """Download knowledge bases associated to cell types (CTYP).

    ``path`` is a path to store the downloaded files.
    """
    _path = f"{path}/cell_types"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = {
        "raw_CellTaxonomy_annotations.txt": "https://download.cncb.ac.cn/celltaxonomy/Cell_Taxonomy_resource.txt",
        "cell_types_ontology.obo": "https://purl.obolibrary.org/obo/cl/cl-basic.obo",
    }
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


def _download_diseases(path: str) -> None:
    """Download knowledge bases associated to diseases (DISE).

    ``path`` is a path to store the downloaded files.
    """
    _path = f"{path}/diseases"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = {
        "MeSH_hierarchy_C.bin": "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/c2025.bin",
        "MeSH_hierarchy_D.bin": "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2025.bin",
        "raw_CTD_annotations.csv.gz": "https://ctdbase.org/reports/CTD_curated_genes_diseases.csv.gz",
    }
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


def _download_gene_groups(path: str) -> None:
    """Download knowledge bases associated to gene groups (GENG).

    ``path`` is a path to store the downloaded files.
    """
    _path = f"{path}/gene_groups"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = {
        "HGNC_labels.csv": "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/family.csv",
        "HGNC_hierarchy.csv": "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/hierarchy.csv",
    }
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            urlretrieve(url, f"{_path}/{file}")
            log.write(f"\t\t {_get_timestamp()}: Downloaded {file} from {url}\n")


def _download_pathways(path: str) -> None:
    """Download knowledge bases associated to pathways (PATH).

    ``path`` is a path to store the downloaded files.
    """
    _path = f"{path}/pathways"
    if not os.path.isdir(_path):
        os.mkdir(_path)

    resources = {
        "Reactome_annotations.txt": "https://reactome.org/download/current/UniProt2Reactome.txt",
        "Reactome_hierarchy.txt": "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
        "Reactome_labels.txt": "https://reactome.org/download/current/ReactomePathways.txt",
    }
    with open("./downloads.log", "a") as log:
        for file, url in resources.items():
            urlretrieve(url, f"{_path}/{file}")
            log.write(f"\t\t {_get_timestamp()}: Downloaded {file} from {url}\n")
