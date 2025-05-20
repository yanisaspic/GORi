"""Functions used by gori.init.load_priors().

    2025/05/20 @yanisaspic"""

import json
import pandas as pd
from pypath.utils import go
from pypath.inputs import hpo, hgnc
from nxontology.imports import from_file
from typing import Any
from gori.src.utils import _prune_hierarchy, _get_uniprot_id


def _load_cell_types(path: str) -> dict[str, Any]:
    """Load the cell types knowledge base.

    ``path`` is the path to the JSON and OBO files containing the cell types knowledge base.

    Returns
        A dict with two keys:
            - "annotations": a dict associating a Uniprot ID to its associated cell types
            - "ontology": a graph associating a cell type to its parents in the hierarchy and its human-readable label.
    """
    ontology = from_file(f"{path}/CTYP_o.obo").graph
    current_cells = set(ontology.nodes)

    with open(f"{path}/CTYP_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        uid: set(cells).intersection(current_cells)
        for uid, cells in annotations.items()
    }

    return {"annotations": annotations, "ontology": ontology}


def _get_roots_diseases() -> set[str]:
    """Get the roots of the diseases hierarchy.

    Returns
        A set of MeSH IDs.
    """
    roots = {
        "D000820",  # Animal Diseases
        "D002318",  # Cardiovascular Diseases
        "D064419",  # Chemically-Induced Disorders
        "D009358",  # Congenital, Hereditary, and Neonatal Diseases and Abnormalities
        "D004066",  # Digestive System Diseases
        "D007280",  # Disorders of Environmental Origin
        "D004700",  # Endocrine System Diseases
        "D005128",  # Eye Diseases
        "D006425",  # Hemic and Lymphatic Diseases
        "D007154",  # Immune System Diseases
        "D007239",  # Infections
        "D009140",  # Musculoskeletal Diseases
        "D009369",  # Neoplasms
        "D009422",  # Nervous System Diseases
        "D009750",  # Nutritional and Metabolic Diseases
        "D009784",  # Occupational Diseases
        "D010038",  # Otorhinolaryngologic Diseases
        "D013568",  # Pathological Conditions, Signs and Symptoms
        "D012140",  # Respiratory Tract Diseases
        "D017437",  # Skin and Connective Tissue Diseases
        "D009057",  # Stomatognathic Diseases
        "D000091642",  # Urogenital Diseases
        "D014947",  # Wounds and Injuries
    }
    return roots


def _load_diseases(path: str) -> dict[str, dict[str, Any]]:
    """Load the diseases knowledge base.

    ``path`` is the path to the JSON files containing the diseases knowledge base.

    Returns
        A dict with three keys:
            - "annotations": a dict associating a Uniprot ID to its associated diseases
            - "hierarchy": a dict associating a MeSH ID to its parents in the hierarchy
            - "translations": a dict associating a MeSH ID to its human-readable label
    """
    with open(f"{path}/DISE_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        _get_uniprot_id(gene): set(diseases) for gene, diseases in annotations.items()
    }

    with open(f"{path}/DISE_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {meshid: set(parents) for meshid, parents in hierarchy.items()}
    hierarchy = _prune_hierarchy(hierarchy, roots=_get_roots_diseases())

    with open(f"{path}/DISE_t.json", "r") as file:
        translations = json.load(file)

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _load_gene_groups(path: str) -> dict[str, dict[str, Any]]:
    """Load the gene groups knowledge base.

    ``path`` is the path to the JSON files containing the gene groups knowledge base.

    Returns
        A dict with three keys:
            - "annotations": a dict associating a Uniprot ID to its associated gene groups
            - "hierarchy": a dict associating a gene group to its parents in the hierarchy
            - "translations": a dict associating a gene group to its human-readable label
    """
    annotations = hgnc.hgnc_genegroups()

    with open(f"{path}/GENG_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {gene_group: set(parents) for gene_group, parents in hierarchy.items()}

    with open(f"{path}/GENG_t.json", "r") as file:
        translations = json.load(file)

    return {
        "annotations": hgnc.hgnc_genegroups(),
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _load_pathways(path: str) -> dict[str, dict[str, Any]]:
    """Load the pathways knowledge base.

    ``path`` is the path to the JSON files containing the pathways knowledge base.

    Returns
        A dict with three keys:
            - "annotations": a dict associating a Uniprot ID to its associated pathways
            - "hierarchy": a dict associating a pathway to its parents in the hierarchy
            - "translations": a dict associating a pathway to its human-readable label
    """
    with open(f"{path}/PATH_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {uid: set(paths) for uid, paths in annotations.items()}

    with open(f"{path}/PATH_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {path: set(parents) for path, parents in hierarchy.items()}

    with open(f"{path}/PATH_t.json", "r") as file:
        translations = json.load(file)

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _get_root_phenotypes() -> set[str]:
    """Get the root of the phenotypes hierarchy.

    Returns
        A set of HPO IDs.
    """
    root = {"HP:0000118"}
    return root


def _load_phenotypes(path: None) -> dict[str, dict[str, Any]]:
    """Load the phenotypes knowledge base.

    ``path`` is a placeholder.

    Returns
        A dict with three keys:
            - "annotations": a dict associating a Uniprot ID to its associated phenotypes
            - "hierarchy": a dict associating a phenotype to its parents in the hierarchy
            - "translations": a dict associating a phenotype to its human-readable label
    """
    annotations = hpo.hpo_annotations()
    hierarchy = hpo.hpo_ontology()["parents"]
    hierarchy = _prune_hierarchy(hierarchy, roots=_get_root_phenotypes())
    translations = hpo.hpo_ontology()["terms"]
    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _load_gene_ontology() -> go.GOAnnotation:
    """Load the gene ontology knowledge base.

    Returns
        A GOAnnotation object containing the gene ontology.
    """
    return go.GOAnnotation()
