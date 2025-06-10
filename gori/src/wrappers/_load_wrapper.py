"""Functions used by load_wrapper().

    2025/05/22 @yanisaspic"""

import json
from typing import Any, Optional
from pypath.inputs import hpo, hgnc
from pypath.utils.go import GOAnnotation
from nxontology.imports import from_file
from gori.src.utils import _prune_hierarchy, _get_uniprot_id


def _load_cell_types(path: str, label: str) -> dict[str, Any]:
    """Load the cell types knowledge base.

    ``path`` is the path to the JSON and OBO files containing the cell types knowledge base.
    ``label`` is one of CellMarker2 or CellTaxonomy.

    Returns
        A dict with two keys:
            `annotations`: a dict associating a Uniprot ID to its associated cell types
            `ontology`: a graph associating a cell type to its parents in the hierarchy and its human-readable label.
    """
    ontology = from_file(f"{path}/CellMarker2_o.obo").graph
    current_cells = set(ontology.nodes)

    with open(f"{path}/{label}_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        uid: set(cells).intersection(current_cells)
        for uid, cells in annotations.items()
    }
    return {"annotations": annotations, "ontology": ontology}


def _load_cellmarker2_cell_types(path: str) -> dict[str, Any]:
    """Load the cell types knowledge base with CellMarker 2.0 annotations.

    ``path`` is the path to the JSON and OBO files containing the cell types knowledge base.

    Returns
    A dict with two keys:
        `annotations`: a dict associating a Uniprot ID to its associated cell types
        `ontology`: a graph associating a cell type to its parents in the hierarchy and its human-readable label.
    """
    return _load_cell_types(path, "CellMarker2")


def _load_celltaxonomy_cell_types(path: str) -> dict[str, Any]:
    """Load the cell types knowledge base with CellMarker 2.0 annotations.

    ``path`` is the path to the JSON and OBO files containing the cell types knowledge base.

    Returns
    A dict with two keys:
        `annotations`: a dict associating a Uniprot ID to its associated cell types
        `ontology`: a graph associating a cell type to its parents in the hierarchy and its human-readable label.
    """
    return _load_cell_types(path, "CellTaxonomy")


def _get_roots_diseases() -> set[str]:
    """Get the roots of the diseases hierarchy.

    Returns
        A set of MeSH IDs.
    """
    return {
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


def _load_diseases(path: str) -> dict[str, dict[str, Any]]:
    """Load the diseases knowledge base.

    ``path`` is the path to the JSON files containing the diseases knowledge base.

    Returns
        A dict with three keys:
            `annotations`: a dict associating a Uniprot ID to its associated diseases
            `hierarchy`: a dict associating a MeSH ID to its parents in the hierarchy
            `translations`: a dict associating a MeSH ID to its human-readable label
    """
    with open(f"{path}/MeSH_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        _get_uniprot_id(gene): set(diseases) for gene, diseases in annotations.items()
    }

    with open(f"{path}/MeSH_h.json", "r") as file:
        hierarchy = json.load(file)
    tmp = {meshid: set(parents) for meshid, parents in hierarchy.items()}
    roots = _get_roots_diseases()
    hierarchy = _prune_hierarchy({"hierarchy": tmp}, roots=roots)

    with open(f"{path}/MeSH_t.json", "r") as file:
        translations = json.load(file)
        translations = {
            k: v for k, v in translations.items() if k in hierarchy.keys()
        }  # filter out non-diseases (e.g. molecules)

    # assign a meta-root to the hierarchy:
    for r in roots:
        hierarchy[r] = {"disease"}
    translations["disease"] = "disease"

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
            `annotations`: a dict associating a Uniprot ID to its associated gene groups
            `hierarchy`: a dict associating a gene group to its parents in the hierarchy
            `translations`: a dict associating a gene group to its human-readable label
    """
    annotations = hgnc.hgnc_genegroups()

    with open(f"{path}/HGNC_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {gene_group: set(parents) for gene_group, parents in hierarchy.items()}

    with open(f"{path}/HGNC_t.json", "r") as file:
        translations = json.load(file)

    # assign a meta-root to the hierarchy:
    roots = set(translations.values()).difference(set(hierarchy.keys()))
    for r in roots:
        hierarchy[r] = {"gene_group"}
    roots = set(translations.values()).difference(set(hierarchy.keys()))
    translations["gene_group"] = "gene_group"

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _load_pathways(path: str) -> dict[str, dict[str, Any]]:
    """Load the pathways knowledge base.

    ``path`` is the path to the JSON files containing the pathways knowledge base.

    Returns
        A dict with three keys:
            `annotations`: a dict associating a Uniprot ID to its associated pathways
            `hierarchy`: a dict associating a pathway to its parents in the hierarchy
            `translations`: a dict associating a pathway to its human-readable label
    """
    with open(f"{path}/Reactome_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {uid: set(paths) for uid, paths in annotations.items()}

    with open(f"{path}/Reactome_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {path: set(parents) for path, parents in hierarchy.items()}

    with open(f"{path}/Reactome_t.json", "r") as file:
        translations = json.load(file)

    # assign a meta-root to the hierarchy:
    roots = set(translations.keys()).difference(set(hierarchy.keys()))
    for r in roots:
        hierarchy[r] = {"pathway"}
    translations["pathway"] = "pathway"

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _get_roots_phenotypes() -> set[str]:
    """Get the roots of the phenotypes hierarchy.

    Returns
        A set of HPO IDs.
    """
    return {"HP:0000118"}


def _load_phenotypes(path: Optional[str]) -> dict[str, dict[str, Any]]:
    """Load the phenotypes knowledge base.

    ``path`` is a placeholder.

    Returns
        A dict with three keys:
            `annotations`: a dict associating a Uniprot ID to its associated phenotypes
            `hierarchy`: a dict associating a phenotype to its parents in the hierarchy
            `translations`: a dict associating a phenotype to its human-readable label
    """
    roots = _get_roots_phenotypes()
    annotations = hpo.hpo_annotations()
    tmp = hpo.hpo_ontology()["parents"]
    hierarchy = _prune_hierarchy({"hierarchy": tmp}, roots=roots)
    for r in roots:
        hierarchy[r] = set()
    translations = hpo.hpo_ontology()["terms"]
    translations = {
        k: v for k, v in translations.items() if k in hierarchy.keys()
    }  # filter out non-phenotypes (e.g. onsets)
    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def _load_gene_ontology(path: Optional[str]) -> GOAnnotation:
    """Load the Gene Ontology knowledge base.

    ``path`` is a placeholder.

    Returns
        A GOAnnotation object.
    """
    return GOAnnotation()
