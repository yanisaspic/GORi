"""Functions called to load the knowledge base resources for the GORi analysis.

    2025/06/10 @yanisaspic"""

import os
import json
import pandas as pd
from typing import Any
from gori.params import get_parameters
from gori.src.utils import _get_uniprot_id


def get_available_resources() -> list[str]:
    """Get a list of resources that can be used in GORi.

    Returns
        A list of resource labels.
    """
    return [
        "GO_BP",
        "GO_CC",
        "GO_MF",
        "HGNC",
        "Reactome",
        "MeSH",
        "HPO",
        "CellMarker2_Abdomen",
        "CellMarker2_Adipose tissue",
        "CellMarker2_Airway",
        "CellMarker2_Amniotic fluid",
        "CellMarker2_Artery",
        "CellMarker2_Articulation",
        "CellMarker2_Belly",
        "CellMarker2_Bile duct",
        "CellMarker2_Biliary tract",
        "CellMarker2_Bladder",
        "CellMarker2_Blood",
        "CellMarker2_Blood vessel",
        "CellMarker2_Bone",
        "CellMarker2_Bone marrow",
        "CellMarker2_Brain",
        "CellMarker2_Breast",
        "CellMarker2_Bronchi",
        "CellMarker2_Bronchus",
        "CellMarker2_Decidua",
        "CellMarker2_Embryo",
        "CellMarker2_Endocrine organ",
        "CellMarker2_Endometrium",
        "CellMarker2_Epidermis",
        "CellMarker2_Epithelium",
        "CellMarker2_Esophagus",
        "CellMarker2_Eye",
        "CellMarker2_Fetal brain",
        "CellMarker2_Fetal liver",
        "CellMarker2_Fetal striatum",
        "CellMarker2_Fetus",
        "CellMarker2_Gall bladder",
        "CellMarker2_Gastrointestinal tract",
        "CellMarker2_Gut",
        "CellMarker2_Heart",
        "CellMarker2_Intervertebral disc",
        "CellMarker2_Intestine",
        "CellMarker2_Joint",
        "CellMarker2_Kidney",
        "CellMarker2_Knee",
        "CellMarker2_Knee joint",
        "CellMarker2_Limb",
        "CellMarker2_Liver",
        "CellMarker2_Lung",
        "CellMarker2_Lymph",
        "CellMarker2_Lymph node",
        "CellMarker2_Lymphoid tissue",
        "CellMarker2_Muscle",
        "CellMarker2_Nasal",
        "CellMarker2_Nasopharynx",
        "CellMarker2_Neck",
        "CellMarker2_Nodular tissue",
        "CellMarker2_Nose",
        "CellMarker2_Oral cavity",
        "CellMarker2_Ovary",
        "CellMarker2_Pancreas",
        "CellMarker2_Periodontium",
        "CellMarker2_Peritoneal fluid",
        "CellMarker2_Peritoneum",
        "CellMarker2_Pharynx",
        "CellMarker2_Placenta",
        "CellMarker2_Prostate",
        "CellMarker2_Salivary gland",
        "CellMarker2_Skeletal muscle",
        "CellMarker2_Skin",
        "CellMarker2_Soft tissue",
        "CellMarker2_Spinal cord",
        "CellMarker2_Spleen",
        "CellMarker2_Sputum",
        "CellMarker2_Stomach",
        "CellMarker2_Suprarenal gland",
        "CellMarker2_Synovial",
        "CellMarker2_Synovium",
        "CellMarker2_Taste bud",
        "CellMarker2_Tendon",
        "CellMarker2_Testis",
        "CellMarker2_Thorax",
        "CellMarker2_Thymus",
        "CellMarker2_Thyroid",
        "CellMarker2_Tonsil",
        "CellMarker2_Tooth",
        "CellMarker2_Trachea",
        "CellMarker2_Umbilical cord",
        "CellMarker2_Urine",
        "CellMarker2_Uterine cervix",
        "CellMarker2_Uterus",
        "CellMarker2_Vagina",
        "CellTaxonomy_Adrenal gland",
        "CellTaxonomy_Alimentary part of gastrointestinal system",
        "CellTaxonomy_Aorta wall",
        "CellTaxonomy_Ascitic fluid",
        "CellTaxonomy_Blood",
        "CellTaxonomy_Blood serum",
        "CellTaxonomy_Blood vessel",
        "CellTaxonomy_Bone marrow",
        "CellTaxonomy_Bone tissue",
        "CellTaxonomy_Brain",
        "CellTaxonomy_Breast",
        "CellTaxonomy_Cardiac muscle tissue",
        "CellTaxonomy_Carotid artery segment",
        "CellTaxonomy_Cerebrospinal fluid",
        "CellTaxonomy_Colon",
        "CellTaxonomy_Colonic mucosa",
        "CellTaxonomy_Colorectum",
        "CellTaxonomy_Common carotid artery plus branches",
        "CellTaxonomy_Duodenum",
        "CellTaxonomy_Endometrium",
        "CellTaxonomy_Entorhinal cortex",
        "CellTaxonomy_Esophagus",
        "CellTaxonomy_Femoral artery",
        "CellTaxonomy_Gall bladder",
        "CellTaxonomy_Gingiva",
        "CellTaxonomy_Heart",
        "CellTaxonomy_Ileum",
        "CellTaxonomy_Intestine",
        "CellTaxonomy_Kidney",
        "CellTaxonomy_Layer of synovial tissue",
        "CellTaxonomy_Liver",
        "CellTaxonomy_Lung",
        "CellTaxonomy_Lymph",
        "CellTaxonomy_Lymph node",
        "CellTaxonomy_Lymphoid tissue",
        "CellTaxonomy_Mammary gland",
        "CellTaxonomy_Meninx",
        "CellTaxonomy_Muscle organ",
        "CellTaxonomy_Nasal cavity mucosa",
        "CellTaxonomy_Nasopharynx",
        "CellTaxonomy_Neural tissue",
        "CellTaxonomy_Omentum",
        "CellTaxonomy_Oral cavity",
        "CellTaxonomy_Ovary",
        "CellTaxonomy_Pancreas",
        "CellTaxonomy_Peritoneal cavity",
        "CellTaxonomy_Peritoneum",
        "CellTaxonomy_Pleura",
        "CellTaxonomy_Pleural fluid",
        "CellTaxonomy_Prefrontal cortex",
        "CellTaxonomy_Prostate gland",
        "CellTaxonomy_Rectus femoris",
        "CellTaxonomy_Respiratory airway",
        "CellTaxonomy_Saliva-secreting gland",
        "CellTaxonomy_Skin epidermis",
        "CellTaxonomy_Skin of body",
        "CellTaxonomy_Stomach",
        "CellTaxonomy_Synovial fluid",
        "CellTaxonomy_Synovial joint",
        "CellTaxonomy_Testis",
        "CellTaxonomy_Thyroid gland",
        "CellTaxonomy_Tonsil",
        "CellTaxonomy_Urinary bladder",
        "CellTaxonomy_Uvea",
        "CellTaxonomy_Vascular system",
        "CellTaxonomy_White adipose tissue",
    ]


def load_resources(
    resources: set[str],
    path: str = "./resources",
    params: dict[str, Any] = get_parameters(),
) -> dict[str, Any]:
    """Load a subset of knowledge bases readily available in the GORi package.

    ``resources`` is a set of strings containing the names of the knowledge bases to load, including:
        `GO_BP` (biological processes),
        `GO_CC` (cellular components),
        `CellMarker2` (cell types with CellMarker 2.0 annotations),
        `CellTaxonomy` (cell types with CellTaxonomy annotations),
        `MeSH` (diseases),
        `HGNC` (gene groups),
        `GO_MF` (molecular functions),
        `Reactome` (pathways) and
        `HPO` (phenotypes).
    ``path`` is the path to the JSON files containing the knowledge bases.
    ``params`` is a dict of parameters.

    Returns
        A dict associating knowledge base names (keys) to their contents (values).
    """
    load_wrapper = params["wrappers"]["load_wrapper"]
    if any(p in resources for p in ["GO_BP", "GO_CC", "GO_MF"]):
        go = load_wrapper["GO_BP"](
            path
        )  # outputs the same GOAnnotation as GO_CC and GO_MF

    knowledge_bases = {}
    for p in resources:
        if p not in load_wrapper.keys():
            if p.startswith(("CellMarker2", "CellTaxonomy")):
                knowledge_bases[p] = load_wrapper["CellMarker2"](path, p)
            else:
                raise ValueError(
                    f"Invalid resource: {p}. Valid resources are: {get_available_resources()}"
                )
        elif p in ["GO_BP", "GO_CC", "GO_MF"]:
            knowledge_bases[p] = go
        else:
            knowledge_bases[p] = load_wrapper[p](path)
    return knowledge_bases


def load_feve(path: str, direction: str = "any") -> dict[str, dict[str, Any]]:
    """Load a custom knowledge base from a fEVE analysis.

    ``path`` is the path to the .xlsx file containing the results of a fEVE analysis.
    ``direction`` is the direction of the genes' expression: one of 'up', 'down' or 'any'.

    Returns
        A dict with three keys:
            ``annotations``: a dict associating a UniProt ID to its related clusters
            ``hierarchy``: a dict associating a cluster to its parent clusters in the hierarchy
            ``translations``: a placeholder for the cluster names
    """
    meta = pd.read_excel(path, sheet_name="meta", index_col=0)
    features = pd.read_excel(path, sheet_name="features", index_col=0)

    genes_wrapper = {
        "up": lambda x: x > 0,
        "down": lambda x: x < 0,
        "any": lambda x: x != 0,
    }
    f = genes_wrapper[direction]

    annotations = {
        _get_uniprot_id(gene): set(features.columns[f(features.loc[gene])])
        for gene in features.index
    }
    annotations = {gene: clusters for gene, clusters in annotations.items() if clusters}
    hierarchy = meta.iloc[1:].groupby(meta.iloc[1:].index).parent.apply(set).to_dict()
    translations = {clu: clu for clu in meta.index}

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }


def load_local(resource: str, path: str) -> dict[str, dict[str, Any]]:
    """Load a local knowledge base.

    ``resource`` is the label of the knowledge base to add.
    ``path`` is the path to the JSON files containing the annotations, the hierarchy
    and the translations of the knowledge base.

    Returns
        A dict with three keys:
            ``annotations``: a dict associating a Uniprot ID to its related concepts
            ``ontology``: a dict associating a concept to its parent concepts in the hierarchy
            ``translations``: a dict associating a concept to its human-readable label
    """
    with open(f"{path}/{resource}_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        _get_uniprot_id(gene): set(cids) for gene, cids in annotations.items()
    }

    with open(f"{path}/{resource}_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {cid: set(parents) for cid, parents in hierarchy.items()}

    if os.path.isfile(f"{path}/{resource}_t.json"):
        with open(f"{path}/{resource}_t.json", "r") as file:
            translations = json.load(file)
    else:
        concepts = {p for parents in hierarchy.values() for p in parents}
        concepts = concepts.union(hierarchy.keys())
        translations = {c: c for c in concepts}

    return {
        "annotations": annotations,
        "hierarchy": hierarchy,
        "translations": translations,
    }
