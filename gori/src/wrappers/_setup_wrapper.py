"""Functions used by setup_wrapper().

    2025/05/20 @yanisaspic"""

import os
import re
import json
import pandas as pd
from typing import Any, Callable
from gori.src.utils import _get_uniprot_id


def _setup_cell_types(dl_path: str, su_path: str) -> None:
    """Set-up the CellTaxonomy and Cell Ontology (CL) cell type annotations.

    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.

    """
    _dl_path = f"{dl_path}/cell_types"
    annotations = pd.read_csv(f"{_dl_path}/CellTaxonomy_annotations.csv")
    # split rows with multiple uniprot ids into multiple rows with a single id:
    annotations = annotations.assign(
        Uniprot=annotations.Uniprot.str.split(",")
    ).explode("Uniprot")
    annotations = annotations.groupby("Uniprot")["Specific_Cell_Ontology_ID"].apply(
        list
    )
    annotations = annotations.to_dict()
    with open(f"{su_path}/CTYP_a.json", "w") as file:
        json.dump(annotations, file)
    os.rename(f"{_dl_path}/cell_types_ontology.obo", f"{su_path}/CTYP_o.obo")


def _get_mesh_template(aspect: str) -> dict[str, str]:
    """Get the information required to parse a MeSH file.

    ``aspect`` is one of 'C' (supplementary concepts) or 'D' (main concepts).

    Returns
        A dict associating template flags (keys) to their values.
    """
    MeSH_TEMPLATE = {
        "RECORD_SEPARATOR": "\n\n\\*NEWRECORD",
        "ENTRY_SEPARATOR": " = ",
        "ID_FLAG": "UI =",
    }
    if aspect == "D":  # main concepts
        MeSH_TEMPLATE["LABEL_FLAG"] = "MH ="
        MeSH_TEMPLATE["TREE_FLAG"] = "MN"
    if aspect == "C":  # supplementary concepts
        MeSH_TEMPLATE["LABEL_FLAG"] = "NM ="
        MeSH_TEMPLATE["TREE_FLAG"] = "HM ="
    return MeSH_TEMPLATE


def _parse_mesh_block(block: str, aspect: str) -> dict[str, dict[str, Any]]:
    """Parse a MeSH record element.

    ``block`` is a record element from a MeSH ASCII file.
    ``aspect`` is one of 'C' (supplementary concepts) or 'D' (main concepts).

    Returns
        A dict associating a MeSH ID to its human-readable label and its IDs in the hierarchy tree.
    """
    template = _get_mesh_template(aspect)
    data = block.split("\n")
    flags_of_interest = [
        template["ID_FLAG"],
        template["LABEL_FLAG"],
        template["TREE_FLAG"],
    ]
    parse_entry = lambda entry: entry.split(template["ENTRY_SEPARATOR"])[1]
    entries = [
        parse_entry(entry)
        for entry in data
        if any(flag in entry for flag in flags_of_interest)
    ]
    return {entries[-1]: {"label": entries[0], "tree_ids": entries[1:-1]}}


def _parse_mesh_file(filename: str, aspect: str) -> dict[str, dict[str, Any]]:
    """Parse a MeSH ASCII file.

    ``filename`` is a path to a MeSH ASCII file.
    ``aspect`` is one of 'C' (supplementary concepts) or 'D' (main concepts).

    Returns
        A dict associating MeSH IDs to their parents and their human-readable labels.
    """
    template = _get_mesh_template(aspect)
    with open(filename, "r") as file:
        mesh = file.read()
    flag = re.compile(template["RECORD_SEPARATOR"])
    mesh_blocks = flag.split(mesh)
    mesh_records = [_parse_mesh_block(block, aspect) for block in mesh_blocks]

    mesh_label_translations = {
        meshid: data["label"] for rec in mesh_records for meshid, data in rec.items()
    }
    mesh_tree_translations = {
        tree_id: meshid
        for rec in mesh_records
        for meshid, data in rec.items()
        for tree_id in data["tree_ids"]
    }

    if aspect == "D":  # main concepts
        get_parent_ids = lambda tree_ids: {  # e.g. tree_ids: D02.355.291.933.125
            mesh_tree_translations.get(tid.rsplit(".", 1)[0])
            for tid in tree_ids
            if "." in tid
        }
    if aspect == "C":  # supplementary concepts
        get_parent_ids = lambda tree_ids: {  # e.g. tree_ids: *D001561-Benzilates
            tid.split("-", 1)[0][1:] for tid in tree_ids
        }
    mesh_hierarchy = {
        meshid: get_parent_ids(data["tree_ids"])
        for rec in mesh_records
        for meshid, data in rec.items()
    }

    return {"translations": mesh_label_translations, "hierarchy": mesh_hierarchy}


def _setup_diseases(dl_path: str, su_path: str) -> None:
    """Set-up the Comparative Toxicogenomics Database (CTD) and Medical Subject Headings (MeSH) disease annotations and hierarchy.

    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    """
    _dl_path = f"{dl_path}/diseases"
    annotations = pd.read_csv(f"{_dl_path}/CTD_annotations.csv", header=None)
    annotations.columns = [
        "GeneSymbol",
        "GeneID",
        "DiseaseName",
        "DiseaseID",
        "DirectEvidence",
        "OmimIDs",
        "PubMedIDs",
    ]
    annotations = annotations[["GeneSymbol", "DiseaseID"]]
    trim_mesh_id = lambda mesh_id: mesh_id.split(":")[
        1
    ]  # remove the 'MESH': prefix to align with MeSH hierarchy

    annotations.DiseaseID = annotations.DiseaseID.apply(trim_mesh_id)
    annotations = annotations.groupby("GeneSymbol")["DiseaseID"].apply(list)
    annotations = annotations.to_dict()
    uniprot_annotations = {
        _get_uniprot_id(gs): mesh_ids for gs, mesh_ids in annotations.items()
    }
    with open(f"{su_path}/DISE_a.json", "w") as file:
        json.dump(uniprot_annotations, file)

    mesh_d = _parse_mesh_file(f"{_dl_path}/MeSH_hierarchy_D.bin", "D")
    mesh_c = _parse_mesh_file(f"{_dl_path}/MeSH_hierarchy_C.bin", "C")
    mesh = {
        "translations": mesh_d["translations"] | mesh_c["translations"],
        "hierarchy": mesh_d["hierarchy"] | mesh_c["hierarchy"],
    }
    with open(f"{su_path}/DISE_t.json", "w") as file:
        json.dump(mesh["translations"], file)
    hierarchy = {meshid: list(parents) for meshid, parents in mesh["hierarchy"].items()}
    with open(f"{su_path}/DISE_h.json", "w") as file:
        json.dump(hierarchy, file)


def _setup_gene_groups(dl_path: str, su_path: str) -> None:
    """Set-up the HGNC gene group hierarchy and translations.

    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    """
    _dl_path = f"{dl_path}/gene_groups"
    hierarchy = pd.read_csv(f"{_dl_path}/HGNC_hierarchy.csv")
    labels = pd.read_csv(f"{_dl_path}/HGNC_labels.csv", index_col=0)
    labels.name = labels.name.apply(lambda n: n.strip())
    translations = labels.name.to_dict()
    with open(f"{su_path}/GENG_t.json", "w") as file:
        json.dump(translations, file)

    get_family_name = lambda id: translations[id]
    hierarchy = hierarchy.applymap(get_family_name)
    hierarchy = hierarchy.groupby("child_fam_id")["parent_fam_id"].apply(list)
    hierarchy = hierarchy.to_dict()
    with open(f"{su_path}/GENG_h.json", "w") as file:
        json.dump(hierarchy, file)


def _setup_pathways(dl_path: str, su_path: str) -> None:
    """Set-up the Reactome pathway annotations and hierarchy.

    ``dl_path`` is a path where downloaded files are stored.
    ``su_path`` is a path to store the set-up files.
    """
    _dl_path = f"{dl_path}/pathways"
    labels = pd.read_csv(
        f"{_dl_path}/Reactome_labels.txt", sep="\t", index_col=0, header=None
    )
    labels.columns = ["name", "species"]
    hierarchy = pd.read_csv(f"{_dl_path}/Reactome_hierarchy.txt", sep="\t", header=None)
    hierarchy.columns = ["parent_reactome_id", "child_reactome_id"]
    annotations = pd.read_csv(
        f"{_dl_path}/Reactome_annotations.txt", sep="\t", header=None
    )
    annotations.columns = [
        "uniprot_id",
        "reactome_id",
        "link",
        "name",
        "curation",
        "species",
    ]

    labels = labels.loc[labels.species == "Homo sapiens"]
    translations = labels.name.to_dict()
    with open(f"{su_path}/PATH_t.json", "w") as file:
        json.dump(translations, file)

    _annotations = annotations.loc[annotations.species == "Homo sapiens"]
    annotations = _annotations.groupby("uniprot_id")["reactome_id"].apply(list)
    annotations = annotations.to_dict()
    with open(f"{su_path}/PATH_a.json", "w") as file:
        json.dump(annotations, file)

    hierarchy = hierarchy.loc[
        hierarchy.child_reactome_id.isin(_annotations.reactome_id)
    ]
    hierarchy = hierarchy.groupby("child_reactome_id")["parent_reactome_id"].apply(list)
    hierarchy = hierarchy.to_dict()
    with open(f"{su_path}/PATH_h.json", "w") as file:
        json.dump(hierarchy, file)
