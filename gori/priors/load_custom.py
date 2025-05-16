"""Function called to load a custom prior for a GORi analysis.

    2025/04/30 @yanisaspic"""

import os
import json
from typing import Any
from gori.utils import _get_uniprot_id


def load_custom(prior: str, path: str) -> dict[str, dict[str, Any]]:
    """Load a custom knowledge base.

    ``prior`` is the name of the knowledge base to add.
    ``path`` is the path to the JSON files containing the annotations, the hierarchy
    and the translations of the knowledge base.

    Returns
        A dict with three keys:
            - "annotations": a dict associating a Uniprot ID to its related concepts
            - "ontology": a dict associating a concept to its parent concepts in the hierarchy
            - "translations": a dict associating a concept to its human-readable label
    """
    with open(f"{path}/{prior}_a.json", "r") as file:
        annotations = json.load(file)
    annotations = {
        _get_uniprot_id(gene): set(cids) for gene, cids in annotations.items()
    }

    with open(f"{path}/{prior}_h.json", "r") as file:
        hierarchy = json.load(file)
    hierarchy = {cid: set(parents) for cid, parents in hierarchy.items()}

    if os.path.isfile(f"{path}/{prior}_t.json"):
        with open(f"{path}/{prior}_t.json", "r") as file:
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
