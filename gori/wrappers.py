"""Wrapper functions used during a GORi analysis.

    2025/05/22 @yanisaspic"""

from typing import Callable
from gori.src.wrappers._ancestors_wrapper import _get_ctyp_ancestors, _get_go_ancestors
from gori.src.wrappers._annotations_wrapper import (
    _get_biop_annotations,
    _get_celc_annotations,
    _get_ctyp_annotations,
    _get_dise_annotations,
    _get_geng_annotations,
    _get_molf_annotations,
    _get_phen_annotations,
)
from gori.src.wrappers._descendants_wrapper import (
    _get_ctyp_descendants,
    _get_go_descendants,
)
from gori.src.wrappers._download_wrapper import (
    _download_cellmarker2_cell_types,
    _download_celltaxonomy_cell_types,
    _download_diseases,
    _download_gene_groups,
    _download_pathways,
)
from gori.src.wrappers._inverse_translate_wrapper import (
    _get_ctyp_inverse_translation,
    _get_geng_inverse_translation,
    _get_go_inverse_translation,
)
from gori.src.wrappers._load_wrapper import (
    _load_cell_types,
    _load_diseases,
    _load_gene_groups,
    _load_gene_ontology,
    _load_pathways,
    _load_phenotypes,
)
from gori.src.wrappers._roots_wrapper import (
    _get_biop_roots,
    _get_celc_roots,
    _get_ctyp_roots,
    _get_molf_roots,
)
from gori.src.wrappers._setup_wrapper import (
    _setup_cellmarker2_cell_types,
    _setup_celltaxonomy_cell_types,
    _setup_diseases,
    _setup_gene_groups,
    _setup_pathways,
)
from gori.src.wrappers._terms_wrapper import (
    _get_biop_terms,
    _get_celc_terms,
    _get_ctyp_terms,
    _get_geng_terms,
    _get_molf_terms,
)
from gori.src.wrappers._translate_wrapper import (
    _get_ctyp_translation,
    _get_geng_translation,
    _get_go_translation,
)
from gori.src.wrappers._urls_wrapper import (
    _get_ctyp_url,
    _get_geng_url,
)


def ancestors_wrapper() -> dict[str, Callable]:
    """A wrapper to get the ancestors of an annotation term.

    It returns a resource-specific function, which expects two parameters:
        ``terms`` is a set of annotation terms. (set[str])
        ``resource`` is the contents of a knowledge base.

    Returns
        A dict associating resources (keys) to their ancestors function.
    """
    return {
        "GO_BP": _get_go_ancestors,
        "GO_CC": _get_go_ancestors,
        "CellMarker2": _get_ctyp_ancestors,  # CellTaxonomy uses the same function
        "GO_MF": _get_go_ancestors,
    }


def annotations_wrapper() -> dict[str, Callable]:
    """A wrapper to get the annotations of a gene.

    It returns a resource-specific function, which expects four parameters:
        ``uid`` is a UniProtID. (str)
        ``resource`` is a resource label. (str)
        ``data`` is a dict associating resources (keys) to their contents (values). (dict[str, Any])
        ``params`` is a dict of parameters. (dict[str, Any])

    Returns
        A dict associating resources (keys) to their annotations function.
    """
    return {
        "GO_BP": _get_biop_annotations,
        "GO_CC": _get_celc_annotations,
        "CellMarker2": _get_ctyp_annotations,  # CellTaxonomy uses the same function
        "MeSH": _get_dise_annotations,
        "HGNC": _get_geng_annotations,
        "GO_MF": _get_molf_annotations,
        "HPO": _get_phen_annotations,
    }


def descendants_wrapper() -> dict[str, Callable]:
    """A wrapper to get the descendants of an annotation term.

    It returns a resource-specific function, which expects two parameters:
        ``terms`` is a set of annotation terms. (set[str])
        ``resource`` is the contents of a knowledge base.

    Returns
        A dict associating resources (keys) to their descendants function.
    """
    return {
        "GO_BP": _get_go_descendants,
        "GO_CC": _get_go_descendants,
        "CellMarker2": _get_ctyp_descendants,  # CellTaxonomy uses the same function
        "GO_MF": _get_go_descendants,
    }


def download_wrapper() -> dict[str, Callable]:
    """A wrapper to download curated knowledge bases (i.e. resources).

    It returns a resource-specific function, which expects two parameters:
        ``path`` is a path to store the downloaded files. (str)
        ``params`` is a dict of parameters. (dict[str, Any])

    Returns
        A dict associating resources (keys) to a download function (values).
    """
    return {
        "CellMarker2": _download_cellmarker2_cell_types,
        "CellTaxonomy": _download_celltaxonomy_cell_types,
        "MeSH": _download_diseases,
        "HGNC": _download_gene_groups,
        "Reactome": _download_pathways,
    }


def headers_wrapper() -> dict[str, str]:
    """A wrapper with url headers to knowledge bases.

    Returns
        A dict associating resources (keys) to their url headers (values).
    """
    s = "%252F"  # slash character \
    headers = {
        "GO_BP": "www.ebi.ac.uk/QuickGO/term/",
        "GO_CC": "www.ebi.ac.uk/QuickGO/term/",
        "CellMarker2": f"www.ebi.ac.uk/ols4/ontologies/cl/classes/http:{s}{s}purl.obolibrary.org{s}obo{s}",  # CellTaxonomy uses the same header
        "MeSH": "meshb.nlm.nih.gov/record/ui?ui=",
        "HGNC": "www.genenames.org/data/genegroup/#!/group/",
        "GO_MF": "www.ebi.ac.uk/QuickGO/term/",
        "Reactome": "reactome.org/content/detail/",
        "HPO": "hpo.jax.org/browse/term/",
    }
    return headers


def inverse_translate_wrapper() -> dict[str, Callable]:
    """A wrapper to get the term associated to a human-readable label.

    It returns a resource-specific function, which expects two parameters:
        ``label`` is a a human-readable label. (str)
        ``data`` is a dict associating resources (keys) to their contents (values). (dict[str, dict[str, Any]])

    Returns
        A dict associating resources (keys) to an inverse_translate function.
    """
    return {  # Missing resources are translated with a generic function in _get_resource_inverse_translation()
        "GO_BP": _get_go_inverse_translation,
        "GO_CC": _get_go_inverse_translation,
        "CellMarker2": _get_ctyp_inverse_translation,  # CellTaxonomy uses the same function
        "HGNC": _get_geng_inverse_translation,
        "GO_MF": _get_go_inverse_translation,
    }


def load_wrapper() -> dict[str, Callable]:
    """A wrapper to load curated knowledge bases (i.e. resources).

    It returns a resource-specific function, which expects one parameter:
        ``path`` is the path to the JSON files containing the knowledge base, or a placeholder. (str)

    Returns
        A dict associating resources (keys) to a load function (values).
    """
    return {  # GO resources are included to generate a valid error message.
        "GO_BP": _load_gene_ontology,
        "GO_CC": _load_gene_ontology,
        "CellMarker2": _load_cell_types,  # CellTaxonomy uses the same function
        "MeSH": _load_diseases,
        "HGNC": _load_gene_groups,
        "GO_MF": _load_gene_ontology,
        "Reactome": _load_pathways,
        "HPO": _load_phenotypes,
    }


def resources_wrapper() -> dict[str, dict[str, str]]:
    """A wrapper with download links for curated knowledge bases (i.e. resources).

    Returns
        A dict associating resource labels (keys) to their download links (values).
    """
    wrapper = {
        "CellMarker2": {
            "raw_CellMarker2_annotations.xlsx": "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx",
            "cell_types_ontology.obo": "https://purl.obolibrary.org/obo/cl/cl-basic.obo",
        },
        "CellTaxonomy": {
            "raw_CellTaxonomy_annotations.txt": "https://download.cncb.ac.cn/celltaxonomy/Cell_Taxonomy_resource.txt",
            "cell_types_ontology.obo": "https://purl.obolibrary.org/obo/cl/cl-basic.obo",
        },
        "MeSH": {
            "MeSH_hierarchy_C.bin": "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/c2025.bin",
            "MeSH_hierarchy_D.bin": "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2025.bin",
            "raw_CTD_annotations.csv.gz": "https://ctdbase.org/reports/CTD_curated_genes_diseases.csv.gz",
        },
        "HGNC": {
            "HGNC_labels.csv": "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/family.csv",
            "HGNC_hierarchy.csv": "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/hierarchy.csv",
        },
        "Reactome": {
            "Reactome_annotations.txt": "https://reactome.org/download/current/UniProt2Reactome.txt",
            "Reactome_hierarchy.txt": "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
            "Reactome_labels.txt": "https://reactome.org/download/current/ReactomePathways.txt",
        },
    }
    return wrapper


def roots_wrapper() -> dict[str, Callable]:
    """A wrapper to get the roots of a resource.

    Returns
        A dict associating resources (keys) to their roots (values).
    """
    return {
        "GO_BP": _get_biop_roots,
        "GO_CC": _get_celc_roots,
        "CellMarker2": _get_ctyp_roots,  # CellTaxonomy uses the same function
        "GO_MF": _get_molf_roots,
    }


def setup_wrapper() -> dict[str, Callable]:
    """A wrapper to set-up curated knowledge bases (i.e. resources).

    It returns a resource-specific function, which expects two parameters:
        ``dl_path`` is a path where downloaded files are stored. (str)
        ``su_path`` is a path to store the set-up files. (str)

    Returns
        A dict associating resources (keys) to a download function (values).
    """
    return {
        "CellMarker2": _setup_cellmarker2_cell_types,
        "CellTaxonomy": _setup_celltaxonomy_cell_types,
        "MeSH": _setup_diseases,
        "HGNC": _setup_gene_groups,
        "Reactome": _setup_pathways,
    }


def terms_wrapper() -> dict[str, Callable]:
    """A wrapper to get every terms in a resource.

    It returns a resource specific function, which expects one parameter:
        ``resource`` is a knowledge base. (Any)

    Returns
        A dict associating resources (keys) to their term functions.
    """
    return {
        "GO_BP": _get_biop_terms,
        "GO_CC": _get_celc_terms,
        "CellMarker2": _get_ctyp_terms,  # CellTaxonomy uses the same function
        "HGNC": _get_geng_terms,
        "GO_MF": _get_molf_terms,
    }


def translate_wrapper() -> dict[str, Callable]:
    """A wrapper to get the human-readable label of a term.

    It returns a resource-specific function, which expects two parameters:
        ``term`` is an annotation term. (str)
        ``data`` is a dict associating resources (keys) to their contents (values). (dict[str, dict[str, Any]])

    Returns
        A dict associating resources (keys) to a translate function.
    """
    return {  # Missing resources are translated with a generic function in _get_resource_translation()
        "GO_BP": _get_go_translation,
        "GO_CC": _get_go_translation,
        "CellMarker2": _get_ctyp_translation,  # CellTaxonomy uses the same function
        "HGNC": _get_geng_translation,
        "GO_MF": _get_go_translation,
    }


def urls_wrapper() -> dict[str, Callable]:
    """A wrapper to get the url of an annotation term.

    It returns a resource-specific function, which expects two parameters:
        ``term`` is an annotation term. (str)
        ``header`` is a url header. (str)
        ``data`` is a dict associating resources (keys) to their contents (values). (dict[str, Any])

    Returns
        A dict associating resources (keys) to their url functions (values).
    """
    return {  # Missing resources' urls are generated with a generic function in _get_resource_url()
        "CellMarker2": _get_ctyp_url,  # CellTaxonomy uses the same function
        "HGNC": _get_geng_url,
    }
