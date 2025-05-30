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
    _download_cell_types,
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
    _setup_cell_types,
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

    It returns a prior-specific function, which expects two parameters:
        ``terms`` is a set of annotation terms. (set[str])
        ``prior`` is the contents of a knowledge base.

    Returns
        A dict associating priors (keys) to their ancestors function.
    """
    return {
        "BIOP": _get_go_ancestors,
        "CELC": _get_go_ancestors,
        "CTYP": _get_ctyp_ancestors,
        "MOLF": _get_go_ancestors,
    }


def annotations_wrapper() -> dict[str, Callable]:
    """A wrapper to get the annotations of a gene.

    It returns a prior-specific function, which expects four parameters:
        ``uid`` is a UniProtID. (str)
        ``prior`` is a prior label. (str)
        ``data`` is a dict associating priors (keys) to their contents (values). (dict[str, Any])
        ``params`` is a dict of parameters. (dict[str, Any])

    Returns
        A dict associating priors (keys) to their annotations function.
    """
    return {
        "BIOP": _get_biop_annotations,
        "CELC": _get_celc_annotations,
        "CTYP": _get_ctyp_annotations,
        "DISE": _get_dise_annotations,
        "GENG": _get_geng_annotations,
        "MOLF": _get_molf_annotations,
        "PHEN": _get_phen_annotations,
    }


def descendants_wrapper() -> dict[str, Callable]:
    """A wrapper to get the descendants of an annotation term.

    It returns a prior-specific function, which expects two parameters:
        ``terms`` is a set of annotation terms. (set[str])
        ``prior`` is the contents of a knowledge base.

    Returns
        A dict associating priors (keys) to their descendants function.
    """
    return {
        "BIOP": _get_go_descendants,
        "CELC": _get_go_descendants,
        "CTYP": _get_ctyp_descendants,
        "MOLF": _get_go_descendants,
    }


def download_wrapper() -> dict[str, Callable]:
    """A wrapper to download curated knowledge bases (i.e. priors).

    It returns a prior-specific function, which expects two parameters:
        ``path`` is a path to store the downloaded files. (str)
        ``params`` is a dict of parameters. (dict[str, Any])

    Returns
        A dict associating priors (keys) to a download function (values).
    """
    return {
        "CTYP": _download_cell_types,
        "DISE": _download_diseases,
        "GENG": _download_gene_groups,
        "PATH": _download_pathways,
    }


def headers_wrapper() -> dict[str, str]:
    """A wrapper with url headers to knowledge bases.

    Returns
        A dict associating priors (keys) to their url headers (values).
    """
    s = "%252F"  # slash character \
    headers = {
        "BIOP": "www.ebi.ac.uk/QuickGO/term/",
        "CELC": "www.ebi.ac.uk/QuickGO/term/",
        "CTYP": f"www.ebi.ac.uk/ols4/ontologies/cl/classes/http:{s}{s}purl.obolibrary.org{s}obo{s}",
        "DISE": "meshb.nlm.nih.gov/record/ui?ui=",
        "GENG": "www.genenames.org/data/genegroup/#!/group/",
        "MOLF": "www.ebi.ac.uk/QuickGO/term/",
        "PATH": "reactome.org/content/detail/",
        "PHEN": "hpo.jax.org/browse/term/",
    }
    return headers


def inverse_translate_wrapper() -> dict[str, Callable]:
    """A wrapper to get the term associated to a human-readable label.

    It returns a prior-specific function, which expects two parameters:
        ``label`` is a a human-readable label. (str)
        ``data`` is a dict associating priors (keys) to their contents (values). (dict[str, dict[str, Any]])

    Returns
        A dict associating priors (keys) to an inverse_translate function.
    """
    return {  # Missing priors are translated with a generic function in _get_prior_inverse_translation()
        "BIOP": _get_go_inverse_translation,
        "CELC": _get_go_inverse_translation,
        "CTYP": _get_ctyp_inverse_translation,
        "GENG": _get_geng_inverse_translation,
        "MOLF": _get_go_inverse_translation,
    }


def load_wrapper() -> dict[str, Callable]:
    """A wrapper to load curated knowledge bases (i.e. priors).

    It returns a prior-specific function, which expects one parameter:
        ``path`` is the path to the JSON files containing the knowledge base, or a placeholder. (str)

    Returns
        A dict associating priors (keys) to a load function (values).
    """
    return {  # GO priors are included to generate a valid error message.
        "BIOP": _load_gene_ontology,
        "CELC": _load_gene_ontology,
        "CTYP": _load_cell_types,
        "DISE": _load_diseases,
        "GENG": _load_gene_groups,
        "MOLF": _load_gene_ontology,
        "PATH": _load_pathways,
        "PHEN": _load_phenotypes,
    }


def resources_wrapper() -> dict[str, dict[str, str]]:
    """A wrapper with download links for curated knowledge bases (i.e. priors).

    Returns
        A dict associating prior labels (keys) to their download links (values).
    """
    wrapper = {
        "CTYP": {
            "raw_CellMarker2_annotations.xlsx": "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx",
            "cell_types_ontology.obo": "https://purl.obolibrary.org/obo/cl/cl-basic.obo",
        },
        "DISE": {
            "MeSH_hierarchy_C.bin": "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/c2025.bin",
            "MeSH_hierarchy_D.bin": "https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2025.bin",
            "raw_CTD_annotations.csv.gz": "https://ctdbase.org/reports/CTD_curated_genes_diseases.csv.gz",
        },
        "GENG": {
            "HGNC_labels.csv": "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/family.csv",
            "HGNC_hierarchy.csv": "https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/hierarchy.csv",
        },
        "PATH": {
            "Reactome_annotations.txt": "https://reactome.org/download/current/UniProt2Reactome.txt",
            "Reactome_hierarchy.txt": "https://reactome.org/download/current/ReactomePathwaysRelation.txt",
            "Reactome_labels.txt": "https://reactome.org/download/current/ReactomePathways.txt",
        },
    }
    return wrapper


def roots_wrapper() -> dict[str, Callable]:
    """A wrapper to get the roots of a prior.

    Returns
        A dict associating priors (keys) to their roots (values).
    """
    return {
        "BIOP": _get_biop_roots,
        "CELC": _get_celc_roots,
        "CTYP": _get_ctyp_roots,
        "MOLF": _get_molf_roots,
    }


def setup_wrapper() -> dict[str, Callable]:
    """A wrapper to set-up curated knowledge bases (i.e. priors).

    It returns a prior-specific function, which expects two parameters:
        ``dl_path`` is a path where downloaded files are stored. (str)
        ``su_path`` is a path to store the set-up files. (str)

    Returns
        A dict associating priors (keys) to a download function (values).
    """
    return {
        "CTYP": _setup_cell_types,
        "DISE": _setup_diseases,
        "GENG": _setup_gene_groups,
        "PATH": _setup_pathways,
    }


def terms_wrapper() -> dict[str, Callable]:
    """A wrapper to get every terms in a prior.

    It returns a prior specific function, which expects one parameter:
        ``prior`` is a knowledge base. (Any)

    Returns
        A dict associating priors (keys) to their term functions.
    """
    return {
        "BIOP": _get_biop_terms,
        "CELC": _get_celc_terms,
        "CTYP": _get_ctyp_terms,
        "GENG": _get_geng_terms,
        "MOLF": _get_molf_terms,
    }


def translate_wrapper() -> dict[str, Callable]:
    """A wrapper to get the human-readable label of a term.

    It returns a prior-specific function, which expects two parameters:
        ``term`` is an annotation term. (str)
        ``data`` is a dict associating priors (keys) to their contents (values). (dict[str, dict[str, Any]])

    Returns
        A dict associating priors (keys) to a translate function.
    """
    return {  # Missing priors are translated with a generic function in _get_prior_translation()
        "BIOP": _get_go_translation,
        "CELC": _get_go_translation,
        "CTYP": _get_ctyp_translation,
        "GENG": _get_geng_translation,
        "MOLF": _get_go_translation,
    }


def urls_wrapper() -> dict[str, Callable]:
    """A wrapper to get the url of an annotation term.

    It returns a prior-specific function, which expects two parameters:
        ``term`` is an annotation term. (str)
        ``header`` is a url header. (str)
        ``data`` is a dict associating priors (keys) to their contents (values). (dict[str, Any])

    Returns
        A dict associating priors (keys) to their url functions (values).
    """
    return {  # Missing priors' urls are generated with a generic function in _get_prior_url()
        "CTYP": _get_ctyp_url,
        "GENG": _get_geng_url,
    }
