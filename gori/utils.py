"""Miscalleneous functions called to conduct a GORi analysis.

    2025/05/15 @yanisaspic"""

import networkx as nx
from typing import Any, Optional
from pypath.utils.go import GOAnnotation
from pypath.utils.mapping import id_from_label, id_from_label0


def _get_uniprot_id(gene: str) -> str:
    """Get a UniProtID.

    ``gene`` is a gene symbol (or a UniProtID).

    Returns
        A UniProtID.
    """
    uids = id_from_label(gene)  # use id_from_label to get a stable output.
    tmp = sorted(list(uids))
    uid = tmp[0]
    if len(uids) > 1:
        print(
            f"Multiple UniProtIDs found for gene {gene}: {uids}. Selecting {uid} (alphanumeric order)."
        )
    return uid


def _get_generic_ancestors(
    terms: set[str], hierarchy: dict[str, set[str]], ancestors: set[str] = set()
) -> set[str]:
    """Get the ancestors of multiple terms.
    An ancestor is defined as an element ascending any of the input terms.

    ``terms`` is a set of annotation terms.
    ``hierarchy`` is a dict associating terms (keys) to their parents (values).
    ``ancestors`` is a set of terms extended recursively.

    Returns
        A set of ancestor terms.
    """
    tmp = [hierarchy.get(e, set()) for e in terms]
    parents = {parent for subset in tmp for parent in subset}
    terms = parents.difference(ancestors)
    ancestors = ancestors.union(parents)
    if len(terms) > 0:
        ancestors = _get_generic_ancestors(terms, hierarchy, ancestors)
    return ancestors


def _get_generic_descendants(
    terms: set[str], hierarchy: dict[str, set[str]]
) -> set[str]:
    """Get the descendants of multiple terms.
    A descendant is defined as an element descending any of the input terms.

    ``terms`` is a set of annotation terms.
    ``hierarchy`` is a dict associating terms (keys) to their parents (values).
    ``ancestors`` is a set of terms extended recursively.

    Returns
        A set of descendant terms.
    """
    ancestry = {e: _get_generic_ancestors({e}, hierarchy) for e in hierarchy.keys()}
    is_descendant = lambda e: len(ancestry.get(e, set()).intersection(terms)) > 0
    descendants = {d for d, ancestors in ancestry.items() if is_descendant(d)}
    return descendants


def _prune_hierarchy(
    hierarchy: dict[str, set[str]],
    roots: Optional[set[str]] = None,
    leaves: Optional[set[str]] = None,
) -> dict[str, set[str]]:
    """Prune a hierarchy, i.e. get a sub-hierarchy with specific roots and leaves.
    If both roots and leaves are indicated, terms associated to both of them are returned.

    ``hierarchy`` is a dict associating terms (keys) to their parents (values).
    ``roots`` is an optional set of terms to use as roots for the new hierarchy.
    ``leaves`` is an optional set of terms to use as leaves for the new hierarchy.

    Returns
        A dict associating terms (keys) to their parents (values).
    """
    ancestors = (
        set(hierarchy.keys())
        if leaves is None
        else _get_generic_ancestors(leaves, hierarchy) | leaves
    )
    descendants = (
        set(hierarchy.keys())
        if roots is None
        else _get_generic_descendants(roots, hierarchy) | roots
    )
    terms = ancestors.intersection(descendants)
    pruned_hierarchy = {e: hierarchy.get(e, set()) for e in terms}
    return pruned_hierarchy


def _get_ctyp_ancestors(terms: set[str], prior: dict[str, dict[str, Any]]) -> set[str]:
    """Get every CTYP ancestor terms.
    A term is ancestor if at least one of the input CTYP ontology term descends from it.

    ``terms`` is a set of CTYP ontology terms.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CTYP ids.
    """
    tmp = [nx.ancestors(prior["ontology"], t) for t in terms]
    ancestors = {ancestor for ancestors0 in tmp for ancestor in ancestors0}
    return ancestors


def _get_ctyp_descendants(
    terms: set[str], prior: dict[str, dict[str, Any]]
) -> set[str]:
    """Get every CTYP descendant terms.
    A term is descendant if at least one of the input CTYP ontology term is its descendant.

    ``terms`` is a set of CTYP ontology terms.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A set of CTYP ids.
    """
    tmp = [nx.descendants(prior["ontology"], t) for t in terms]
    descendants = {descendant for descendants0 in tmp for descendant in descendants0}
    return descendants


def _get_go_ancestors(terms: set[str], prior: GOAnnotation) -> set[str]:
    """Get every GO ancestor terms.
    A term is ancestor if at least one of the input GO ontology term descends from it.

    ``terms`` is a set of GO ontology terms.
    ``prior`` is a GOAnnotation object.

    Returns
        A set of GO ids.
    """
    _get_ancestors0 = lambda t: prior.ontology.get_all_ancestors(
        t, relations="is_a"
    ).difference({t})
    tmp = [_get_ancestors0(t) for t in terms]
    ancestors = {ancestor for ancestors0 in tmp for ancestor in ancestors0}
    return ancestors


def _get_go_descendants(terms: set[str], prior: GOAnnotation) -> set[str]:
    """Get every GO descendant terms.
    A term is descendant if at least one of the input GO ontology term is its ancestor.

    ``terms`` is a set of GO ontology terms.
    ``prior`` is a GOAnnotation object.

    Returns
        A set of GO ids.
    """
    _get_descendants0 = lambda t: prior.ontology.get_all_descendants(
        t, relations="is_a"
    ).difference({t})
    tmp = [_get_descendants0(t) for t in terms]
    descendants = {descendant for descendants0 in tmp for descendant in descendants0}
    return descendants


def _get_prior_ancestors(terms: set[str], prior: str, data: dict[str, Any]) -> set[str]:
    """A wrapper to get the ancestor annotations of a term from a prior.

    ``terms`` is a set of annotation terms.
    ``prior`` is the name of a knowledge base to use.
    ``data`` is a dict associating priors (keys) to their contents (values).

    Returns
        A set of ancestor terms.
    """
    if prior in ["BIOP", "CELC", "MOLF"]:
        ancestors = _get_go_ancestors(terms, data[prior])
    elif prior == "CTYP":
        ancestors = _get_ctyp_ancestors(terms, data[prior])
    else:
        ancestors = _get_generic_ancestors(terms, data[prior]["hierarchy"])
    return ancestors


def _get_prior_descendants(
    terms: set[str], prior: str, data: dict[str, Any]
) -> set[str]:
    """A wrapper to get the descendant annotations of a term from a prior.

    ``terms`` is a set of annotation terms.
    ``prior`` is the name of a knowledge base to use.
    ``data`` is a dict associating priors (keys) to their contents (values).

    Returns
        A set of descendant terms.
    """
    if prior in ["BIOP", "CELC", "MOLF"]:
        descendants = _get_go_descendants(terms, data[prior])
    elif prior == "CTYP":
        descendants = _get_ctyp_descendants(terms, data[prior])
    else:
        descendants = _get_generic_descendants(terms, data[prior]["hierarchy"])
    return descendants


def _get_generic_translation(term: str, prior: dict[str, dict[str, Any]]) -> str:
    """Get the human-readable label of an annotation.

    ``term`` is an annotation.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A human-readable label.
    """
    return prior["translations"][term]


def _get_ctyp_translation(term: str, prior: dict[str, Any]) -> str:
    """Get the human-readable label of a CTYP annotation.

    ``term`` is a CTYP annotation.
    ``prior`` is a dict with two keys: "annotations" and "ontology".

    Returns
        A human-readable label.
    """
    return prior["ontology"].nodes[term]["name"]


def _get_go_translation(term: str, prior: GOAnnotation) -> str:
    """Get the human-readable label of a GO annotation.

    ``term`` is a GO annotation.
    ``prior`` is a GOAnnotation object.

    Returns
        A human-readable label.
    """
    return prior.get_name(term)


def _get_geng_translation(term: str, prior: dict[str, dict[str, Any]]) -> str:
    """Get the human-readable label of a GENG annotation.

    ``term`` is a GENG annotation.
    ``prior`` is a dict with three keys: "annotations", "hierarchy" and "translations".

    Returns
        A human-readable label.
    """
    return term


def _get_prior_translation(term: str, prior: str, data: dict[str, Any]) -> str:
    """Get the human-readable label of an annotation.

    ``terms`` is an annotation term.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).

    Returns
        A human-readable label.
    """
    tmp = term.split(":", 1)[1]
    if prior in ["BIOP", "CELC", "MOLF"]:
        return _get_go_translation(tmp, data[prior])
    if prior == "CTYP":
        return _get_ctyp_translation(tmp, data[prior])
    if prior == "GENG":
        return _get_geng_translation(tmp, data[prior])
    return _get_generic_translation(tmp, data[prior])


def _get_url_headers() -> dict[str, str]:
    """Get url headers specific to priors.

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


def _get_prior_url(term: str, prior: str, data: dict[str, Any]) -> Optional[str]:
    """Get an url leading to an annotation definition.

    ``terms`` is an annotation term.
    ``prior`` is a prior label.
    ``data`` is a dict associating priors (keys) to their contents (values).

    Returns
        A url.
    """
    headers = _get_url_headers()
    if prior in headers.keys():
        tmp = term.split(":", 1)[1]
        if prior == "CTYP":
            tmp = tmp.replace(":", "_")
        if prior == "GENG":
            translations = {
                label: id for id, label in data[prior]["translations"].items()
            }
            tmp = translations[tmp]
        url = headers[prior] + tmp
        return url
    return None
