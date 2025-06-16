"""Functions used by roots_wrapper().
    Roots are manually defined for GO and Cell Ontology due to their specific data structures.

    2025/05/26 @yanisaspic"""


def _get_biop_roots() -> set[str]:
    """Get the roots of the GO_BP resource.

    Returns
        A set of GO ids.
    """
    return {"GO:0008150"}


def _get_celc_roots() -> set[str]:
    """Get the roots of the GO_CC resource.

    Returns
        A set of GO ids.
    """
    return {"GO:0005575"}


def _get_ctyp_roots() -> set[str]:
    """Get the roots of the Cell Ontology resource.

    Returns
        A set of CL ids.
    """
    return {"CL:0000000"}


def _get_molf_roots() -> set[str]:
    """Get the roots of the MOLC resource.

    Returns
        A set of GO ids.
    """
    return {"GO:0003674"}
