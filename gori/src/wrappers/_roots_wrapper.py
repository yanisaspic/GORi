"""Functions used by roots_wrapper().
    Roots are manually defined for GO and CTYP due to their specific data structures.

    2025/05/26 @yanisaspic"""


def _get_biop_roots() -> set[str]:
    """Get the roots of the BIOP prior.

    Returns
        A set of GO ids.
    """
    return {"GO:0008150"}


def _get_celc_roots() -> set[str]:
    """Get the roots of the CELC prior.

    Returns
        A set of GO ids.
    """
    return {"GO:0005575"}


def _get_ctyp_roots() -> set[str]:
    """Get the roots of the CTYP prior.

    Returns
        A set of CL ids.
    """
    return {"CL:0000000"}


def _get_molf_roots() -> set[str]:
    """Get the roots of the MOLC prior.

    Returns
        A set of GO ids.
    """
    return {"GO:0003674"}
