"""
Functions related to ontology managing.

@ May 2021
@ Asloudj Yanis
"""

import math
from goatools import obo_parser

def save_as_obo(dictio, filename, header):
    """
    # Description
    Saves a dictionary of terms as a lighter v1.2 .obo file. \n 
    A term is a Python object with id, name, namespace, and _parents attributes. \n
    A term can also have a relationship attribute (e.g. part_of).

    # Arguments
    ``dictio`` (dict): GOTerms or REACTerms.
    ``filename`` (string): the name of the .obo file.
    ``header`` (string): header added to the .obo file.

    # Usage
    >>> oboDag = goatools.obo_parser.GODag("go.obo")
    >>> oboDag.pop(x)
    >>> save_as_obo(oboDag, "smaller_go.obo")
    """

    # .obo file header:
    dictio_obo_file = open(filename, 'wt')
    dictio_obo_file.write("format-version: 1.2\n")
    dictio_obo_file.write("%s\n\n" % header)

    for k in dictio.keys():
        term = dictio[k]

        parents = set()
        for par in term._parents:
            parents.add("is_a: %s ! %s\n" % (par, dictio[par].name))

        try:
            relations = set()
            for rel in term.relationship.keys():
                # e.g. rel = 'part_of'
                for cousin in term.relationship[rel]:
                    relations.add("relationship: %s %s ! %s\n" % (rel, cousin.id, cousin.name))
        except AttributeError:
            relations = ""

        # adds the .obo term:
        entry = "".join(
            ["[Term]\n",
            "id: %s\n" % term.id,
            "name: %s\n" % term.name,
            "namespace: %s\n" % term.namespace,
            "%s" % "".join(parents),
            "%s\n" % "".join(relations)]
        )
        dictio_obo_file.write(entry)
    print("WROTE: %s" % filename)

def get_term_ontology(term_id, ontologies):
    """
    # Description
    Returns the ontology corresponding to an ID using its first 2 characters.

    # Arguments
    ``term_id`` (string): the unique and stable identifier of an ontology term. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters 
    as keys.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> print(len(get_term_ontology("GO:0005524", onto)))
    ... 47284
    """
    return ontologies[(term_id[:2])]

def get_term_ends(term_id, ontologies, roots = False):
    """
    # Description
    Returns the ID of all the leaf or root nodes for a given ontology term.

    # Arguments 
    ``term_id`` (string): the unique and stable identifier of an ontology term. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters 
    as keys.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> print(get_term_ends('GO:0031424', onto))
    ... {'GO:1905716', 'GO:1905717'}
    >>> print(get_term_ends('R-HSA-975634', onto, roots=True))
    ... {'R-HSA-9709957', 'R-HSA-1430728'}
    """
    extremes = set()
    source = get_term_ontology(term_id, ontologies)
    term = source[term_id]

    # get all the parents or all the children of a term (according to roots value)
    if roots:
        tree_terms_id = term.get_all_upper()
    else:
        tree_terms_id = term.get_all_lower()

    # root = a term without upper terms
    # leaf = a term without lower terms
    for tree_tid in tree_terms_id:
        tree_term = source[tree_tid]
        if roots:
            if len(tree_term.parents) + len(tree_term.relationship) == 0:
                extremes.add(tree_tid)
        else:
            if len(tree_term.children) + len(tree_term.relationship_rev) == 0:
                extremes.add(tree_tid)
    return extremes

def get_one_intrinsic_IC(term_id, ontologies, max_leaves = {}):
    """
    # Description
    Returns the Information Content of an ontology term t according to its features in the ontology.
    numerator = ( leaves(t)/subsumers(t) ) + 1
    denominator = leaves(root) + 1
    IC(t) = -log(numerator/denominator) \n
    @ Ontology-based information content computation, Sánchez et al. (2011)
    @ Knowledge Based Systems 24 (297-303)

    # Arguments
    ``term_id`` (string): the unique and stable identifier of an ontology term. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters 
    as keys. \n
    ``max_leaves`` (dict): root term ids as keys and their number of leaves as values.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> print(get_one_intrinsic_IC('GO:0009913', onto))
    ... 8.3, {'GO:0008150': 12086}
    >>> print(get_one_intrinsic_IC('GO:1905716', onto))
    ... 9.4, {'GO:0008150': 12086}
    """
    # numerator
    n_term_leaves = len(get_term_ends(term_id, ontologies))
    source = get_term_ontology(term_id, ontologies)
    term = source[term_id]
        # include the term itself with +1 (e.g. if term is root, subsumers(term) = 1)
    n_term_subsumers = len(term.get_all_upper()) + 1
    numerator = n_term_leaves / n_term_subsumers + 1

    # denominator
    term_roots = get_term_ends(term_id, ontologies, roots=True)
    # if the term is a root
    if len(term_roots) < 1:
        term_roots.add(term_id)
    n_root_leaves = 0
    for root in term_roots:
        try:
            n_root_leaves = max_leaves[root]
        except KeyError:
            max_leaves[root] = len(get_term_ends(root, ontologies))
            n_root_leaves += max_leaves[root]
    denominator = n_root_leaves + 1

    # intrinsic IC
    IC = round(-math.log(numerator/denominator), 2)
    return IC, max_leaves

def get_all_intrinsic_IC(terms_ids, ontologies):
    """
    # Description
    Returns the relative IC content of ontology terms according to their features in the
    ontology. The IC are returned in a dict with terms ids as keys and values in [0, 1].

    # Arguments
    ``terms_ids`` (list of strings): the unique and stable identifier of an ontology term. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters 
    as keys.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> terms = [
        "GO:0031424", "GO:0009888", "GO:0008150", 
        "R-HSA-9709957", "R-HSA-975634", "R-HSA-6806667"]
    >>> print(get_all_intrinsic_IC(terms, onto))
    ... {'GO:0031424': 1.0, 'GO:0009888': 0.44, 'GO:0008150': -0.0, 'R-HSA-9709957': -0.0, 
    'R-HSA-975634': 0.46, 'R-HSA-6806667': 1.0}
    """
    all_max_leaves = {}
    max_int_IC = {}
    int_IC = {}
    for term_prefix in ontologies.keys():
        max_int_IC[term_prefix] = 0
    
    # calculate the IC values
    for term in terms_ids:
        term_prefix = term[:2]
        term_IC, max_leaves = get_one_intrinsic_IC(term, ontologies, all_max_leaves)
        all_max_leaves.update(max_leaves)
        int_IC[term] = term_IC

        if term_IC > max_int_IC[term_prefix]:
            max_int_IC[term_prefix] = term_IC

    # divide them by the maximal value of the ontology
    for term in int_IC.keys():
        term_prefix = term[:2]
        int_IC[term] = round(int_IC[term] / max_int_IC[term_prefix], 2)
    return int_IC