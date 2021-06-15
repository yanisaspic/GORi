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
    >>> oboDag = goatools.obo_parser.GODag("go-basic.obo")
    >>> oboDag.pop(x)
    >>> save_as_obo(oboDag, "smaller_go-basic.obo")
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
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters as keys.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go-basic.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> print(len(get_term_ontology("GO:0005524", onto)))
    ... 47284
    """
    return ontologies[(term_id[:2])]

def term_weight(term_id, ontologies):
    """
    # Description
    Returns a weight between 0 and 1 corresponding to the Information Content (IC) of an ontology term.
    IC = p/(p+c) with *p = number of parent terms and *c = number of children terms.

    # Arguments
    ``term_id`` (string): the unique and stable identifier of an ontology term. \n
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters as keys.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go-basic.obo", optional_attrs = "relationship"),
        'R-': obo_parser.GODag("reactome.obo", optional_attrs = "relationship")}
    >>> print(term_weight("GO:0005524", onto))
    ... 1
    """
    source = get_term_ontology(term_id, ontologies)
    term = source[term_id]
    top_lineage = len(term.get_all_parents())
    bottom_lineage = len(term.get_all_children())
    return (top_lineage + 1) / (top_lineage + bottom_lineage + 1)