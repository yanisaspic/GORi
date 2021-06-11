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
    Saves a dictionary of terms as a light v1.2 .obo file. \n 
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

def get_ancestry_id(gene_id, annotation, ontology, relationship = False):
    """
    # Description
    Returns the set of ids corresponding to a gene's leaf, branch and root nodes using an ontology. \n
    Nodes which aren't parents (i.e. 'is a') but share a relationship (e.g. 'regulates') can be excluded.

    # Arguments
    ``gene_id`` (set of strings): the unique id corresponding to your gene on UniProtKB. \n
    ``annotation`` (dict): genes and their corresponding leaf nodes. \n
    ``ontology```(GODag object): ontology corresponding to the terms used. \n
    ``relationship`` (boolean): should non-parents be included ?

    # Usage
    >>> ancestry_with_cousins = get_ancestry_goid(gene_id = "P08069", gene_go_annotation, go_onto, relationship = True)
    >>> print(len(ancestry_with_cousins))
    ... 321
    >>> print(ancestry_with_cousins))
    ... {'GO:0042326', 'GO:0051896', ...}
    >>> ancestry_without_cousins = get_ancestry_goid(gene_id = "P08069")
    >>> print(len(ancestry_without_cousins))
    ... 270
    """
    leaf_nodes = annotation[gene_id]
    tree = set()
    for t_id in leaf_nodes:
        tree.add(t_id)
        term = ontology.query_term(t_id)
        if relationship:
            tree.update(term.get_all_upper())
        else:
            tree.update(term.get_all_parents())
    return tree

def get_term_ontology(term_id, ontologies):
    """
    # Description
    Returns the ontology corresponding to an ID using its first 2 characters.

    # Arguments
    ``term_id`` (string): the unique and stable identifier of an ontology term.
    ``ontologies`` (dict of GODag objects): ontologies as values and their terms' first 2 characters as keys.

    # Usage
    >>> onto = {
        'GO': obo_parser.GODag("go-basic.obo"),
        'R-': obo_parser.GODag("reactome.obo")}
    >>> print(len(get_term_ontology("GO:0005524", onto)))
    ... 47284
    """
    return ontologies[term_id[:2]]

def term_weight(termID, ontology):
    source = get_term_ontology(termID, ontology)
    term = source[termID]

    param = {}
    param['name'] = term.name
    param['level'] = term.level
    param['depth'] = term.depth
    param['paths_to_top'] = len(source.paths_to_top(termID))

    weight = math.exp( param['depth'] - param['level']) * param['level'] * param['depth'] * param['paths_to_top']

    return weight