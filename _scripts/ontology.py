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

def get_extrinsic_IC(annotations):
    """
    # Description
    Returns a dict of terms as keys and their Information Content (IC) as values. The IC is calculated
    according to the frequencies of each term in the annotations.

    # Arguments
    ``annotations`` (dict of dict): a dict of genes as keys and their terms as values. The terms must 
    be grouped according to their source ontology.

    # Usage
    >>> annotations = {
        'GeneA': {
            'GO': ["GO:0009987", "GO:0008150"], 
            'R-': ["R-HSA-73857", "R-HSA-74160"],
            'HP': set()
            },
        'GeneB': {
            'GO': ["GO:0008150"], 
            'R-': ["R-HSA-74160"],
            'HP': ["HP:0000001"]
            },
        'GeneC': {
            'GO': ["GO:0009987", "GO:0008150", "GO:0030029"], 
            'R-': ["R-HSA-73857", "R-HSA-74160", "R-HSA-1643685"],
            'HP': ["HP:0000118", "HP:0000001"]
            }
        }
    >>> print(get_extrinsic_IC(annotations))
    ... {'GO:0009987': 0.37, 'GO:0008150': -0.0, 'R-HSA-73857': 0.37, 'R-HSA-74160': -0.0, 
    'HP:0000001': -0.0, 'GO:0030029': 1.0, 'R-HSA-1643685': 1.0, 'HP:0000118': 0.63}
    """
    ## Reads the annotations to get the absolute frequencies of the terms and the number
    ## genes annotated by each ontology
    annotated_genes_count = {}
    terms_count = {}
    for gene in annotations.keys():
        one_gene_annotations_from_all_onto = annotations[gene]
        for ontology in one_gene_annotations_from_all_onto.keys():
            one_gene_annotations_from_one_onto = one_gene_annotations_from_all_onto[ontology]
            if len(one_gene_annotations_from_one_onto) > 0:
                annotated_genes_count[ontology] = annotated_genes_count.get(ontology, 0) + 1
                for term in one_gene_annotations_from_one_onto:
                    terms_count[term] = terms_count.get(term, 0) + 1

    ## Reads the absolute frequencies and the ontology counters to get the relative frequencies
    min_freq = 1
    for term in terms_count.keys():
        term_prefix = term[:2]
        term_freq = terms_count[term] / annotated_genes_count[term_prefix]
        terms_count[term] = term_freq
        if term_freq <= min_freq:
            min_freq = term_freq

    ## Reads the relative frequencies to get the relative Information Content of a term
    ## IC = -log(f)
    max_IC = -math.log(min_freq)
    for term in terms_count.keys():
        terms_count[term] = round(-math.log(terms_count[term]) / max_IC, 2)
    return terms_count