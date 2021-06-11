"""
Functions called to create a Reactome Ontology.

@ May 2021
@ Asloudj Yanis
"""

from _scripts import common as cmn

REAC_ANNO_FIELDS = [
    "DB_Object_ID",
    "Path_ID",
    "URL",
    "Path_Label",
    "Evidence",
    "Species"
    ]

REAC_LABL_FIELDS = [
    "Path_ID",
    "Path_Label",
    "Species"
]

REAC_HIER_FIELDS = [
    "Parent_Path_ID",
    "Child_Path_ID"
]


class REACTerm:
    """
    Reactome term created in order to write the .obo file.
    """

    def __init__(self, id, name, namespace, _parents):
        """
        ``id``: unique and stable identifier of the Reactome pathway. \n
        ``name``: human-readbable name of the Reactome pathway. \n
        ``namespace``: name of the root pathway. \n
        ``_parents``: _parents pathways of the Reactome pathway.
        """
        self.id = id
        self.name = name
        self.namespace = namespace
        self._parents = _parents

    def __str__(self):
        """Displays the Reactome ID only."""
        return "REACTerm('%s')\t %s" % (self.id, self.name)
    
    def __repr__(self):
        """Displays all attributes."""
        ret = [
            "REACTerm('%s'):\n" % self.id,
            "name:\t%s\n" % self.name,
            "_parents:\t{%s items}\t%s" % (len(self._parents), self._parents)]
        return "".join(ret)

def load_reactome_annotation(anno_path):
    """
    # Description
    Returns a dictionary with genes ids as keys and their respective leaf Reactome nodes as a set of values.

    # Arguments
    ``anno_path`` (string): path leading to the Reactome annotation file.

    # Usage
    >>> anno = load_reactome_annotation("data/raw/reactome_annotation.txt")
    >>> print(anno['P08069']) 
    ... {'R-HSA-2404192', 'R-HSA-9009391', 'R-HSA-2428928', 'R-HSA-2428933'}
    """
    reactome_anno = {}
    with open(anno_path, 'rt') as ra:
        for rec in cmn.record_iterator(ra, REAC_ANNO_FIELDS):
            try:
                reactome_anno[rec['DB_Object_ID']].add(rec['Path_ID'])
            except KeyError:
                reactome_anno[rec['DB_Object_ID']] = set([rec['Path_ID']])
    return reactome_anno

def load_reactome_hierarchy(hier_path):
    """
    # Description
    Returns a dict with REACTerms as keys and their direct _parents as a set of values.

    # Arguments
    ``hier_path`` (string): path leading to the Reactome hierarchy file.

    # Usage
    >>> hier = load_reactome_hierarchy("data/raw/reactome_hierarchy.txt")
    >>> print(hier['R-HSA-198753'])
    ... {'R-HSA-198725', 'R-HSA-450282'}
    """
    reactome_hierarchy = {}
    with open(hier_path, 'rt') as rh:
        for rec in cmn.record_iterator(rh, REAC_HIER_FIELDS):
            try:
                reactome_hierarchy[rec['Child_Path_ID']].add(rec['Parent_Path_ID'])
            except KeyError:
                reactome_hierarchy[rec['Child_Path_ID']] = set([rec['Parent_Path_ID']])
    return reactome_hierarchy

def load_reactome_label(labl_path):
    """
    # Description
    Returns a dict with REACTerms ID as keys and their human-readable label as values.

    # Arguments
    ``labl_path`` (string): path leading to the Reactome label file.

    # Usage
    >>> labl = load_reactome_label("data/raw/reactome_label.txt")
    >>> print(labl['R-HSA-198753'])
    ... ERK/MAPK targets
    """
    reactome_labl = {}
    with open(labl_path, 'rt') as rl:
        for rec in cmn.record_iterator(rl, REAC_LABL_FIELDS):
            reactome_labl[rec['Path_ID']] = rec['Path_Label']
    return reactome_labl

def get_pathway_roots(path_id, hierarchy):
    """
    # Description
    Returns the IDs of the root nodes of a Reactome pathway.

    # Arguments
    ``path_id`` (string): unique and stable identifier of a Reactome pathway. \n
    ``hierarchy`` (dict): Reactome pathways and their respective parents.

    # Usage
    >>> hier = load_reactome_hierarchy("data/raw/reactome_hierarchy.txt")
    >>> print(get_reactome_namespace('R-HSA-198753', hier))
    ... {'R-HSA-168256', 'R-HSA-162582'}
    """
    namespace = set()
    try:
        for parent_id in hierarchy[path_id]:
            namespace.update(get_pathway_roots(parent_id, hierarchy))
    except KeyError:
        return [path_id]
    return namespace

def get_path_namespace(path_id, hierarchy, label):
    """
    # Description
    Returns the namespaces (i.e. names of the root nodes) of a Reactome pathway in a single string.
    
    # Arguments
    ``path_id`` (string): unique and stable identifier of a Reactome pathway. \n
    ``hierarchy`` (dict): Reactome pathways and their respective parents. \n
    ``label`` (dict): Reactome pathways and their respective human-readable names.

    # Usage
    >>> hier = load_reactome_hierarchy("data/raw/reactome_hierarchy.txt")
    >>> labl = load_reactome_label("data/raw/reactome_label.txt")
    >>> print(get_path_namespace('R-HSA-198753', hier, labl))
    ... Immune System & Signal Transduction
    """
    namespace_ids = get_pathway_roots(path_id, hierarchy)
    namespace = cmn.get_values_from_keys(namespace_ids, label)
    return " & ".join(namespace)

def get_reacterm_dict(hierarchy, label):
    """
    # Description
    Returns a dictionary of REACTerms.

    # Arguments
    ``hierarchy`` (dict): Reactome pathways and their respective parents. \n
    ``label`` (dict): Reactome pathways and their respective human-readable names.

    # Usage
    >>> hier = load_reactome_hierarchy("data/raw/reactome_hierarchy.txt")
    >>> labl = load_reactome_label("data/raw/reactome_label.txt")
    >>> pseudo_reactome_onto = get_reactome_dict(hier, labl)
    >>> print(pseudo_reactome_onto['R-HSA-198753'])
    ... REACTerm('R-HSA-198753')         ERK/MAPK targets
    """
    pseudo_reactome_onto = {}
    for k in label.keys():
        try:
            parents = hierarchy[k]
        except KeyError:
            parents = set()
        pseudo_reactome_onto[k] = REACTerm(k, label[k], get_path_namespace(k, hierarchy, label), parents)
    return pseudo_reactome_onto

def load_reacterm_dict(hier_path, labl_path):
    """
    # Description
    Returns a dictionary of REACTerms after loading the necessary files.

    # Arguments
    ``hier_path`` (string): path leading to the Reactome hierarchy file.
    ``labl_path`` (string): path leading to the Reactome label file.

    # Usage
    >>> reacterm = load_reacterm_dict("data/raw/reactome_hierarchy.txt", "data/raw/reactome_label.txt")
    >>> print(reacterm['R-HSA-198753'])
    ... REACTerm('R-HSA-198753')         ERK/MAPK targets
    """
    hier = load_reactome_hierarchy(hier_path)
    labl = load_reactome_label(labl_path)
    return get_reacterm_dict(hier, labl)