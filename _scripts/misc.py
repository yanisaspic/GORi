"""
Common uses and miscalleneous functions.

@ May 2021
@ Asloudj Yanis
@ Arthur CARON
"""

import json
import requests
import os
import shutil
import csv
import pandas as pd
import mygene as mg

HPO_ANNO_FIELDS = [
    'Gene_ID',
    'Gene_Symbol',
    'HPO_Term_ID',
    'HPO_Term_Label',
    'Raw_Frequency',
    'HPO_Frequency',
    'Additional_Info',
    'Source',
    'Disease_ID'
]

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

def download_url(url, save_path, chunk_size=128):
    """
    # Description
    Downloads a file from an online resource.

    # Arguments
    ``url`` (string): an url leading to a file online. \n
    ``save_path`` (string): a path where the file can be saved.

    # Usage
    >>> download_url(url = "http://current.geneontology.org/ontology/go.obo", 
        save_path = "data/raw/go.obo")
    """
    if not os.path.exists(save_path):
        r = requests.get(url, stream=True)
        with open(save_path, 'wb') as fd:
            for chunk in r.iter_content(chunk_size=chunk_size):
                fd.write(chunk)
        print("WROTE: %s" % save_path)

def replace_first_line(file_path, new_line):
    """
    # Description
    Replaces the first line of a text file.

    # Arguments
    ``file_path`` (string): a path leading to a text file. \n
    ``new_line`` (string): the new contents of the line.

    # Usage
    >>> replace_first_line(file_path = "data/raw/hpo_annotation.txt", new_line = "")
    """
    with open(file_path, 'rt') as old_f:
        with open("%s.tmp" % file_path, 'wt') as new_f:
            old_f.readline()
            new_f.write("%s" % new_line)
            shutil.copyfileobj(old_f, new_f)
    os.replace("%s.tmp" % file_path, file_path)

def download_replace_first(url, save_path, new_line):
    """
    # Description
    Downloads a file and replaces its first line if the file does not exist already.

    # Arguments
    ``url`` (string): an url leading to a file online. \n
    ``save_path`` (string): a path where the file can be saved. \n
    ``new_line`` (string): the new contents of the line.

    # Usage
    >>> download_replace_first(url = "http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt", 
    save_path = "data/raw/hpo_annotation.txt", new_line = "")
    """
    if not os.path.exists(save_path):
        download_url(url, save_path)
        replace_first_line(save_path, new_line)

def get_values_from_keys(keys, dictio):
    """
    # Description
    Returns a list of dictionary values corresponding to a set of keys.

    # Arguments
    ``keys`` (list of strings): keys of the dictionary.
    ``dictio``: a dictionary with values associated to the keys.

    # Usage
    >>> alphabet = {1: 'a', 2: 'b', 3: 'c'}
    >>> print(get_values_from_keys([1, 2], alphabet))
    ... ['a', 'b']
    """
    values = []
    for k in keys:
        try:
            values.append(dictio[k])
        except KeyError:
            pass
    return values

def save_as_json(dictio, filename):
    """
    # Description
    Saves a complex dictionary as a JSON file.

    # Arguments
    ``dictio`` (dict): can include complex and nested Python objects.
    ``filename`` (string): the name of the JSON file.

    # Usage
    >>> small = {'a': 1, 'b':2, 'c': {3, 4}}
    >>> save_as_json(small, "small.json")
    """

    class setEncoder(json.JSONEncoder):
        """Also converts sets to lists and returns the values of a Python object."""
        def default(self, obj):
            if isinstance(obj, set):
                return list(obj)
            try:
                return obj.__dict__
            except AttributeError:
                return json.JSONEncoder.default(self, obj)

    dictio_json_str = json.dumps(dictio, cls=setEncoder)
    dictio_json_file = open(filename, 'wt')
    dictio_json_file.write(dictio_json_str)
    dictio_json_file.close()
    print("WROTE: %s" % filename)

def record_iterator(handle, fields):
    """
    # Description
    Iterate over records in a file.

    # Argument
    ``handle``: a handle corresponding to an open text file. \n
    ``fields`` (list): the corresponding ordered columns names.  

    # Usage
    >>> REAC_PATH_FIELDS = [
    "DB_Object_ID",
    "Path_ID",
    "URL",
    "Path_Label",
    "Evidence",
    "Species"
    ]
    >>> with open("reactome_annotation.txt", 'rt') as ra:
        for path in record_iterator(ra, REAC_ANNO_FIELDS):
            print(path)
    ... {'DB_Object_ID': 'A0A023GPK8', 
    'Path_ID': 'R-DME-373753', 
    'URL': 'https://reactome.org/PathwayBrowser/#/R-DME-373753', 
    'Path_Label': 'Nephrin family interactions', 
    'Evidence': 'IEA', 
    'Species': 'Drosophila melanogaster'}
    """
    for inline in handle:
        inrec = inline.rstrip("\n").split("\t")
        yield dict(zip(fields, inrec))

def read_csv(x):
    """
    # Description :
    Will read a file .csv and put its data in a set

    # Argument :
    ``x`` (str): Name of the file .csv from which we want to extract the data

    # Usage : 
    >>> getGenes = fichier_de_base("GSE31684_Intersect_genes.csv")
    """
    with open(x,'r') as f:
        genes = []
        read = csv.reader(f)

        for i in read : 
            bop = " ".join(i)
            genes.append(bop)
        del genes[0]

    return set(genes)

def split_dict_values(mixed_dict, key_prefix_length = 2):
    """
    # Description
    Splits a dict into subdicts with values splitted according to group of keys.

    # Arguments
    ``mixed_dict`` (dict): the values you want to split according to their key prefixes.
    ``key_prefix_length`` (int): the length of the prefixes used to identify the groups of items.

    # Usage
    >>> melt = {'A:12': 12, 'A:24':32, 'B:10': 0.9, 'B:8': 7.5, 'C:2034': "XYZ"}
    >>> print(split_dict_values(melt, 1))
    ... {'A': [12, 32], 'B': [0.9, 7.5], 'C': ['XYZ']}
    """
    splitted_dict = {}
    for k in mixed_dict.keys():
        k_prefix = k[:key_prefix_length]
        try:
            splitted_dict[k_prefix].append(mixed_dict[k])
        except KeyError:
            splitted_dict[k_prefix] = [mixed_dict[k]]
    return splitted_dict

def truncate(v, n):
    """
    # Description
    Truncates a float and returns it as a string

    # Arguments
    ``v`` (float): the value you want to truncate.
    ``n`` (int): the number of decimals to keep.

    # Usage
    >>> print(truncate(2.1093, 2))
    ... 2.10
    >>> print(truncate(1.0, 4))
    ... 1.0
    """
    s = str(v)
    main, decimals = s.split(".")
    return main + "." + decimals[:n]

def load_hpo_annotation(anno_path):
    """
    # Description
    Returns a dictionary with genes symbols as keys and their respective leaf HPO nodes as a set of values.

    # Arguments
    ``anno_path`` (string): path leading to the HPO annotation file.

    # Usage
    >>> anno = load_hpo_annotation("data/hpo_annotation.txt")
    >>> print(len(anno['IGF1R'])) 
    ... 53
    """
    hpo_anno = {}
    with open(anno_path, 'rt') as ha:
        for rec in record_iterator(ha, HPO_ANNO_FIELDS):
            try:
                hpo_anno[rec['Gene_Symbol']].add(rec['HPO_Term_ID'])
            except KeyError:
                hpo_anno[rec['Gene_Symbol']] = set([rec['HPO_Term_ID']])
    return hpo_anno

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
        for rec in record_iterator(ra, REAC_ANNO_FIELDS):
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
        for rec in record_iterator(rh, REAC_HIER_FIELDS):
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
        for rec in record_iterator(rl, REAC_LABL_FIELDS):
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
    namespace = get_values_from_keys(namespace_ids, label)
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

def translation(liste, input, output, species):
    """
    # Description
    Uses the module MyGene (https://docs.mygene.info/projects/mygene-py/en/latest/) to translate gene information 
    (nomenclature), into a list of dictionary, with different keys such as the 'query', the '_id', the '_score' and its 
    translation (ex. of key : 'uniprot'). For some translation, it is possible to have several results because of the 
    different databases of Uniprot (exemple : from 'symbol' to 'uniprot', in the key 'uniprot', there is another 
    dictionnary with key : 'Swiss-Prot' and 'TrEMBL' in which there is either one or several translation)

    # Arguments
    ``liste`` : list or set of gene that are to be translated to another type of gene representation
    ``input`` : the type of gene representation ('uniprot' / 'symbol' or other, see the link of the module)
    ``output`` : in what the gene is translated ('uniprot' / 'symbol' or other, see the link of the module)

    # Usage 
    >>> print(translation(listOfGene, 'symbol', 'uniprot'))
    """
    out = mg.MyGeneInfo().querymany(liste, scopes = input, fields = output, species = species)
    return out

def writing(liste, retour):
    """
    # Description 
    From the use of translation(), extract the uniprot ID from the Swiss-Prot database keys in dictionnary. 
    Return a dictionnary of queries (keys) and results (item) : {'Symbol':'uniprotID' , ...}

    # Arguments 
    ``liste`` (list): list of dictionary with the queries (contain 'query', '_id', '_score' and translation)
    ``retour`` (str): type of nomenclature to be kept ('symbol' or 'uniprot') in the liste

    # Usage 
    >>> liste = translation(listOfGene, 'symbol', 'uniprot')
    >>> writing(liste, 'uniprot')
    """
    dico = {}

    for i in liste :

        if i.get(retour) is None :
            continue
        
        else :
            if 'Swiss-Prot' in i[retour].keys():
                dico[i['query']] = i[retour]['Swiss-Prot']

            else :
                continue
            
    return dico

def get_symbol_dict(liste, input = 'symbol', output = 'uniprot', species = 'human'):
    """
    # Description
    Calls translation() and writing() to return a dictionary of inputs as keys and outputs as values.
    """
    out = translation(liste, input, output, species)
    return writing(out, output) 