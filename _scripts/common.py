"""
Common uses and miscalleneous functions.

@ May 2021
@ Asloudj Yanis
"""

import json
import requests
import os
import shutil
import csv

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