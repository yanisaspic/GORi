"""
Gene HPO annotations parser and loader.

@ May 2021
@ Asloudj Yanis
"""

from _scripts import common as cmn

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

def load_hpo_annotation(anno_path):
    """
    # Description
    Returns a dictionary with genes symbols as keys and their respective leaf HPO nodes as a set of values.

    # Arguments
    ``anno_path`` (string): path leading to the HPO annotation file.

    # Usage
    >>> anno = load_hpo_annotation("data/raw/hpo_annotation.txt")
    >>> print(len(anno['IGF1R'])) 
    ... 53
    """
    hpo_anno = {}
    with open(anno_path, 'rt') as ha:
        for rec in cmn.record_iterator(ha, HPO_ANNO_FIELDS):
            try:
                hpo_anno[rec['Gene_Symbol']].add(rec['HPO_Term_ID'])
            except KeyError:
                hpo_anno[rec['Gene_Symbol']] = set([rec['HPO_Term_ID']])
    return hpo_anno