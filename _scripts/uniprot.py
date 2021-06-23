"""
Goal >  Translate a set of gene symbols to Uniprot IDs with Mygene, and put them in a dictionnary. 

@ June 2021
@ Arthur CARON
"""

import mygene as mg
import time as time

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