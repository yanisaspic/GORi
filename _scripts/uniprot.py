"""
Goal >  Translate a set of gene symbols to Uniprot IDs with Mygene, and put them in a dictionnary. 

@ June 2021
@ Arthur CARON
"""

import mygene as mg
import time

def get_genes(x):
    with open(x,'r') as f:
        genes = []
        read = csv.reader(f)
        for i in read : 
            bop = " ".join(i)
            genes.append(bop)
        del genes[0]
    return set(genes)

def translation(liste, input, output, species):
    out = mg.MyGeneInfo().querymany(liste, scopes = input, fields = output, species = species)
    return out

def writing(liste, retour):
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
    out = translation(liste, input, output, species)
    return writing(out, output)