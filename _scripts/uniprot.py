"""
Goal >  Translate a set of gene symbols to Uniprot IDs with Mygene, and put them in a dictionnary. 

@ June 2021
@ Arthur CARON
"""

import mygene as mg
import time as time

#exemple ID Uniprot to use MyGene
xli = ['A0A023GPK8','P08069','HJHKKK','P08069']
#exemple symbols to use MyGene
xlu = ['fred','IGF1R','IGF1R','EEF1AKNMT']


begin = time.time()


def translation(liste, input, output, species = 'human'):

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

end = time.time()

print(writing(translation("""set of genes""",'symbol','uniprot'),'uniprot'))

print(f"total time = {end - begin}s")