"""
> run a Kolmogorov-Smirnov test on two-sample over a list of data-set and return a colors if 
they are significantly similars.

@ Arthur CARON
@ June 2021
"""

from scipy import stats
import numpy as np
import copy 

# data to use as an example 
test = [1,2,3,4]
test2 = [1,2,3,4]
test3 = [5,6,7,8]
test4 = [5,6,7,8]
listedetest = [test,test3,test4,test2]

def comparison_to_colors(listeTest):

    """
    # Description 
    Return a list of colors, up to 4 colors from a list of data-set. The aim is to realise a Kolmogorov
    Smirnov test on each data-set, two-by-two, and affect the same colors to data-sets that are 
    significantly similars, and different colors to data-sets that are not. It returns a list of colors
    in the same order as the data is ordered in the list.

    # Arguments
    ``listeTest`` (list): list of data-set, to be compared and affected colors. 

    # Usage 
    >>> listedetest = [test,test3,test4,test2]
    >>> print(comparison_to_colors(listedetest))
    >>> ['blue', 'green', 'green', 'blue']

    """

    forColor = copy.deepcopy(listeTest)
    listOfColors = []

    dicoCorresp = {}
    dicoCorresp['blue'] = []
    dicoCorresp['green'] = []
    dicoCorresp['red'] = []
    dicoCorresp['yellow'] = []

    n = 1
    toCopy = listeTest

    for key in dicoCorresp.keys() :

        dicoCorresp[key] = toCopy
        toCopy = []

        while n < len(dicoCorresp[key]) :

            test = stats.ks_2samp(dicoCorresp[key][0],dicoCorresp[key][n])

            if test.pvalue <= 0.05 :
                toCopy.append(dicoCorresp[key][n])
                del dicoCorresp[key][n]
                
            else:
                n += 1 

    for j in forColor : 
        for key in dicoCorresp.keys():
            if j in dicoCorresp[key] :
                listOfColors.append(key)

    return listOfColors
