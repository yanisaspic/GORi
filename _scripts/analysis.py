import numpy as np
import matplotlib.artist as artists
import matplotlib.pyplot as plt
import statistics as stc
from scipy import stats

def four_histograms(data, listeTitles = ["GO","Reactome","HPO","MIX"], listeOfColors = ['blue','blue','blue','blue']):
    """
    # Description 
    Return four histograms from a list of data, with its Mean, its Median and its Standard Deviation

    # Arguments
    ``data`` (list, set, array...): data used to create the histogram 
    ``listeTitles`` (list of str): Titles of the graphs
    ``listeOfColors`` (list of str) : list with the colors of the graphs if statistical analysis has
    been done on the data before hand. Graphs with the same colors are significantly similars, by a 
    Kolmogorov-Smirnov test.

    # Usage 
    >>> four_histograms(liste,listeTitres)
    >>> plt.tight_layout()
    >>> plt.show

    """
    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    compt = 0

    for i in range(1, 5):
        ax = fig.add_subplot(2, 2, i) # i est la position dans le subplot
        median = round(stc.median(data[compt]),2)
        mean = round(stc.mean(data[compt]),2)
        ecart = round(stc.pstdev(data[compt]),2)
        text = f'Mean = {mean}\nMedian = {median}\nStandard Deviation = {ecart}'
        props = dict(boxstyle='square', alpha=0.5   , color = 'grey')

        ax.title.set_text(f"{listeTitles[compt]}")
        ax.grid(axis='y', alpha=0.75)
        ax.set_xlabel('Value')
        ax.set_ylabel('Frequency') 
        ax.hist(x = data[compt], bins = 10, color=listeOfColors[compt], alpha=0.5, rwidth=0.50)
        ax.text(0.05,0.95, text, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

        compt += 1

    plt.tight_layout()
    plt.show