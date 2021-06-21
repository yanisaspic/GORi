"""
Goal > From datas of GPO, HPO and Reactome, create 1 histogram from each and one more with the three combined 

@ June 2021
@ Arthur CARON
"""

import numpy as np
import matplotlib.artist as artists
import matplotlib.pyplot as plt
import statistics as stc

# Tests as an example :
# test = np.random.sample((1000,))
# test2= np.random.normal(0,0.1,1000)
# test3= np.random.normal(0,1,1000)
# test4= np.random.normal(0,0.5,1000)

def histogram(ax, list, titre):
    """
    # Description 
    Return a Histogram from a list of data, with its Mean, its Median and its Standard Deviation

    # Arguments
    ``ax`` (axes object):   plt.subplot return a tuple containing a figure for an image file (fig) and axes objects (ax). 
                            We use only the axe object to create the graph, to modify the X and Y axis, the title ...
    ``list`` (list, set, array...): data used to create the histogram 
    ``titre`` (str): Title of the graph

    # Usage 
    >>> test1 = np.random.sample((1000,))
    >>> fig, ax = plt.subplot()
    >>> histogram(ax, test1, 'Uniform distribution')
    >>> plt.tight_layout()
    >>> plt.show

    or with several graph : 
    >>> test1 = np.random.sample((1000,))
    >>> test2 = np.random.normal(0,0.1,1000)
    >>> fig, (ax1,ax2) = plt.subplot(nrows = 1, ncols = 2)
    >>> histogram(ax1, test1, 'Uniform distribution')
    >>> histogram(ax2, test2, 'Normal distribution')
    >>> plt.tight_layout()
    >>> plt.show

    """
    median = round(stc.median(list),2)
    mean = round(stc.mean(list),2)
    ecart = round(stc.pstdev(list),2)
    text = f'Mean = {mean}\nMedian = {median}\nStandard Deviation = {ecart}'
    props = dict(boxstyle='square', alpha=0.5   , color = 'grey')

    ax.title.set_text(f'{titre}')
    ax.grid(axis='y', alpha=0.75)
    ax.set_xlabel('Value')
    ax.set_ylabel('Frequency') 
    ax.hist(x = list, bins = 10, color='#0504aa', alpha=0.5, rwidth=0.50)
    ax.text(0.05,0.95, text, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

#example of an application
# fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2)
# histogram(ax1, test, 'oui')
# histogram(ax2, test2, 'non')
# histogram(ax3, test3, 'pourquoi')
# histogram(ax4, test4, 'okay')

# plt.tight_layout()
# plt.show()