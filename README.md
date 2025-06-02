# GORi: An Alternative Multi-Scale, Multi-Resolution Approach for Genesets Analyses

GORi is an algorithm developed to conduct enrichment analyses on a pool of genesets, with multiple knowledge bases, simultaneously. By employing the frequent itemset mining framework, GORi efficiently predicts annotation terms for a pool of genesets. It also takes into account their hierarchical relationships (*e.g.* genesets of clusters and their sub-clusters). Finally, it automatically generates an interactive HTML report to explore the results of the analysis.

Currently, eight types of annotations can be leveraged with GORi:

- **biological processes (BIOP)**, from the Gene Ontology (GO).
- **cellular components (CELC)**, from the Gene Ontology (GO).
- **cell types (CTYP)**, from CellMarker 2.0 and the Cell Ontology (CL).
- **diseases (DISE)**, from the Comparative Toxicogenomics Database (CTD) and the Medical Subject Headings (MeSH).
- **gene groups (GENG)**, from the HUGO Gene Nomenclature Committee (HGNC).
- **molecular functions (MOLF)**, from the Gene Ontology (GO).
- **pathways (PATH)**, from Reactome.
- **phenotypes (PHEN)**, from the Human Phenotype Ontology (HPO).

## Installation

You can install GORi from Github with:

```{shell}
pip3 install git+https://github.com/yanisaspic/GORi.git
```

## Overview of the GORi algorithm

- A tutorial on how to use GORi is available as a Jupyter Notebook on this repostiory: `./docs/gori.ipynb`.
- A demo dataset for the tutorial is also available on this repository: `./data/Darmanis_HumGBM.xlsx`  

## Troubleshooting

GORi uses `pypath`, which generates cache files to reduce loading times when running an analysis. Sometimes, cache files can be corrupted (*e.g.* when server requests are slowed and timed out). This can lead to errors downstream during the analysis. If you have any error, please make sure to wipe your cache:

```{python3}
from pypath.share import settings
settings.get('cachedir')    # your cache is stored in this directory
'home/yanisaspic/.cache/pypath'
```

```{shell}
rm /home/yanisaspic/.cache/pypath/*
```