# GORi: An Alternative Approach for Genesets Enrichment Analyses

GORi is an algorithm developed to conduct enrichment analyses on a pool of genesets, with multiple knowledge bases, simultaneously. By employing the frequent itemset mining framework, GORi efficiently predicts annotation terms.
It also leverages the hierarchical relationships between genesets (*e.g.* genesets of clusters and their sub-clusters)
to improve its predictions. Finally, it automatically generates an interactive HTML report, to explore the results of the analysis.

Currently, GORi leverages eight semantically distinct concepts, spread across nine distinct *resources* (cf. GORi's tutorial):

- **biological processes (resource: GO_BP)**, from the Gene Ontology (GO).
- **cellular components (resource: GO_CC)**, from the Gene Ontology (GO).
- **cell types, with:**
    - **resource: CellMarker2**, from CellMarker 2.0 and the Cell Ontology (CL).
    - **resource: CellTaxonomy**, from CellTaxonomy and the Cell Ontology (CL).
- **diseases (resource: MeSH)**, from the Comparative Toxicogenomics Database (CTD) and the Medical Subject Headings (MeSH).
- **gene groups (resource: HGNC)**, from the HUGO Gene Nomenclature Committee (HGNC).
- **molecular functions (resource: GO_MF)**, from the Gene Ontology (GO).
- **pathways (resource: Reactome)**, from Reactome.
- **phenotypes (resource: HPO)**, from the Human Phenotype Ontology (HPO).

## Installation

You can install GORi from Github with:

```{shell}
pip3 install git+https://github.com/yanisaspic/GORi.git
```

## Overview of the GORi algorithm

- A tutorial on how to use GORi is available as a Jupyter Notebook in this repository: `./docs/gori.ipynb`.
- A demo dataset is also stored in this repository: `./data/Darmanis_HumGBM.xlsx`  

## Troubleshooting

GORi uses `pypath`, which writes cache files to reduce its loading times. Sometimes, cache files can be corrupted (*e.g.* when server requests are slowed and timed out). This can lead to errors downstream during the analysis. If you have any error, please make sure to wipe your cache:

```{python3}
from pypath.share import settings
settings.get('cachedir')    # your cache is stored in this directory
'home/yanisaspic/.cache/pypath'
```

```{shell}
rm /home/yanisaspic/.cache/pypath/*
```