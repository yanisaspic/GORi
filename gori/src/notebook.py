"""Functions called to generate the Jupyter Notebook for GORi.

    2025/05/29 @yanisaspic"""


import os
import nbformat
from typing import Any
from nbconvert import HTMLExporter
from traitlets.config import Config
from nbconvert.preprocessors import ExecutePreprocessor


def _get_notebook_template() -> str:
    """Get a Jupyter Notebook template for GORi.

    Returns
        A string with the Jupyter Notebook template, with placeholder values to replace.
    """
    template = """{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "c1ea0c03",
   "metadata": {},
   "source": [
    "Summary generated from: **PLACEHOLDER_SHEETS_PATH**.\n",
    "\n",
    "This summary includes interactive tables and plots: ***if you have any issue, please try refreshing or changing your browser.***"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "8c76f94f",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "\n",
    "# import itables to interact with data tables\n",
    "from itables import init_notebook_mode\n",
    "from itables import show\n",
    "\n",
    "init_notebook_mode(all_interactive=True)\n",
    "\n",
    "# import Plotly to interact with plots\n",
    "import plotly.express as px\n",
    "\n",
    "results = pd.read_excel('PLACEHOLDER_SHEETS_PATH', sheet_name=None)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "8e9e501f",
   "metadata": {},
   "source": [
    "## Gene annotations\n",
    "\n",
    "The boxplot and the table below report the number of annotations retrieved for each gene and resource.\n",
    "They include both direct and indirect annotations (*i.e.* annotations inferred from the hierarchical relationships of a direct annotation).\n",
    "\n",
    "**Note** that you can reach a webpage describing a gene by clicking on it in the table."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2f43b930",
   "metadata": {},
   "outputs": [],
   "source": [
    "# set-up the table of annotations\n",
    "annotations = results['annotations_counter'].copy()\n",
    "if PLACEHOLDER_USE_GENE_SYMBOL:\n",
    "    urls = [\n",
    "        f'https://pubchem.ncbi.nlm.nih.gov/gene/{g}/Homo_sapiens' for g in annotations.gene\n",
    "    ]\n",
    "else:\n",
    "    urls = [\n",
    "        f'https://www.uniprot.org/uniprotkb/{g}/entry' for g in annotations.gene\n",
    "    ]\n",
    "annotations.gene = [\n",
    "    '<a href=\\'{}\\' target=\\'_blank\\'>{}</a>'.format(u, g)\n",
    "    for u, g in zip(urls, annotations.gene)\n",
    "]\n",
    "\n",
    "init_notebook_mode(all_interactive=True)\n",
    "# data table\n",
    "show(\n",
    "    annotations,\n",
    "    layout={'topStart': 'search', 'topEnd': None},\n",
    "    classes='display nowrap cell-border compact',\n",
    "    columnDefs=[{'className': 'dt-left', 'targets': '_all'}],\n",
    "    style='table-layout:auto;width:100%;margin:auto',\n",
    "    allow_html=True\n",
    ")"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "0b81e48c",
   "metadata": {},
   "outputs": [],
   "source": [
    "# set-up the boxplot data\n",
    "boxplot_data = results['annotations_counter'].copy()\n",
    "boxplot_data = pd.melt(\n",
    "    boxplot_data, id_vars='gene', var_name='resource', value_name='# annotations'\n",
    ")\n",
    "boxplot = px.box(\n",
    "    boxplot_data,\n",
    "    x='resource',\n",
    "    y='# annotations',\n",
    "    color='resource',\n",
    "    points='all',\n",
    "    hover_data=['gene'],\n",
    "    color_discrete_sequence=px.colors.qualitative.Dark2,\n",
    ")\n",
    "\n",
    "boxplot.show(renderer='notebook')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "fb7e5c26",
   "metadata": {},
   "source": [
    "## Pairwise associations \n",
    "\n",
    "The barplot below report the number of association identified with each resource during the GORi analysis. These associations are progressively filtered out during the analysis, and can be split in four types:\n",
    "\n",
    "- **target associations** are associations between two resources; namely, a geneset you are trying to annotate, and an annotation term from a specific knowledge base.\n",
    "- **strong associations** are **target associations** with a high lift value (according to a heuristic analysis); they are likely to be statistically significant.\n",
    "- **unique associations** are **strong associations** after filtering out the redundant annotation terms; *e.g* `cell` is redundant with `basal cell`.\n",
    "- **corrected associations** are **unique associations** with a significant p-value after applying the Benjamini-Hochberg correction on unique associations; only corrected associations are returned by GORi.\n",
    "\n",
    "**Note** that the y-axis of the barplot is log10-transformed."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2f9f2410",
   "metadata": {},
   "outputs": [],
   "source": [
    "# set-up the barplot data\n",
    "barplot_data = results['associations_counter'].copy()\n",
    "\n",
    "total = barplot_data[barplot_data.columns[1:]].sum(axis=0).to_dict()\n",
    "tmp = {c: f'{c} ({v})' for c,v in total.items()}\n",
    "barplot_data = barplot_data.rename(tmp, axis=1)\n",
    "\n",
    "barplot_data = pd.melt(\n",
    "    barplot_data,\n",
    "    id_vars='resource',\n",
    "    var_name='association type',\n",
    "    value_name='# associations',\n",
    ")\n",
    "barplot = px.bar(\n",
    "    barplot_data,\n",
    "    x='association type',\n",
    "    y='# associations',\n",
    "    color='resource',\n",
    "    text_auto=True,\n",
    "    barmode='group',\n",
    "    log_y=True,\n",
    "    color_discrete_sequence=px.colors.qualitative.Dark2,\n",
    ")\n",
    "\n",
    "barplot.show(renderer='notebook')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "94eefd14",
   "metadata": {},
   "source": [
    "## Words of associations\n",
    "\n",
    "The table below reports the most characteristic words of each group annotated. For each group, the words are ranked according to their importance in the associations predicted by GORi. Each column is a group, and the rows are ranks, so that the first row includes the most characteristic word of each group, the second row the second most characteristic word, etc..."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "fec2d397",
   "metadata": {},
   "outputs": [],
   "source": [
    "# set-up the table of words\n",
    "words = results['words'].copy()\n",
    "words = words.set_index('word')\n",
    "_tmp = {}   # type: dict[str, pd.Series]\n",
    "for group in words.columns:\n",
    "    words = words.sort_values(by=[group, 'word'], ascending=[False, True])\n",
    "    _tmp[group] = words[group].loc[words[group] > 0].index.tolist()\n",
    "top_words = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in _tmp.items()]))\n",
    "top_words = top_words.fillna(\"\")\n",
    "\n",
    "# data table\n",
    "show(top_words,\n",
    "     layout={'topStart': 'search', 'topEnd': None},\n",
    "     paging=False, scrollY='350px',\n",
    "     classes='display nowrap cell-border compact',\n",
    "     columnDefs=[{'className': 'dt-left', 'targets': '_all'}])"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "3fd83faf",
   "metadata": {},
   "source": [
    "## List of associations\n",
    "\n",
    "The table below reports the associations predicted by GORi. Each row represents an association, and there are multiple columns:\n",
    "\n",
    "- **antecedents** and **consequents** are two concepts associated together.\n",
    "- **lift**, **n_genes** and **fdr** are metrics used to quantify the strength of an association.\n",
    "- **genes** is a list of comma-separated genes that are annotated by the **antecedents** and **consequents** simultaneously.\n",
    "- **resource_c** indicates the resource of the annotation terms in the column **consequents**.\n",
    "\n",
    "**Note** that you can reach a webpage defining an annotation term by clicking on it in the table."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "7eccdb17-b113-4e5d-ac08-128b0e430ba7",
   "metadata": {},
   "outputs": [],
   "source": [
    "# set-up the table of associations (in full)\n",
    "associations = results['associations'].copy()\n",
    "associations = associations.sort_values(by=['antecedents', 'lift'], ascending=[True, False])\n",
    "associations = associations.reset_index()\n",
    "associations.antecedents = [\n",
    "    a\n",
    "    if not isinstance(u, str)\n",
    "    else '<a href=\\'https://{}\\' target=\\'_blank\\'>{}</a>'.format(u, a)\n",
    "    for a, u in zip(associations.antecedents, associations.url_a)\n",
    "]\n",
    "associations.consequents = [\n",
    "    '<a href=\\'https://{}\\' target=\\'_blank\\'>{}</a>'.format(u, c)\n",
    "    for c, u in zip(associations.consequents, associations.url_c)\n",
    "]\n",
    "\n",
    "# data table\n",
    "show(\n",
    "    associations[['antecedents', 'consequents', 'lift', 'n_genes', 'fdr', 'genes',  'resource_c']],\n",
    "    layout={'top1': 'searchBuilder'},\n",
    "    style='table-layout:auto',\n",
    "    classes='display nowrap cell-border compact',\n",
    "    columnDefs=[{'className': 'dt-left', 'targets': '_all'}],\n",
    "    rowGroup={'dataSrc': 0},\n",
    "    allow_html=True\n",
    ")"
   ]
  }
 ],
 "metadata": {},
 "nbformat": 4,
 "nbformat_minor": 5
}"""
    template = template.replace('\n"', '\\n"')
    template = template.replace("\\'", "\\\\'")
    template = template.replace('""', '\\"\\"')
    return template


def _get_config() -> Config:
    """Get an HTMLExporter for a Jupyter Notebook.

    Returns
        A Config object."""
    c = Config()
    c.TemplateExporter.exclude_input = True
    c.TemplateExporter.exclude_output_prompt = True
    c.TemplateExporter.exclude_input_prompt = True
    return c


def _write_notebook(params: dict[str, Any]) -> None:
    """Execute a Jupyter Notebook for GORi.

    ``params`` is a dict of parameters.
    """
    template = _get_notebook_template()
    template = template.replace("PLACEHOLDER_SHEETS_PATH", params["sheets_path"])
    template = template.replace(
        "PLACEHOLDER_USE_GENE_SYMBOL", str(params["use_gene_symbol"])
    )

    hashcode = str(hash(params["sheets_path"] + params["report_path"]))
    ipynb_path = f"./.{hashcode}.ipynb"
    with open(ipynb_path, "w", encoding="utf-8") as f:
        f.write(template)

    with open(ipynb_path) as f:
        nb_in = nbformat.read(f, nbformat.NO_CONVERT)
    ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
    nb_out = ep.preprocess(nb_in)

    html_exporter = HTMLExporter(template_name="lab", config=_get_config())
    (body, _) = html_exporter.from_notebook_node(nb_out[0])

    with open(params["report_path"], "w", encoding="utf-8") as f:
        f.write(body)
    os.remove(ipynb_path)
