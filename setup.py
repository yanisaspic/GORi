from setuptools import setup, find_packages

NAME = "GORi-KB"
VERSION = "0.1"
DESCRIPTION = "Context-informed genesets annotation."
LONG_DESCRIPTION = "An algorithm to annotate related genesets, with any knowledge base."
KEYWORDS = [
    "geneset annotation",
    "frequent itemset mining",
    "multi-scale multi-resolution",
    "semantic enrichment",
    "single-cell",
]
AUTHOR = "Yanis Asloudj"
AUTHOR_EMAIL = "yasloudj@u-bordeaux.fr"
URL = "https://github.com/yanisaspic/GORi-KB"

REQUIREMENTS = [
    "itables>=2.2.4",
    "mlxtend>=0.23.1",
    "networkx>=3.2",
    "nxontology>=0.5.0",
    "pandas>=1.4.4",
    "pipreqs>=0.5.0",
    "plotly>=5.19.0",
    "pypath>=0.1",
    "pypath_common>=0.2.5",
    "pypath-omnipath>=0.16.17",
    "scipy>=1.9.1",
    "setuptools>=69.2.0",
    "statsmodels>=0.13.2",
]

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    keywords=KEYWORDS,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    url=URL,
    python_requires=">=3.9",
    install_requires=REQUIREMENTS,
    packages=find_packages(),
)
