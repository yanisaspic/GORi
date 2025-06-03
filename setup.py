from setuptools import setup, find_packages

NAME = "gori"
VERSION = "0.1"
DESCRIPTION = "Context-informed genesets annotations."
LONG_DESCRIPTION = "An algorithm to annotate related genesets, with any knowledge base."
AUTHOR = "Yanis Asloudj"
AUTHOR_EMAIL = "yasloudj@u-bordeaux.fr"
URL = "https://github.com/yanisaspic/GORi"
KEYWORDS = [
    "geneset annotation",
    "frequent itemset mining",
    "multi-scale multi-resolution",
    "semantic enrichment",
    "single-cell",
]
PYTHON = ">=3.9"
REQUIREMENTS = [
    "beautifulsoup4>=4.12.3",
    "itables>=2.2.4",
    "mlxtend>=0.23.1",
    "networkx>=3.2",
    "notebook>=7.4.3",
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
    python_requires=PYTHON,
    install_requires=REQUIREMENTS,
    packages=find_packages(),
)
