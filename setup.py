from setuptools import setup, find_packages

NAME = "GORi-KB"
VERSION = "0.0.9"
DESCRIPTION = "Context-informed geneset annotation."
LONG_DESCRIPTION = (
    "An algorithm to annotate hierarchies of genesets, with any knowledge base."
)
KEYWORDS = [
    "geneset annotation",
    "frequent itemset mining",
    "multiple resolutions",
    "knowledge base",
    "generic",
]
AUTHOR = "Yanis Asloudj"
AUTHOR_EMAIL = "yasloudj@u-bordeaux.fr"
URL = "https://github.com/yanisaspic/GORi-KB"

REQUIREMENTS = [
    "datetime",
    "gzip",
    "itables",
    "mlxtend",
    "networkx",
    "nxontology",
    "os",
    "pandas",
    "pipreqs",
    "plotly",
    "pypath",
    "pypath-omnipath",
    "re",
    "scipy",
    "setuptools",
    "shutil",
    "statsmodels",
    "time",
    "typing",
    "urllib.request",
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
    install_requires=REQUIREMENTS,
    packages=find_packages(),
)
