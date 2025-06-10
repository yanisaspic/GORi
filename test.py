import gori

# we will also load two other dependencies required for this tutorial
import session_info
from urllib.request import urlretrieve


from gori.loaders import load_feve

feve_resource = load_feve(path="./data/Darmanis_HumGBM.xlsx)

from gori.init import download_resources, setup_resources
from gori.loaders import load_resources

resources = {"GO_BP", "GO_CC", "CellMarker2", "CellTaxonomy", "MeSH", "HGNC", "GO_MF", "Reactome", "HPO"}
resources_requiring_init = {"CellMarker2", "CellTaxonomy", "MeSH", "HGNC", "Reactome"}

download_resources(resources_requiring_init)
setup_resources(resources_requiring_init)
data = load_resources(resources)


from gori.loaders import load_local

dise_resource = load_local(resource="MeSH", path="./resources")


from gori.run import gori

data["fEVE"] = feve_resource
geneset = data["fEVE"]["annotations"].keys()
results = gori(geneset=geneset, antecedent_resource="fEVE", consequent_resources=resources, data=data, save=False)
print(results)



from gori.run import gorilon

results = gorilon(path)
print(results)




from gori.params import get_parameters

params = get_parameters()
params["n_genes_threshold"] = 10
params["pvalue_threshold"] = 0.001
params["use_heuristic"] = False
params["use_gene_symbol"] = False
params["stopwords"] = params["stopwords"] | {"bind, binding"}
results = gori(geneset=geneset, antecedent_resource="fEVE", consequent_resources={"CellMarker2", "GO_BP", "Reactome"}, data=data, params=params, save=False)
print(results)
