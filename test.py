from gori.params import get_parameters
from gori.loaders import load_priors
from gori.init import setup_iics
from gori.src.utils import _get_prior_inverse_translation
from gori.etc.lemmas import _get_lemmas_importances, _get_top_lemmas
from gori.init import setup_lemmas


priors = {"BIOP", "CTYP", "GENG", "PATH"}
params = get_parameters()
data = load_priors(priors)

import pandas as pd
setup_lemmas(data)
associations = pd.read_excel("./GORi.xlsx", sheet_name="associations")
test = _get_lemmas_importances(associations)
x = _get_top_lemmas(test)
print(x)