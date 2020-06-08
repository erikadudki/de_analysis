# # -*- coding: utf-8 -*-
# __version__ = '0.1.0-dev'
#
# """This file contains the initialization.
# """
#
# import pandas as pd
# import numpy as np
# import scipy.stats
# import time
# import math
# import os
from itertools import combinations as comb
# noqa = "no quality assurance"
from .get_perm_array_ind import get_perm_array_ind  # noqa:F401
from .filtering_cell_numbers import filtering_cell_numbers
from .create_patient_list import create_patient_list
from .anndata_to_tsv import anndata_to_tsv
from .de_analysis import de_analysis
from .main_wilc_test import main_wilc_test
from .bool_run_on_grid import bool_run_on_grid
# import de_analysis.main_wilc_test
from de_analysis import as_numpy

# from os.path import dirname, basename, isfile, join
# import glob
# modules = glob.glob(join(dirname(__file__), "*.py"))
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]

# print(__all__)

# __all__ = ['de_analysis.py','as_numpy','create_patient_list',
#            'filtering_cell_numbers','get_perm_array_ind',
#            'main_wilc_test']