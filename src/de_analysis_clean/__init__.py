# -*- coding: utf-8 -*-
__version__ = '0.1.0-dev'

__title__ = 'de-analysis'
__description__ = 'Distribution-free differential expression analysis for ' \
                  'scRNA-seq data across patient groups'
__url__ = 'https://github.com/erikadudki/de_analysis_clean'
__author__ = 'Erika Dudkin'
__email__ = 'erikadudkin@gmx.de'
__license__ = 'MIT License'
__copyright__ = 'Copyright (c) 2020 Erika Dudkin'

"""This file contains the initialization.

This package has several goals:
-...
-...

"""


# noqa = "no quality assurance"
from .get_perm_array_ind import get_perm_array_ind  # noqa:F401
from .filtering_cell_numbers import filtering_cell_numbers
from .create_patient_list import create_patient_list
from .anndata_to_tsv import anndata_to_tsv
import anndata
import pandas as pd
import numpy as np
import scipy.stats
import time
import math