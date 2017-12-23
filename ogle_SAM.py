#!/usr/bin/env python

"""
Attempt to port Ogle's ANPP OpenBUGS example from Appendix 2. See Box 2 in the
main text.

Reference
---------
* Ogle et al. (2015) Quantifying ecological memory in plant and ecosystem
  processes. Ecology Letters, 18: 221–235
"""

import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

__author__  = "Martin De Kauwe"
__version__ = "1.0 (23.12.2017)"
__email__   = "mdekauwe@gmail.com"

N = 52
# Number of past years, including the current year for which the antecedent conditions are computed
Nlag = 5
Nyrs = 91
Nblocks = 38

# the time block that each month is assigned to such that for 60 different
# months, we are only estimating 38 unique monthly weights
block = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,\
                  18, 19, 20, 21, 22, 23, 24, 25, 25, 26, 26, 27, 27, 28, 28,\
                  29, 29, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, 33, 34, 34,\
                  34, 35, 35, 35, 36, 36, 36, 37, 37, 37, 38,\
                  38, 38]).reshape(5,12)




#with pm.Model() as model:
