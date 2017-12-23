import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
import sys

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
