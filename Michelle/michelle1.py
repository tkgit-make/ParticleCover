import numpy as np 
import matplotlib.pyplot as plt
from data import *
from cover import *
from testing import *

env = Environment()
data = DataSet(env, n_points = 150)
cover = Cover(env, data)
#data.plot()
#print(data.array[0])
'''



cover.michelle()
print(cover.n_patches)

cover.plot()
'''

numCovers()
pointRepetitionFactor()