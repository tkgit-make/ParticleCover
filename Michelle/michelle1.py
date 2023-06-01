import numpy as np 
import matplotlib.pyplot as plt
from data import *
from cover import *
from testing import *
import time as t

env = Environment()
data = DataSet(env, n_points = 150)
cover = Cover(env, data)

start_time = t.time()
cover.michelle()
print('Runtime:', t.time()-start_time)

#data.plot()
#print(data.array[0])
'''



cover.michelle()
print(cover.n_patches)

cover.plot()
'''

#numCovers()
#pointRepetitionFactor()