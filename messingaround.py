import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import *
from test_modules import *
import math
import cv2 
import os 
import glob

env = Environment()
data = DataSet(env, n_points=150) 
cover = Cover(env, data)
#cover.solveS()
cover.solveQ()
cover.plot(data=True, lines=True, patches=True)
#numCovers(clustering='', lining='solveQ')
final = cover.patches

colorset = ['orange', 'k', 'b', 'r', 'y', 'm']
shapeset = ['o', '^']
print(len(final), 'patches')
'''
for j in range(len(final)):
    for i in range(5):
    	points_to_plot = final[j].superpoints[i].points
    	plt.scatter(points_to_plot, np.full_like(points_to_plot, (i+1)*5), s = 15, color = colorset[j%6], marker = shapeset[j%2])
plt.xlim(-1.1, 1.1)
plt.ylim(0, 26)
plt.show()
'''