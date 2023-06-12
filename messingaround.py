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
data = DataSet(env, n_points=150, equal_spacing = False) 
cover = Cover(env, data)
#data.plot()

#cover.solve(lining='solveS_relaxed_both', nlines=100, z0 = -15)
#cover.solve(lining='solveS_center2', nlines=100, z0 = -15)
#cover.plot(data=True, lines=True, patches=True)
'''
lGen = LineGenerator(env, -15)
fitting_lines = lGen.generateEvenGrid(100)
for line in fitting_lines:
    line.plot('r')

plt.ylim(-5, 30)
plt.show()
'''
def all_test(lining, ideal = False):

    acceptSlopePlot(lining=lining, savefig = True, ideal = ideal)
    numCovers(lining=lining, savefig = True, ideal = ideal)
    pointRepetitionFactor(lining=lining, savefig = True, ideal = ideal)


#duplicates('solveQ_relaxed_end', z0 = [-10, 0, 10], events = 1000)
fourTests(lining = 'solveQ_relaxed_both', solve_at = [-10, 0, 10], z0 = np.arange(-15, 15.5, 0.5), savefig = True)
#all_test('solveS', False)



#final = cover.patches
#print(len(final), 'patches')
#print(cover.n_patches)
