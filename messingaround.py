import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import *
from test_modules import *
import math
import cv2 
import os 
import glob
import matplotlib
import ast
env = Environment()
data = DataSet(env, n_points=150, equal_spacing = True) 
cover = Cover(env, data)
#data.plot()

#cover.solve(lining='solveS', nlines=100, z0 = 0, n = 16)
#print(cover.patches[0].superpoints[4].points)
#print(cover.patches[1].superpoints[4].points)
#cover.plot()
#print(cover.patches)
#print(cover.patches[0].superpoints[4].points)
#cover.solve(lining='solveS', nlines=100, z0 = 0, ideal = True)
#cover.plot(data=True, lines=True, patches=True)
#acceptSlopePlot(lining = 'solveS', ideal = True)
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
fourTests(lining = 'solveQ_relaxed_both', solve_at =  [-10, 0, 10], n = 32, z0 = np.arange(-15, 15.5, 0.5), savefig = True)
#all_test('solveS', False)

#wedge_data = np.loadtxt('Wedged_Spacepoints.txt', max_rows=1, usecols = np.arange(1699)+1)

#line = np.genfromtxt(r'Wedge_Data.txt', delimiter=',', max_rows = 1)

file = open("Wedge_Data.txt", 'r')
counter = np.zeros(5)
for i, line in enumerate(file.readlines()):
    d = np.array(ast.literal_eval(line))
    for point in d:
        layer = int(point[0])
        counter[layer-1] = counter[layer-1]+1
    if i == 127:
        break
print(counter)
print(counter/128)
'''
plt.scatter(d[:, 3],d[:, 1])
plt.plot([15, np.max(d[:, 3])], [0, 25], 'r')
plt.plot([-15, np.min(d[:, 3])], [0, 25], 'r')
plt.ylim(-2, 27)
plt.xlim(-200, 200)

plt.show()
'''