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
import time
from reader import *
from converter import *
from wedgecover import *
'''
env = Environment()
data = DataSet(env, n_points=150, equal_spacing = True) 
events = readFile('wedgeData_v2_128.txt', 5)
wedge1 = convertToDataset(events[0])
#cover = Cover(env, data)
cover = wedgeCover(env, wedge1)
cover.solve('solveS', z0 = [-10,0,10], show = False)
cover.plot()
#cover.plot()
#data.plot()
'''
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

#duplicates('solveQ_relaxed_end', z0 = [-10, 0, 10], events = 1000)
#fourTests(lining = 'solveS_center2', solve_at =  [-10, 0, 10], n = 16, z0 = np.arange(-15, 15.5, 0.5), savefig = False)
#all_test('solveS', False)

#wedge_data = np.loadtxt('Wedged_Spacepoints.txt', max_rows=1, usecols = np.arange(1699)+1)

#line = np.genfromtxt(r'Wedge_Data.txt', delimiter=',', max_rows = 1)

'''
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

file.readline()
line = file.readline()
d = np.array(ast.literal_eval(line))

#plt.show()
data.input_data(d)
cover = Cover(env, data)
cover.solve(lining = 'solveS_center2')
cover.plot()
acceptSlopePlot(events=1, custom = d)


for i, line in enumerate(file.readlines()):
    z0test = []
    count = []
    acc = []
    d = np.array(ast.literal_eval(line))
    data.input_data(d)
    cover = Cover(env, data)
    cover.solve(lining = 'solveS_center2', z0 =0, show = True)
    count.append(cover.n_patches)
    cover.plot()
    
    acc.append(acceptSlopePlot(events=1, custom = d, show = False))
    for z0 in np.arange(-15, 15.5, 0.5):
        z0test.append(acceptSlopePlot(events=1, custom = d, show = False, solve_at=0, z0 = z0))
    print(f'acc: {np.mean(acc)}')
    print(f'mean: {np.mean(count)}')
    plt.scatter(np.arange(-15, 15.5, 0.5), z0test, color = 'r')
    plt.plot(np.arange(-15, 15.5, 0.5), z0test, color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance Rate',  fontsize = 16)
    plt.ylim(0, 1.0)
    plt.title(f'z0 Offset vs Acceptance Rate for Wedge 1', fontsize = 16)
    plt.show()
    plt.clf()
    
    if i == 5:
        break


line = file.readline()
d = np.array(ast.literal_eval(line))
data.input_data(d)
cover = Cover(env, data)
cover.solve(lining = 'solveS', z0 = 0,  show = False)
cover.plot()
print(acceptSlopePlot(events=1, custom = d, show = False))
file = open("Wedge_Datav2.txt", 'r')

events2= readFile('wedgeData_v2.1_128.txt', 5)
events1 = readFile('wedgeData_v2_128.txt', 5)
first_wedge2 = convertToDataset(events2[0]).array

first_wedge1 = convertToDataset(events1[0]).array
for i in range(len(first_wedge1)):
    for j in range(len(first_wedge1[i])):
        if first_wedge1[i][j].z != first_wedge2[i][j].z:
            print('not same')
            break
print('same')

file = open('wedgeData_v2.1_128.txt', 'r')
line1 = file.readline()
d = np.array(ast.literal_eval(line1))
file = open('wedgeData_v2_128.txt', 'r')
line1 = file.readline()
f = np.array(ast.literal_eval(line1))
print(np.all(d == f))
events2= readFile('wedgeData_v2.1_128.txt', 5)




for i in range(128):
    first_wedge3 = convertToDataset(events3[i])
    print(first_wedge3.array[0][0].phi)

events3 = readFile('wedgeData_v3_128.txt', 128)
first_wedge3 = convertToDataset(events3[1])
#first_wedge3.plot(True)
#first_wedge3.add()
#first_wedge3.plot(True)
'''

wedge_test(lining = 'solveS_center2', solve_at = [-10,0,10], n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,2], savefig = False, v = 'v3')
wedge_test(lining = 'solveS', solve_at = 0, n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,1], savefig = False, v = 'v3')
#wedge_test(lining = 'solveS_reverse', solve_at = 0, n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,30], savefig = False, v = 'v3')
#wedge_test(lining = 'solveQ', solve_at = [-10,0,10], n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,128], savefig = True, v = 'v3')