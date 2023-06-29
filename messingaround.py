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
f = open('wedgeData_v2_128.txt')
env = Environment()
#f.readline()
for i in range(7):
    line = f.readline()
d = np.array(ast.literal_eval(f.readline()))
data = DataSet(env = env, n_points = 150)
data.input_data(d, add = True)
cover = Cover(env, data)
cover.solve(z0 = 0, lining='solveS', n = 16, show = True)
cover.plot()
'''

env = Environment()
data = DataSet(env, n_points=150, equal_spacing = True) 
events = readFile('wedgeData_v2_128.txt', 8)
wedge1 = convertToDataset(events[2])
#cover = Cover(env, data)
cover = wedgeCover(env, wedge1)
cover.solve('solveQ', z0 = 0, show = True)
cover.plot()
#cover.plot()
#data.plot()


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

env = Environment()
data = DataSet(enqv = env, n_points = 150)
file = open('wedgeData_v2_128.txt')

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
first_wedge3.add()
first_wedge3.plot(True)
'''
#wedgeSlopePlot(lining = 'solveQ', z0 =-10)
#wedge_test_old(lining = 'solveS_reverse', solve_at = 0, n = 16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,128], savefig = False, v = 'v3')
#wedge_test_old(lining = 'solveQ', solve_at = [-10,0,10], n = 16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,128], savefig = False, v = 'v2')
#wedge_test_old(lining = 'solveS_center2', solve_at = 0, n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,128], savefig = False, v = 'v2')
#wedge_test(lining = 'solveS_reverse', solve_at = 0, n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,30], savefig = False, v = 'v3')
#wedge_test(lining = 'solveQ', solve_at = [-10,0,10], n =16, z0 = np.arange(-15, 15.5, 0.5), wedges = [0,128], savefig = True, v = 'v3')

def odd_loop(lining = 'solveS', v = 'v2'):
    accep = 0.0
    zvalues = 3
    while accep < 0.999:
        accep = wedge_test_old(lining = lining, solve_at = np.linspace(-15, 15, zvalues), n = 16, wedges = [0, 128], v = v)
        zvalues +=2
    print(accep)

def even_loop(lining = 'solveS', v = 'v2'):
    accep = 0.0
    zvalues = 2
    while accep < 0.999:
        accep = wedge_test_old(lining = lining, solve_at = np.linspace(-15, 15, zvalues), n = 16, wedges = [0, 128], v = v)
        zvalues +=2
    print(accep)
    

#odd_loop('solveS_reverse', v = 'v3')
#even_loop('solveS_reverse', v = 'v3')