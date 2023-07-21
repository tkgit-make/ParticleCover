from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.coverers.wedgecover import *
from src.testers.test import *
import matplotlib.pyplot as plt 
import time 
import ast

filepath = "data/wedgeData_v3_128.txt"
f = open('data/v3_patches.txt')
filedata = readFile(filepath, stop=128, performance=False)
wedge_test(lining = 'c',wedges=[0,1], apexZ0=[-10,0,10], uniform_N_points=122)

'''
env, points = filedata[0]
ds = DataSet(env)
ds.generateUniform([150, 150, 150, 150, 150])
ds.add()
#ds.plot(True, True)
cov = wedgeCover(env, ds)
cov.solve('makePatches_Projective_center', leftRight=False, apexZ0 = [-10, 0, 10])
cov.plot()

for i in range(128):
    d = ast.literal_eval(f.readline())
    d = np.array(d)
    d =d[:, :, :, 3].flatten()
    env, points = filedata[i] 
    ds = DataSet(env)
    ds.importData(points)
    ds.add()
    #ds.plot(True, True)
    cov = wedgeCover(env, ds)
    cov.solve('makePatches_Projective_center', apexZ0 = [-10, 0, 10])
    iter = 0
    d2 = []
    for patch in cov.patches:
        for sp in patch.superpoints:
            for p in sp.points:
                d2.append(p.z)
    if np.all(d == d2) == False:
        print('something is wrong')
    else:
        print('all good')
'''