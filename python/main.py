from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.coverers.wedgecover import *
from src.testers.test import *
import matplotlib.pyplot as plt 
import time 
import ast

filepath = "data/wedgeData_v3_128.txt"
f = open("data/wedgeData_v3_128.txt")
#filedata = readFile(filepath, stop=128, performance=False)
hist_data = [[],[],[],[],[]]
plt.figure(figsize=(15, 12))
for i in range(1):
    d = ast.literal_eval(f.readline())
    d = np.array(d[:,3])
    for j in range(5):
        hist_data[j].append(d[:,j] <= 50)
for j in range(5):
    plt.subplot(1, 5, j+1)
    plt.hist(hist_data[j])
plt.show()

#unaccepted_lines(z0_cutoff=50,  line_origin = [13], apexZ0= [-10, 0, 10], uniform_points=85)
#wedge_test(lining = 'c',wedges=[0,1], apexZ0=[-10,0,10], z0_cutoff=50., uniform_N_points = 85)

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

env, points = filedata[0]
ds = DataSet(env)
num_points = 85
ds.generateUniform([num_points, num_points, num_points, num_points, num_points])
ds.add()
#ds.plot(True, True)
cov = wedgeCover(env, ds)
cov.solve('makePatches_Projective_center', leftRight=False, apexZ0 = [-10, 0, 10])
ends = []
for patch in cov.patches:
    ends.append(patch.end_layer)
plt.hist(ends, bins = np.arange(-0.5, 6.5, 1))
plt.show()
'''