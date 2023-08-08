from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.coverers.wedgecover import *
from src.testers.test import *
from src.debug import * 
import matplotlib.pyplot as plt 
import time 
import ast

# filepath = "python/data/wedgeData_v3_128.txt"
# f = open("python/data/wedgeData_v3_128.txt")
#unaccepted_lines(apexZ0=[ 0], line_origin=[15], wedge_number=56, z0_cutoff=100.)
#filedata = readFile(filepath, stop=1280, performance=False)
minimal_cover_binary_search(lining = 'makePatches_Projective_Leftright', z0_spacing=0.5, ppl =16, z_5=50, wedges=6400, accept=0.999)
#minimal_cover_linear_search(lining = 'c',wedges=[0,128], accept=0.999)
# filedata = readFile(filepath, stop=1280, performance=False)
#wedge_test(lining = 'c',apexZ0=[-15,-10, 0,10, 15], top_layer_cutoff=50, wedges=[0,1], show_acceptance_of_cover=True)
'''

#filedata = readFile(filepath, stop=128, performance=False)
hist_data = [[],[],[],[],[]]
env = Environment(50)
plt.figure(figsize=(15, 7))
for i in range(1280):
    d = ast.literal_eval(f.readline())
    d = np.array(d)
    x_edges = np.array(env.radii)*(env.top_layer_lim-env.beam_axis_lim)/(env.radii[-1]) + env.beam_axis_lim
    count = np.zeros(5)
    for point in d:
        if np.abs(point[3]) <= x_edges[int(point[0]-1)]:
            count[int(point[0]-1)] += 1
    for j in range(5):
        hist_data[j].append(count[j])
for j in range(5):
    plt.subplot(1, 5, j+1)
    plt.hist(hist_data[j], bins = np.arange(20, 110, 10),edgecolor='black', rwidth=0.8, label = f"Average: {np.round(np.mean(hist_data[j]), 2)}\n std: {np.round(np.std(hist_data[j]), 2)}\n Density: {np.round(np.mean(hist_data[j])/x_edges[j], 2)}")
    plt.ylabel('Count', fontsize = 13)
    plt.xlabel('Number of Hits', fontsize = 13)
    plt.title(f'Layer {j+1}', fontsize = 15)
    plt.legend(fontsize = 13)
    print(np.mean(hist_data[j])/x_edges[j])
plt.suptitle('Points Between z = [-50, 50]', fontsize = 20)
plt.show()


#unaccepted_lines(z0_cutoff=50,  line_origin = [13], apexZ0= [-10, 0, 10], uniform_points=85)
#wedge_test(lining = 'c',wedges=[0,1], apexZ0=[-10,0,10], z0_cutoff=50., uniform_N_points = 85)

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
env = Environment(50)
ds = DataSet(env)
ds.importData(points)
ds.add()
#ds.plot(True)
cov = wedgeCover(env, ds)
cov.solve('makePatches_Projective_center', leftRight=False, apexZ0 = [-10, 0, 10])
#cov.plot()
'''
#pointRepetitionFactor()
#pointRepetitionFactorLayer(wedges=128, z_5=50, apexZ0=[-10,0,10])
# pointRepetitionFactorLayer(wedges=128, z_5=50)
# wedge_test(apexZ0=0)
# wedge_test(apexZ0 = [-10, 0, 10], wedges = [56, 57])

#wedge_test(lining= "makePatches_Projective_Rightleft", apexZ0 = [-10,0,10], wedges = [0,1280], top_layer_cutoff = 100, show_acceptance_of_cover=False, ppl = 16)
