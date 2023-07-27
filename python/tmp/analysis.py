from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.coverers.wedgecover import *
from src.testers.test import *
import matplotlib.pyplot as plt 
import time 

z0 = np.arange(-15, 15.5, 0.5)
v3_patch = np.loadtxt('data/solveS_center2_Acceptv3_transpose.csv', delimiter=',')
mask = np.abs(z0) <= 10.
a = np.mean(v3_patch[:, mask], axis = 1)
event = []

for i in range(7, 10):
    print(len(a[128*i:128*(i+1)]))
    #event.append(np.mean(a[128*i:128*(i+1)]))
    #event.append((a[128*i:128*(i+1)]))
    plt.hist(a[128*i:128*(i+1)], bins = np.arange(0.8, 1.02, 0.02), alpha = 0.5)
    #plt.hist(range(1280), a, s = 10)
    
    #plt.hist( np.mean(v3_patch, axis = 1), bins = 40)
#plt.scatter(range(10), event)
#plt.hist(a, bins = 40, label = f"mean: {np.round(np.mean(a), 3)*100}%\n stdev: {np.round(np.std(a)*100,3)}%")
plt.legend()
plt.show()