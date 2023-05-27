import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import * 
import time 

# x = [] 
# for i in range(10000): 
    
#     env = Environment()
#     data = DataSet(env, n_points=150) 
#     cover = Cover(env, data) 
#     cover.solve(100) 
#     x.append(cover.n_patches)

# plt.hist(x)
# plt.show() 

env = Environment()
data = DataSet(env, n_points=150) 
cover = Cover(env, data) 
cover.solve(100) 

num = 1
for patch in cover.patches: 
    patch.plot() 
    plt.savefig(f"Muchang/images3/{num}.png")
    plt.clf()
    num += 1