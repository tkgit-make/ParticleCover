import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import * 
import time 


env = Environment()
data = DataSet(env, n_points=150) 
cover1 = Cover(env, data) 
cover1.cluster("LR")
cover2 = Cover(env, data) 
cover2.cluster("C")
cover2.solveGridLR()
cover2.plot()


# print(cover1.superPoints[0])

# for i in range(5): 
#     print(len(cover1.superPoints[i]), len(cover2.superPoints[i]))


# for j in range(len(cover1.superPoints[1])): 
#     print(cover1.superPoints[1][j].min, cover1.superPoints[1][j].max) 
#     print(cover2.superPoints[1][j].min, cover2.superPoints[1][j].max) 



# cover.solveCenterGrid(100) 

# cover.plot(data=True, patches=False, lines=True)