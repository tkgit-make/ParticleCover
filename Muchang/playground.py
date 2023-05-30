import numpy as np 
from data import * 
from cover import * 

env = Environment()
data = DataSet(env, n_points=150) 
# print(data.array)
# data.plot(show_lines=True)
# plt.show()
cover = Cover(env, data) 
print(cover.n_patches) 
print(cover.patches)
cover.solve(100) 
print(cover.n_patches) 
print(cover.patches)
# cover.plot()