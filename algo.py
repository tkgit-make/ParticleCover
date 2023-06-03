import numpy as np 
from data import * 
from cover import * 

env = Environment()
data = DataSet(env, n_points=11, equal_spacing=False) 
# cover = Cover(env, data) 
# cover.solve(clustering="", lining="SlopeCenterStack1")
# cover.plot() 

