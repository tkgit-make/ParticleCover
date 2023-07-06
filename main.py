from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.cover import * 
import matplotlib.pyplot as plt 
import time 

filepath = "data/wedgeData_v3_128.txt"


filedata = readFile(filepath, stop=10, performance=False)
# filedata is a list of wedges 
# Each wedge is represented by tuple (env, pnts)
# env is an Environment, and pnts is a list of Points


env, points = filedata[0] 
ds = DataSet(env)
ds.importData(points) 
ds.plot() 
plt.show() 

# cov = Cover(ds) 

# cov.solveS() 