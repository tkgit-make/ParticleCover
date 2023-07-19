from src.readers.reader import * 
from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.coverers.wedgecover import *
from src.testers.test import *
import matplotlib.pyplot as plt 
import time 

filepath = "data/wedgeData_v3_128.txt"

filedata = readFile(filepath, stop=10, performance=False)
env, points = filedata[8] 
ds = DataSet(env)
ds.importData(points)
ds.add()
#ds.plot(True, True)
cov = wedgeCover(env, ds)
cov.solve('makePatches_Projective_quartile')
#cov.plot()
wedge_test(lining = 'q')