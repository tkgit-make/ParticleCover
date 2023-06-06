import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import *
from test_modules import *
import math
import cv2 
import os 
import glob

env = Environment()
data = DataSet(env, n_points=150) 
cover = Cover(env, data)

#cover.solve(lining='solveQ_relaxed_end', nlines=100)
#cover.solveS_center1()

#cover.plot(data=True, lines=True, patches=True)

acceptSlopePlot(lining='solveQ_relaxed_end', savefig = True)
numCovers(lining='solveQ_relaxed_end', savefig=True)
pointRepetitionFactor(lining='solveQ_relaxed_end', savefig=True)
pointRepetitionFactor(lining='solveQ_relaxed_end', savefig=True)


final = cover.patches
#print(len(final), 'patches')
#print(cover.n_patches)
