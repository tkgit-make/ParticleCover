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
data = DataSet(env, n_points=150, equal_spacing = False) 
cover = Cover(env, data)

#cover.solve(lining='solveS_improved', nlines=100)
#cover.solveS_center1()

#cover.plot(data=True, lines=True, patches=True)


acceptSlopePlot(lining='solveS_center2', savefig = True, ideal = False)
numCovers(lining='solveS_center2', savefig = True, ideal = False)
pointRepetitionFactor(lining='solveS_center2', savefig = True, ideal = False)

#pointRepetitionFactor(lining='solveS_relaxed_end', ideal = True, savefig=True)
#pointRepetitionFactor(lining='solveQ_relaxed_end', ideal = True, savefig=True)
#pointRepetitionFactor(lining='solveS_relaxed_gaps', ideal = True, savefig=True)


#final = cover.patches
#print(len(final), 'patches')
#print(cover.n_patches)
