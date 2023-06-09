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
#data.plot()

cover.solve(lining='solveS', nlines=100, z0 = 0)

cover.plot(data=True, lines=True, patches=True)



def all_test(lining, ideal = False):

    acceptSlopePlot(lining=lining, savefig = True, ideal = ideal)
    numCovers(lining=lining, savefig = True, ideal = ideal)
    pointRepetitionFactor(lining=lining, savefig = True, ideal = ideal)


#all_test('solveS', False)
#acceptSlopePlotL(lining = 'solveS_center2')
#pointRepetitionFactor(lining='solveS_relaxed_end', ideal = True, savefig=True)
#pointRepetitionFactor(lining='solveQ_relaxed_end', ideal = True, savefig=True)
#pointRepetitionFactor(lining='solveS_relaxed_gaps', ideal = True, savefig=True)


final = cover.patches
print(len(final), 'patches')
#print(cover.n_patches)
