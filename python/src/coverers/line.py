from src.coverers.data_structs import * 
from operator import attrgetter
from cv2 import imread, imshow, waitKey, destroyAllWindows
import numpy as np 
import matplotlib.pyplot as plt 
import math, os, glob

class Line: 
    
    def __init__(self, env:Environment, start:float, slope:float): 
        self.env = env 
        self.slope = slope
        #self.points = (1 / slope) * np.array([0] + self.env.radii) + start 
        self.points = (1 / slope) * np.array([0]+self.env.radii) + start 

    def plot(self, color): 
        # # Plot grey lines 
        # for radius in self.env.radii: 
        #     plt.plot([-self.env.top_layer_lim, self.env.top_layer_lim], 
        #              [radius, radius], 
        #              color=(0.5, 0.5, 0.5, 0.5), 
        #              linewidth=1)
        
        plt.plot(self.points, [0] + self.env.radii, c=color)

class LineGenerator(): 
    
    def __init__(self, env:Environment, start:float): 
        
        self.env = env 
        self.start = start 
        
        if not(-env.beam_axis_lim <= start <= env.beam_axis_lim): 
            raise Exception("Start point is out of range. ")
        
        max_height = env.radii[-1]
        self.slope_ll = max_height / (-env.top_layer_lim - start)
        self.slope_ul = max_height / (env.top_layer_lim - start)
        
    def generateGridLines(self, n=100): 
        
        angle_ll = math.atan(self.slope_ul)
        angle_ul = math.atan(self.slope_ll) + np.pi
        
        theta = np.linspace(angle_ul, angle_ll, n) 
        
        slopes = np.tan(theta) 
        return [Line(self.env, self.start, slope) for slope in slopes] 

    def generateEvenGrid(self, n=100):

        Rcoor = self.env.radii[-1]
        Zcoor = np.linspace(-self.env.top_layer_lim, self.env.top_layer_lim, n)
        
        slopes = Rcoor/(Zcoor-self.start)
        return [Line(self.env, self.start, slope) for slope in slopes] 
    
    def generateRandomGrid(self, n=100): 
        Rcoor = self.env.radii[-1]
        Zcoor = np.random.uniform(low=-self.env.top_layer_lim, high=self.env.top_layer_lim, size=n)
        
        slopes = Rcoor/(Zcoor-self.start)
        return [Line(self.env, self.start, slope) for slope in slopes] 
    
    def generateRandomLines(self, n=100): 
        
        angle_ll = math.atan(self.slope_ul)
        angle_ul = math.atan(self.slope_ll) + np.pi
        # print(angle_ll, angle_ul)
        
        theta = np.random.uniform(low=angle_ll, high=angle_ul, size=n) 
        
        slopes = np.tan(theta) 
        
        return [Line(self.env, self.start, slope) for slope in slopes] 
    
    def generateCenterSpreadLines(self, n=100): 
        
        angle_ll = math.atan(self.slope_ul)
        angle_ul = math.atan(self.slope_ll) + np.pi
        
        ranges = angle_ul - angle_ll
        
        angles = np.empty(shape=0)
        
        for exp in range(1, int(np.log2(n)) + 2):
            denom = 2 ** exp 
            angles = np.concatenate((angles, np.arange(1, denom, 2) / denom))
        
        angles = angles * ranges + angle_ll 
        
        slopes = np.tan(angles)
        
        return [Line(self.env, self.start, slope) for slope in slopes] 
    
    def generateCenterGridLines(self, n=100): 
        
        if n%2 == 0: 
            n += 1 
        
        angle_ll = math.atan(self.slope_ul)
        angle_ul = math.atan(self.slope_ll) + np.pi
        ranges = angle_ul - angle_ll 
        
        angles = np.linspace(-0.5, 0.5, n)
        angles = ranges * (np.array(sorted(angles, key=np.abs)) + 0.5) + angle_ll
        
        slopes = np.tan(angles)
        
        return [Line(self.env, self.start, slope) for slope in slopes]         