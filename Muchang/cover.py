import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
import math

class Line: 
    
    def __init__(self, env:Environment, start:float, slope:float): 
        self.env = env 
        self.slope = slope
        self.points = start + (env.radii / slope) *np.arange(1, env.layers + 1)

    def plot(self, color): 
        plt.plot(self.points, self.env.radii * np.arange(1, self.env.layers + 1), c=color)

class LineGenerator(): 
    
    def __init__(self, env:Environment, start:float): 
        
        self.env = env 
        self.start = start 
        
        if not(-env.bottom_layer_lim < start < env.bottom_layer_lim): 
            raise Exception("Start point is out of range. ")
        
        max_height = env.radii * env.layers
        # right_pnt = (env.top_layer_lim, max_height) 
        # left_pnt = (-env.top_layer_lim, max_height) 
        self.slope_ll = max_height / (-env.top_layer_lim - start)
        self.slope_ul = max_height / (env.top_layer_lim - start)
        # print(self.slope_ll, self.slope_ul) 
        
    def generateGridLines(self, n=100): 
        
        angle_ll = math.atan(self.slope_ul)
        angle_ul = math.atan(self.slope_ll) + np.pi
        # print(angle_ll, angle_ul) 
        
        theta = np.linspace(angle_ul, angle_ll, n) 
        
        slopes = np.tan(theta) 
        
        return [Line(self.env, self.start, slope) for slope in slopes] 
    
    def generateRandomLines(self, n=100): 
        
        angle_ll = math.atan(self.slope_ul)
        angle_ul = math.atan(self.slope_ll) + np.pi
        # print(angle_ll, angle_ul) 
        
        theta = np.random.uniform(low=angle_ll, high=angle_ul, size=n) 
        
        slopes = np.tan(theta) 
        
        return [Line(self.env, self.start, slope) for slope in slopes] 
    
        
class SuperPoint(): 
    
    def __init__(self, points:np.ndarray):
        
        if np.size(points) != 16: 
            raise Exception("This patch does not have 16 points in each layer")
        
        self.points = points 
        self.min = np.min(points) 
        self.max = np.max(points) 
        
    def contains(self, p:float): 
        return self.min < p < self.max
    
    def __eq__(self, other): 
        return (self.min, self.max) == (other.min, other.max)

class Patch(): 
    
    # Should be hashable (nvm we can't make it hashable) 
    
    def __init__(self, env:Environment, superpoints:tuple): 
        self.env = env 
        
        if len(superpoints) != env.layers: 
            raise Exception("The patch layers does not match environment layers. ")
        
        self.superpoints = superpoints
        # first superpoint in array should be the 1st layer 
        
    def contains(self, line:Line): 
        
        for i in range(len(self.superpoints)): 
            if not self.superpoints[i].contains(line.points[i]): 
                return False 
            
        return True
    
    def plot(self): 
        heights = np.arange(1, self.env.layers + 1) * self.env.radii 
        
        for i in range(self.env.layers): 
            sp = self.superpoints[i] 
            
            max_height = self.env.layers * self.env.radii
            
            plt.plot([sp.min, sp.max], [heights[i], heights[i]], c="g")
        
        plt.xticks(np.arange(-self.env.top_layer_lim, self.env.top_layer_lim, 0.1))
        plt.yticks(np.arange(0, max_height + 1, self.env.layers))
        # plt.show()
        
class Cover(): 
    
    def __init__(self, env:Environment, data:DataSet): 
        self.n_patches = 0 
        self.patches = [] 
        self.env = env 
        self.data = data 
        
        out = []
        array = self.data.array 
        
        for layer in array: 
            lyr = [] 
            end = len(layer) 
            
            i = 0 
            while i + 16 < end: 
                lyr.append(SuperPoint(layer[i:i+16]))
                i += 15 
            if i < end: 
                lyr.append(SuperPoint(layer[-16:]))
            
            out.append(lyr)
        
        self.superPoints = out
        
    
    def add_patch(self, curr_patch:Patch): 
        if self.n_patches == 0: 
            self.patches.append(curr_patch) 
            self.n_patches += 1 
        else:
            prev_patch = self.patches[-1] 
            prev_sp = prev_patch.superpoints 
            curr_sp = curr_patch.superpoints 
            
            for l in range(len(prev_sp)): 
                if (prev_sp[l].min != curr_sp[l].min) or (prev_sp[l].max != curr_sp[l].max): 
                    self.patches.append(curr_patch) 
                    self.n_patches += 1 
                    break 
                
            
        

    def solve(self, N=100): 
        lGen = LineGenerator(self.env, 0.0) 
        
        # num = 1
        for line in lGen.generateGridLines(N): 
            
            patch_ingredients = []
            
            indices = [0 for _ in range(self.env.layers)]
            
            for i in range(len(line.points)): 
                for j in range(indices[i], len(self.superPoints[i])): 
                    sp = self.superPoints[i][j]
                    if sp.contains(line.points[i]): 
                        patch_ingredients.append(sp) 
                        indices[i] = j
                        break 
                    
                    
            if len(patch_ingredients) == 5: 
                # line.plot(color="r") 
                pt = Patch(self.env, tuple(patch_ingredients))
                # pt.plot()
                self.add_patch(pt)
                # plt.savefig(f"Muchang/images2/{num}.png", dpi='figure', format=None)
                # plt.clf()
                # num+=1 
            else: 
                # line.plot(color="r")
                pass
        
        plt.show() 
                
            
    def plot(self): 
        for patch in self.patches: 
            patch.plot() 
        plt.show() 
        
