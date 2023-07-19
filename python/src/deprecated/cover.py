import numpy as np 
import matplotlib.pyplot as plt 
import math
import cv2 
import os 
import glob 
from data import * 

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
        
        if not(-env.bottom_layer_lim <= start <= env.bottom_layer_lim): 
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

    def generateEvenGrid(self, n=100):

        ycoor = self.env.layers*self.env.radii #set to radius of outermost layer when list of radii is avaliable
        xcoor = np.linspace(-self.env.top_layer_lim, self.env.top_layer_lim, n)
        
        slopes = ycoor/(xcoor-self.start)
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
        
class SuperPoint(): 
    
    def __init__(self, points:np.ndarray):
        
        if np.size(points) != 16:
            if np.size(points) != 32:
                raise Exception("This patch does not have 16 or 32 points in each layer")
        
        self.points = points 
        self.min = np.min(points) 
        self.max = np.max(points) 
        
    def contains(self, p:float): 
        return self.min <= p <= self.max
    
    def __eq__(self, other): 
        return (self.min, self.max) == (other.min, other.max)

class Patch(): 
    
    # Should be hashable (nvm we can't make it hashable) 
    
    def __init__(self, env:Environment, superpoints:tuple): 
        self.env = env
        self.end = 0
        
        if len(superpoints) != env.layers: 
            raise Exception("The patch layers does not match environment layers. ")
        
        self.superpoints = superpoints
        # first superpoint in array should be the 1st layer 
        
    def contains(self, line:Line): 
        
        for i in range(len(self.superpoints)): 
            if not self.superpoints[i].contains(line.points[i]): 
                return False 
            
        return True

    def contains_p(self, point:float, layer:int): 
        
        sp = self.superpoints[layer] 
        return sp.contains(point)

    def end_layer(self, layer = 0):
        self.end = int(layer)
    
    def __eq__(self, other): 
        if not isinstance(other, Patch): 
            return NotImplemented 
        
        for i in range(len(self.superpoints)): 
            if self.superpoints[i] != other.superpoints[i]: 
                return False 
            
        return True 
        
    def plot(self, color='g'): 
        heights = np.arange(1, self.env.layers + 1) * self.env.radii 
        
        for i in range(self.env.layers): 
            sp = self.superpoints[i] 
            
            max_height = self.env.layers * self.env.radii
            
            plt.plot([sp.min, sp.max], [heights[i], heights[i]], c=color)
        
        plt.xticks(np.arange(-self.env.top_layer_lim, self.env.top_layer_lim, 10))
        plt.yticks(np.arange(0, max_height + 1, self.env.layers))
        # plt.show()
        
def two_digit(number:int): 
    if 0 <= number < 10: 
        return f"0{number}" 
    else: 
        return f"{number}"
    
class Cover(): 
    
    def __init__(self, env:Environment, data:DataSet): 
        self.n_patches = 0 
        self.patches = [] 
        self.env = env 
        self.data = data 
        self.fitting_lines = [] 
        self.superPoints = [] 
        
    def cluster(self, type): 
        # type can be LR or C
        
        out = []
        array = self.data.array 
        
        if type == "LR": 
            
            for layer in array: 
                lyr = [] 
                end = len(layer) 
                
                i = 0 
                while i + 16 < end: 
                    lyr.append(SuperPoint(layer[i:i+16]))
                    i += 14 
                if i < end: 
                    lyr.append(SuperPoint(layer[-16:]))
                
                
                out.append(lyr)
        
        elif type == "C": 
            
            for layer in array: 
                
                lyr = []
                end = len(layer)
                
                first_min = 0
                for i in range(len(layer)): 
                    if layer[i] < 0: 
                        first_min += 1
                    elif layer[i] > 0: 
                        break 
                first_min = max(first_min - 8, 0)
                
                j = first_min 
                while j + 16 < end: 
                    lyr.append(SuperPoint(layer[j:j+16]))
                    j += 14
                if j < end: 
                    lyr.append(SuperPoint(layer[-16:]))
                    
                j = first_min 
                while j - 16 >= 0: 
                    lyr.append(SuperPoint(layer[j-14:j+2]))
                    j -= 14 
                if j > 0: 
                    lyr.append(SuperPoint(layer[:16]))
                
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
            
    def solveGridLR(self, nlines=100): 
        # Decent, with 38 patch performance on average
        
        if self.superPoints == []: 
            raise Exception("You must cluster the superpoints first by running the cluster method.")
        
        lGen = LineGenerator(self.env, 0.0) 
        self.fitting_lines = lGen.generateGridLines(nlines)
        
        for line in self.fitting_lines: 
            
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
                self.add_patch(Patch(self.env, tuple(patch_ingredients)))
        
    def solveRandomized(self, nlines=100): 
        
        if self.superPoints == []: 
            raise Exception("You must cluster the superpoints first by running the cluster method.")
        
        lGen = LineGenerator(self.env, 0.0) 
        self.fitting_lines = lGen.generateRandomLines(nlines)
        
        for line in self.fitting_lines: 
            
            patch_ingredients = []
            
            for i in range(len(line.points)): 
                for j in range(len(self.superPoints[i])): 
                    sp = self.superPoints[i][j]
                    if sp.contains(line.points[i]): 
                        patch_ingredients.append(sp) 
                        break 
                    
                    
            if len(patch_ingredients) == 5: 
                pt = Patch(self.env, tuple(patch_ingredients))
                
                contains = False 
                for patch in self.patches: 
                    if pt == patch: 
                        contains = True 
                        break 
                
                if contains == False: 
                    self.add_patch(pt)
        
    def solveCenterSpread(self, nlines=100): 
        
        if self.superPoints == []: 
            raise Exception("You must cluster the superpoints first by running the cluster method.")
        
        lGen = LineGenerator(self.env, 0.0) 
        self.fitting_lines = lGen.generateCenterSpreadLines(nlines)
        
        for line in self.fitting_lines: 
            
            patch_ingredients = []
            
            for i in range(len(line.points)): 
                for j in range(len(self.superPoints[i])): 
                    sp = self.superPoints[i][j]
                    if sp.contains(line.points[i]): 
                        patch_ingredients.append(sp) 
                        break 
                    
                    
            if len(patch_ingredients) == 5: 
                pt = Patch(self.env, tuple(patch_ingredients))
                
                contains = False 
                for patch in self.patches: 
                    if pt == patch: 
                        contains = True 
                        break 
                
                if contains == False: 
                    self.add_patch(pt)
                    
    def solveCenterGrid(self, nlines=100): 
        
        if self.superPoints == []: 
            raise Exception("You must cluster the superpoints first by running the cluster method.")
        
        lGen = LineGenerator(self.env, 0.0) 
        self.fitting_lines = lGen.generateCenterGridLines(nlines)
        
        for line in self.fitting_lines: 
            
            patch_ingredients = []
            
            for i in range(len(line.points)): 
                for j in range(len(self.superPoints[i])): 
                    sp = self.superPoints[i][j]
                    if sp.contains(line.points[i]): 
                        patch_ingredients.append(sp) 
                        break 
                    
                    
            if len(patch_ingredients) == 5: 
                pt = Patch(self.env, tuple(patch_ingredients))
                
                contains = False 
                for patch in self.patches: 
                    if pt == patch: 
                        contains = True 
                        break 
                
                if contains == False: 
                    self.add_patch(pt)

    def solveS_old1(self, stop = 1): #NOT USED ANYMORE

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            init_patch.append(SuperPoint(self.data.array[row, 0:16]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_repeated(stop)
        return
    
    def solveS_reverse_old1(self, stop = -1): #NOT USED ANYMORE

        init_patch = []

        #pick last 16 points in each layer
        for row in range(5):
            init_patch.append(SuperPoint(self.data.array[row, -16:]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_repeated_reverse(stop = stop)
        return
        
    def S_repeated(self, stop = 1): #NOT USED ANYMORE

        loops = self.n_patches - 1
        mins = []
        stop_count = 0

        #find point in patch with largest slope
        for i in range(5):

            #last_patch looks at the previously generated patch
            last_patch = self.patches[loops].superpoints

            #finds the rescaled point for each layer
            amin = last_patch[i].points[15]/((i+1)/20)
            mins.append(amin)

        #find the layer with the lowest re-scaled point value
        min_index = np.argmin(np.array(mins))
        patch_ingredients = []

        #lowest re-scaled point value
        min_value = last_patch[min_index].points[15]/(min_index+1)*20

        #loops through each layer again
        for i in range(5):

            #Goes through the re-scaled value of each layer and find the re-scaled value that is closest to min_value
            closest_index = np.argmin(np.abs(self.data.array[i]/((i+1)/20) - min_value))

            #picks the index of the stop value
            #stop_rescaled = stop*(0.85*(i+1)/5+0.15)
            #stop_index = np.argmin(np.abs(self.data.array[i] - stop_rescaled))
            stop_rescaled = stop*(i+1)*20
            stop_index = np.argmin(np.abs(self.data.array[i] - stop_rescaled))
            stop_value = self.data.array[i][stop_index]

            #if the closest index is the first point then pick the first 16 points for a superpoint
            if closest_index < 1:
                patch_ingredients.append(SuperPoint(self.data.array[i, 0:16]))

            #if there isn't enough point remaining pick the last 16 points
            elif closest_index > (self.data.n_points) - 16:
                patch_ingredients.append(SuperPoint(self.data.array[i, -16:]))
            
            #otherwise pick 16 points starting from the closest_index-1
            else:
                patch_ingredients.append(SuperPoint(self.data.array[i, closest_index-1:closest_index+15]))

            if np.any(patch_ingredients[i].points > stop_value):
                stop_count += 1

        #generate new patch from those new superpoints   
        new_patch = Patch(self.env, tuple(patch_ingredients))
        #check if the new patch is the same as the last patch; if so, terminate
        if np.array_equal(new_patch.superpoints, last_patch):
            return

        #otherwise add the patch
        else:
            self.add_patch(new_patch)
            if stop_count == 5:
                return
            else:
                return self.S_repeated(stop)

    def S_repeated_reverse(self, stop = -1): #NOT USED ANYMORE

        loops = self.n_patches - 1
        mins = []
        stop_count = 0

        #find point in patch with largest slope
        for i in range(5):

            #last_patch looks at the previously generated patch
            last_patch = self.patches[loops].superpoints

            #finds the rescaled point for each layer
            amin = last_patch[i].points[0]/((i+1)/20)
            mins.append(amin)

        #find the layer with the lowest re-scaled point value
        min_index = np.argmax(np.array(mins))
        patch_ingredients = []

        #lowest re-scaled point value
        min_value = last_patch[min_index].points[0]/(min_index+1)*20

        #loops through each layer again
        for i in range(5):

            #Goes through the re-scaled value of each layer and find the re-scaled value that is closest to min_value
            closest_index = np.argmin(np.abs(self.data.array[i]/((i+1)/20) - min_value))

            #stop_rescaled = stop*(0.85*(i+1)/5+0.15)
            #stop_index = np.argmin(np.abs(self.data.array[i] - stop_rescaled))
            stop_rescaled = stop*(i+1)*20
            stop_index = np.argmin(np.abs(self.data.array[i] - stop_rescaled))
            stop_value = self.data.array[i][stop_index]

            #if the closest index is the first point then pick the first 16 points for a superpoint
            if closest_index > (self.data.n_points) - 2:
                patch_ingredients.append(SuperPoint(self.data.array[i, -16:]))

            #if there isn't enough point remaining pick the last 16 points
            elif closest_index < 15:
                patch_ingredients.append(SuperPoint(self.data.array[i, :16]))
            
            #otherwise pick 16 points starting from the closest_index-1
            else:
                patch_ingredients.append(SuperPoint(self.data.array[i, closest_index-14:closest_index+2]))

            if np.any(patch_ingredients[i].points < stop_value):
                stop_count += 1

        #generate new patch from those new superpoints   
        new_patch = Patch(self.env, tuple(patch_ingredients))
        #check if the new patch is the same as the last patch; if so, terminate
        if np.array_equal(new_patch.superpoints, last_patch):
            return

        #otherwise add the patch
        else:
            self.add_patch(new_patch)
            if stop_count == 5:
                return
            else:
                return self.S_repeated_reverse(stop = stop)

    def solveS_old2(self, stop = 1): #NOT USED ANYMORE

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            start_index = np.argmin(np.abs(self.data.array[row] + (row+1)*20))
            if start_index != 0:
                start_index -= 1
            init_patch.append(SuperPoint(self.data.array[row, start_index:start_index+16]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_loop_improved(stop = 1)
        return

    def solveS_reverse_old2(self, stop = -1): #NOT USED ANYMORE

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            start_index = np.argmin(np.abs(self.data.array[row] - (row+1)*20))
            if start_index != 149:
                start_index += 1
            init_patch.append(SuperPoint(self.data.array[row, start_index-15:start_index+1]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_rloop_improved(stop = -1)
        return

    def solveS(self, z0 = 0, stop = 1, n = 16):

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            y = 5*(row + 1)
            start_index = np.argmin(np.abs(self.data.array[row] - (((-z0-100)*y)/25+z0)))
            if start_index != 0:
            #if self.data.array[row][start_index] > (((-z0-100)*y)/25+z0):
                start_index -= 1
            init_patch.append(SuperPoint(self.data.array[row][start_index:start_index+n]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_loop15(z0=z0, stop=stop, n=n)
        return

    def solveS_reverse(self, z0 = 0, stop = -1, n = 16):

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            y = 5*(row+1)
            start_index = np.argmin(np.abs((self.data.array[row] - ((100-z0)*y/25 + z0))))

            if start_index != len(self.data.array[row])-1:
            #if self.data.array[row][start_index] < ((100-z0)*y/25 + z0):
                start_index += 1
            init_patch.append(SuperPoint(self.data.array[row][start_index-n+1:start_index+1]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_rloop15(z0 = z0, stop = stop, n = n)
        return

    def solveS_center1(self, n = 16):
        init_patch = []

        #pick center 16 points based on index
        for row in range(5):
            init_patch.append(SuperPoint(self.data.array[row, int(self.data.n_points/2)-int(n/2):int(self.data.n_points/2)+int(n/2)]))
        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_loop_improved()
        self.add_patch(Patch(self.env, tuple(init_patch)))
        self.S_rloop_improved()
        self.patches = self.patches[1:]
        self.n_patches = self.n_patches - 1
        return

    def solveS_center2_old(self, center = 0, stop = 'none', n = 16): #NOT USED ANYMORE
        init_patch = []

        #pick center 16 points based on 
        for row in range(5):
            center_index = np.argmin(np.abs(self.data.array[row] - center*(row+1)*20))
            init_patch.append(SuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))
        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_loop_improved(stop = 0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_rloop_improved(n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
            #starts right of center
            if center > 0:
                n_patch_start = self.n_patches
                self.S_loop_improved(n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_rloop_improved(stop = 0, n =n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            #run main algorithm
            n_patch_start = self.n_patches
            self.S_loop_improved(n = n)
            self.add_patch(Patch(self.env, tuple(init_patch)))
            self.S_rloop_improved(n = n)
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 1
            return

    def solveS_center2(self, center = 0, z0 = 0, stop = 'none', n = 16):
        init_patch = []

        #pick center 16 points based on 
        for row in range(5):
            y = 5*(row+1)
            center_index = np.argmin(np.abs(self.data.array[row] - ((y*(center-z0)/25)+z0)))
            
            if (center_index-int(n/2)) < 0:
                center_index = int(n/2)
            elif (center_index+int(n/2)) > len(self.data.array[row]):
                center_index = len(self.data.array[row]) - int(n/2)
            '''
            if self.data.array[row][center_index] >= 0:
                init_patch.append(SuperPoint(self.data.array[row][center_index-int(n/2):center_index+int(n/2)]))
            else:
                init_patch.append(SuperPoint(self.data.array[row][center_index-int(n/2)+1:center_index+int(n/2)+1]))
            '''

            init_patch.append(SuperPoint(self.data.array[row][center_index-int(n/2):center_index+int(n/2)]))
        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, stop = 0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_rloop15(z0 = z0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
            #starts right of center
            if center > 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_rloop15(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            #run main algorithm
            n_patch_start = self.n_patches
            self.S_loop15(z0 = z0, n = n)
            self.add_patch(Patch(self.env, tuple(init_patch)))
            self.S_rloop15(z0 = z0, n = n)
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 1
            return
        
    def solveQ(self, z0 = 0, n = 16):
        self.solveS_center2(center = -50, stop = 'center', z0 = z0, n = n) 
        self.solveS_center2(center = 50, stop = 'center', z0 = z0, n = n)
        return

    def solveS_relaxed_end(self, center = 0, z0 = 0, stop = 'none', n = 16):
        init_patch = []

        #pick center 16 points based on 
        for row in range(5):
            y = 5*(row+1)
            center_index = np.argmin(np.abs(self.data.array[row] - ((y*(center-z0)/25)+z0)))
            if (center_index-int(n/2)) < 0:
                center_index = int(n/2)
            elif (center_index+int(n/2))+1 > len(self.data.array[row]):
                center_index = len(self.data.array[row])-int(n/2)
            init_patch.append(SuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))
        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, stop = 0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_rloop15(z0 = z0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 2
                del self.patches[self.n_patches]
                return
            #starts right of center
            if center > 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, n = n)
                del self.patches[self.n_patches-1]
                self.n_patches = self.n_patches - 1
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_rloop15(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1 
                return
        else:

            #run main algorithm
            n_patch_start = self.n_patches
            self.S_loop15(z0 = z0, n = n)
            del self.patches[self.n_patches-1]
            self.n_patches = self.n_patches - 1
            self.add_patch(Patch(self.env, tuple(init_patch)))
            self.S_rloop15(z0 = z0, n = n)
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 2
            del self.patches[self.n_patches]
            return

    def solveS_relaxed_gap(self, center = 0, z0 = 0, stop = 'none', n = 16):
        init_patch = []

        for row in range(5):
            y = 5*(row+1)
            center_index = np.argmin(np.abs(self.data.array[row] - ((y*(center-z0)/25)+z0)))
            if (center_index-int(n/2)) < 0:
                center_index = int(n/2)
            elif (center_index+int(n/2))+1 > len(self.data.array[row]):
                center_index = len(self.data.array[row])-int(n/2)
            init_patch.append(SuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))

        self.add_patch(Patch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_relaxed_gap(z0 = z0, stop = 0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_reverse_relaxed_gap(z0 = z0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
            #starts right of center
            if center > 0:
                n_patch_start = self.n_patches
                self.S_relaxed_gap(z0 = z0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_reverse_relaxed_gap(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            #run main algorithm
            n_patch_start = self.n_patches
            self.S_relaxed_gap(z0 = z0, n = n)
            self.add_patch(Patch(self.env, tuple(init_patch)))
            self.S_reverse_relaxed_gap(z0 = z0, n = n)
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 1
            return

    def solveS_relaxed_both(self, center = 0, z0 = 0, stop = 'none', n = 16):
        init_patch = []

        for row in range(5):
            y = 5*(row+1)
            center_index = np.argmin(np.abs(self.data.array[row] - ((y*(center-z0)/25)+z0)))
            if (center_index-int(n/2)) < 0:
                center_index = int(n/2)
            elif (center_index+int(n/2))+1 > len(self.data.array[row]):
                center_index = len(self.data.array[row])-int(n/2)
            init_patch.append(SuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))

        self.add_patch(Patch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_relaxed_gap(z0 = z0, stop = 0, n = n)
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_reverse_relaxed_gap(z0 = z0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 2
                del self.patches[self.n_patches]
                return
            #starts right of center
            if center > 0:
                n_patch_start = self.n_patches
                self.S_relaxed_gap(z0 = z0, n = n)
                del self.patches[self.n_patches-1]
                self.n_patches = self.n_patches - 1
                self.add_patch(Patch(self.env, tuple(init_patch)))
                self.S_reverse_relaxed_gap(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            n_patch_start = self.n_patches
            self.S_relaxed_gap(z0 = z0, n = n)
            del self.patches[self.n_patches-1]
            self.add_patch(Patch(self.env, tuple(init_patch)))
            self.n_patches = self.n_patches - 1
            self.S_reverse_relaxed_gap(z0 = z0, n = n)
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 2
            del self.patches[self.n_patches]
            return

    def solveQ_relaxed_end(self, z0 = 0, n = 16):
        self.solveS_relaxed_end(center = -50, z0 = z0, stop = 'center', n = n) 
        self.solveS_relaxed_end(center = 50,  z0 = z0, stop = 'center', n = n)
        return

    def solveQ_relaxed_gap(self, z0 = 0, n = 16):
        self.solveS_relaxed_gap(center = -50,  z0 = z0, stop = 'center', n = n) 
        self.solveS_relaxed_gap(center = 50,  z0 = z0, stop = 'center', n = n)
        return

    def solveQ_relaxed_both(self,z0 = 0, n = 16):
        self.solveS_relaxed_both(center = -50,  z0 = z0, stop = 'center', n = n) 
        self.solveS_relaxed_both(center = 50,  z0 = z0, stop = 'center', n = n)
        return

    def S_relaxed_gap(self, z0 = 0, stop = 1, n = 16):
        
        loops = self.n_patches - 1
        mins = []

        for i in range(5):
            y = 5*(i+1)
            last_patch = self.patches[loops].superpoints

            amin = (last_patch[i].points[n-1]-z0)/(y/100)
            mins.append(amin)

        min_index = np.argmin(np.array(mins))
        patch_ingredients = []

        min_value = (last_patch[min_index].points[n-1]-z0)/((min_index+1)/20)

        for i in range(5):
            y = 5*(i+1)
            closest_index = np.argmin(np.abs((self.data.array[i]-z0)/((i+1)/20) - min_value))
            stop_index = np.argmin(np.abs(self.data.array[i] - (stop*(100-z0)*y/25 + z0)))

            if stop_index != len(self.data.array[i])-1:
                stop_index += 1
            stop_value = self.data.array[i][stop_index]

            if closest_index > stop_index - n + 1:
                if (stop_index - n +1) < 0:
                    stop_index = n-1
                patch_ingredients.append(SuperPoint(self.data.array[i][stop_index-n+1:stop_index+1]))

            
            else:
                patch_ingredients.append(SuperPoint(self.data.array[i][closest_index:closest_index+n]))

        new_patch = Patch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_relaxed_gap(z0=z0, stop=stop, n = n)

    def S_reverse_relaxed_gap(self, z0= 0, stop = -1, n = 16):

        loops = self.n_patches - 1
        mins = []
        stop_count = 0

        for i in range(5):

            last_patch = self.patches[loops].superpoints

            amin = (last_patch[i].points[0]-z0)/((i+1)/20)
            mins.append(amin)

        min_index = np.argmax(np.array(mins))
        patch_ingredients = []

        min_value = (last_patch[min_index].points[0]-z0)/((min_index+1)/20)

        for i in range(5):
            y= 5*(i+1)
            closest_index = np.argmin(np.abs((self.data.array[i]-z0)/((i+1)/20) - min_value))

            stop_index = np.argmin(np.abs(self.data.array[i] - (((stop*(z0+100)*y)/25+z0))))
            if stop_index != 0:
                stop_index -= 1
            stop_value = self.data.array[i][stop_index]

            if closest_index < stop_index + n - 1:
                if (stop_index + n) > len(self.data.array[i]):
                    stop_index = int(len(self.data.array[i]) - stop_index-1)
                patch_ingredients.append(SuperPoint(self.data.array[i][stop_index:stop_index+n]))
            
            else:
                patch_ingredients.append(SuperPoint(self.data.array[i][closest_index-n+1:closest_index+1]))

        new_patch = Patch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_reverse_relaxed_gap(z0=z0, stop=stop, n = n)

    def S_loop_improved(self, stop = 1, n = 16): #NOT USED ANYMORE
        
        loops = self.n_patches - 1
        mins = []

        for i in range(5):

            last_patch = self.patches[loops].superpoints

            amin = last_patch[i].points[n-1]/((i+1)/20)
            mins.append(amin)

        min_index = np.argmin(np.array(mins))
        patch_ingredients = []

        min_value = last_patch[min_index].points[n-1]/((min_index+1)/20)

        for i in range(5):

            closest_index = np.argmin(np.abs(self.data.array[i]/((i+1)/20) - min_value))

            stop_index = np.argmin(np.abs(self.data.array[i] - stop*(i+1)*20))
            if stop_index != len(self.data.array[i])-1:
                stop_index += 1
            stop_value = self.data.array[i][stop_index]

            if closest_index > stop_index - n + 1:
                patch_ingredients.append(SuperPoint(self.data.array[i, stop_index-n+1:stop_index+1]))
            
            else:
                patch_ingredients.append(SuperPoint(self.data.array[i, closest_index-1:closest_index+n-1]))

        new_patch = Patch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_loop_improved(stop, n = n)

    def S_rloop_improved(self, stop = -1, n = 16): #NOT USED ANYMORE

        loops = self.n_patches - 1
        mins = []

        for i in range(5):

            last_patch = self.patches[loops].superpoints

            amin = last_patch[i].points[0]/((i+1)/20)
            mins.append(amin)

        min_index = np.argmax(np.array(mins))
        patch_ingredients = []

        min_value = last_patch[min_index].points[0]/((min_index+1)/20)

        for i in range(5):

            closest_index = np.argmin(np.abs(self.data.array[i]/((i+1)/20) - min_value))

            stop_index = np.argmin(np.abs(self.data.array[i] - stop*(i+1)*20))
            if stop_index != 0:
                stop_index -= 1
            stop_value = self.data.array[i][stop_index]

            if closest_index < stop_index + n - 1:
                patch_ingredients.append(SuperPoint(self.data.array[i, stop_index:stop_index+n]))
            
            else:
                patch_ingredients.append(SuperPoint(self.data.array[i, closest_index-n+2:closest_index+2]))

        new_patch = Patch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_rloop_improved(stop)

    def S_loop15(self, z0 = 0, stop = 1, n = 16):
        
        loops = self.n_patches - 1
        mins = []
        last_patch = self.patches[loops].superpoints
        for i in range(5):
            y = 5*(i+1)
            amin = (last_patch[i].points[n-1]-z0)/y
            mins.append(amin)

        patch_ingredients = []

        min_index = np.argmin(np.array(mins))
        self.patches[loops].end_layer(min_index+1)
        #min_value = (last_patch[min_index].points[n-1]-z0)/(5*(min_index+1))
        min_value = min(mins)

        term = 0

        for i in range(5):
            y = 5*(i+1)
            closest_index = np.argmin(np.abs((self.data.array[i]-z0)/(y) - min_value))
            #if (self.data.array[i][closest_index]-z0)/(y/100) >= min_value:
            #    closest_index -=1
            stop_index = np.argmin(np.abs(self.data.array[i] - (stop*(100-z0)*y/25 + z0)))
            layer_size = len(self.data.array[i])

            if stop_index != len(self.data.array[i])-1:
                stop_index += 1

            if closest_index + n - 1 > stop_index:
            #if closest_index + n > stop_index:
                term += 1

            #if closest_index > stop_index - n + 1:
            if closest_index + n - 1 > layer_size:
                #if (stop_index - n + 1) < 0:
                #    stop_index = n - 1
                #patch_ingredients.append(SuperPoint(self.data.array[i][stop_index-n+1:stop_index+1]))
                patch_ingredients.append(SuperPoint(self.data.array[i][layer_size-n:]))
            
            else:
                if closest_index == 0:
                    closest_index = 1
                patch_ingredients.append(SuperPoint(self.data.array[i][closest_index-1:closest_index + n - 1]))
                #patch_ingredients.append(SuperPoint(self.data.array[i][closest_index:closest_index + n]))

        new_patch = Patch(self.env, tuple(patch_ingredients))

        #if np.array_equal(new_patch.superpoints, last_patch):
            #return
        if term == 5:
            self.add_patch(new_patch)
            return

        else:
            self.add_patch(new_patch)
            return self.S_loop15(z0, stop, n = n)

    def S_rloop15(self, z0 = 0, stop = -1, n = 16):

        loops = self.n_patches - 1
        mins = []
        last_patch = self.patches[loops].superpoints

        for i in range(5):
            y = 5*(i+1)

            amin = (last_patch[i].points[0]-z0)/(y)
            mins.append(amin)

        patch_ingredients = []

        min_index = np.argmax(np.array(mins))
        self.patches[loops].end_layer(min_index+1)
        #min_value = (last_patch[min_index].points[0]-z0)/((min_index+1)/20)
        min_value = max(mins)

        term = 0

        for i in range(5):
            y= 5*(i+1)
            closest_index = np.argmin(np.abs((self.data.array[i]-z0)/(y) - min_value))
            #if (self.data.array[i][closest_index]-z0)/(y/100) <= min_value:
            #    closest_index +=1
            stop_index = np.argmin(np.abs(self.data.array[i] - (stop*(100+z0)*y/25 + z0)))

            stop_index -= 1

            try:
                if self.data.array[i][closest_index] == self.data.array[i][closest_index+1]:
                    closest_index = closest_index +1
            except:
                pass

            if closest_index - n + 2 <= stop_index:
            #if closest_index - n + 1 < stop_index:
                term += 1

            if stop_index != 0:
                stop_index -= 1

            #if closest_index < stop_index + n:
                #if (stop_index + n) > len(self.data.array[i]):
                #    stop_index = int(len(self.data.array[i]) - n)
            if closest_index + 2 < n:
                #patch_ingredients.append(SuperPoint(self.data.array[i][stop_index:stop_index+n]))
                patch_ingredients.append(SuperPoint(self.data.array[i][:n]))

            else:
                if closest_index == len(self.data.array[i]) - 1:
                    closest_index -=1
                patch_ingredients.append(SuperPoint(self.data.array[i][closest_index-n+2:closest_index+2]))
                #patch_ingredients.append(SuperPoint(self.data.array[i][closest_index-n+1:closest_index+1]))


        new_patch = Patch(self.env, tuple(patch_ingredients))
        if term == 5:
            self.add_patch(new_patch)
            return

        #if np.array_equal(new_patch.superpoints, last_patch):
            #return

        else:
            self.add_patch(new_patch)
            return self.S_rloop15(z0, stop, n = n)         

    def solve(self, clustering:str = "", lining:str = "SolveS", z0=0, n = 16, nlines:int=100, show = True): 
        if show == True:
            lGen = LineGenerator(self.env, z0)
            #self.fitting_lines = lGen.generateGridLines(nlines)
            self.fitting_lines = lGen.generateEvenGrid(nlines)

        if lining =="solveS":
            try:
                for s in z0:
                    self.solveS(z0=s, n = n)
            except:
                self.solveS(z0=z0, n = n)
            return

        elif lining =="solveS_reverse":
            try:
                for s in z0:
                    self.solveS_reverse(z0=s, stop = -1, n = n)
            except:
                self.solveS_reverse(z0=z0, stop = -1, n = n)
            return

        elif lining == "solveS_center1": 
            self.solveS_center1(n = n)
            return 

        elif lining == "solveS_center2":
            try:
                for s in z0:
                    self.solveS_center2(z0=s, n = n)
            except:
                self.solveS_center2(z0=z0, n = n)
            return 

        elif lining == "solveQ":
            try:
                for s in z0:
                    self.solveQ(z0=s, n = n)
            except:
                self.solveQ(z0=z0, n = n)
            return 

        elif lining == "solveS_relaxed_end":
            try:
                for s in z0:
                    self.solveS_relaxed_end(z0=s, n = n)
            except:
                self.solveS_relaxed_end(z0=z0, n = n)
            return 

        elif lining == "solveS_relaxed_gap":
            try:
                for s in z0:
                    self.solveS_relaxed_gap(z0=s, n = n)
            except:
                self.solveS_relaxed_gap(z0=z0, n = n)
            return 

        elif lining == "solveS_relaxed_both":
            try:
                for s in z0:
                    self.solveS_relaxed_both(z0=s, n = n)
            except:
                self.solveS_relaxed_both(z0=z0, n = n)
            return 

        elif lining == "solveQ_relaxed_gap":
            try:
                for s in z0:
                    self.solveQ_relaxed_gap(z0=s, n = n)
            except:
                self.solveQ_relaxed_gap(z0=z0, n = n)
            return

        elif lining == "solveQ_relaxed_end":
            try:
                for s in z0:
                    self.solveQ_relaxed_end(z0=s, n = n)
            except:
                self.solveQ_relaxed_end(z0=z0, n = n)
            return


        elif lining == "solveQ_relaxed_both":
            try:
                for s in z0:
                    self.solveQ_relaxed_both(z0=s, n = n)
            except:
                self.solveQ_relaxed_both(z0=z0, n = n)
            return
        
        if clustering == "LeftRight": 
            self.cluster("LR") 
        elif clustering == "Center": 
            self.cluster("C") 
        else: 
            raise Exception("That is not a valid clustering type. ")
            
        if lining == "LeftRight": 
            self.solveGridLR(nlines) 
        elif lining == "Randomized": 
            self.solveRandomized(nlines) 
        elif lining == "CenterGrid": 
            self.solveCenterGrid(nlines) 
        elif lining == "CenterSpread": 
            self.solveCenterSpread(nlines) 
        else: 
            raise Exception("That is not a valid lining type. ")
             
    def plot(self, data=True, lines=True, patches=True): 
        
        if self.n_patches == 0: 
            raise Exception("You must run the solve method first to generate the patches. ")

        if patches == False: 
            if lines == False and data == False: 
                pass 
            elif lines == False and data == True: 
                self.data.plot("g") 
                plt.show() 
            elif lines == True and data == False: 
                for line in self.fitting_lines: 
                    line.plot('r') 
                plt.show() 
            elif lines == True and data == True: 
                self.data.plot("g") 
                for line in self.fitting_lines: 
                    line.plot('r') 
                plt.show() 
            
        elif patches == True: 
            
            if lines == False and data == False: 
                name = 0
                for patch in self.patches: 
                    patch.plot("b")
                    plt.savefig(f"temp_image_dir/{two_digit(name)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == False and data == True: 
                name = 0
                for patch in self.patches: 
                    self.data.plot("g")
                    patch.plot("b")
                    plt.savefig(f"temp_image_dir/{two_digit(name)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == True and data == False: 
                name = 0
                for patch in self.patches: 
                    for line in self.fitting_lines: 
                        if patch.contains(line): 
                            line.plot("r")
                    patch.plot("b")
                
                    plt.savefig(f"temp_image_dir/{two_digit(name)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == True and data == True: 
                name = 0
                for patch in self.patches: 
                    self.data.plot("g")
                    for line in self.fitting_lines: 
                        if patch.contains(line): 
                            line.plot("r")
                    patch.plot("b")
                    plt.savefig(f"temp_image_dir/{two_digit(name)}.png")
                    plt.clf() 
                    name += 1 
   
                
                        
            image_files = glob.glob("temp_image_dir/*.png")
            
            image_files.sort()

            # index to keep track of current image
            current_image_idx = 0
                    
            while True:
                # read and display the current image
                image = cv2.imread(image_files[current_image_idx])
                cv2.imshow("Slideshow", image)

                # wait for a key press and check its value
                key = cv2.waitKey(0) & 0xFF

                # if the 'e' key is pressed, go to the next image
                if key == ord('e'):
                    current_image_idx += 1
                    # wrap around if we're at the end of the list
                    if current_image_idx == len(image_files):
                        current_image_idx = 0

                # if the 'w' key is pressed, go to the previous image
                elif key == ord('w'):
                    current_image_idx -= 1
                    # wrap around if we're at the start of the list
                    if current_image_idx < 0:
                        current_image_idx = len(image_files) - 1

                # if the 'q' key is pressed, break the loop and end the slideshow
                elif key == ord('q'): 
                    break
                
            for file in image_files: 
                os.remove(file) 
                print(f"Deleted File: {file}")
                
            cv2.destroyAllWindows()
                    

