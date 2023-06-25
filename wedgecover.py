import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from wedgedata import *
import math
import cv2 
import os 
import glob
from cover import *
        
class wedgeSuperPoint(): 
    
    def __init__(self, points):
        
        if np.size(points) != 16:
            if np.size(points) != 32:
                raise Exception("This patch does not have 16 or 32 points in each layer")
        z_list = np.array([points[x].z for x in range(len(points))])
        self.z_values = z_list
        self.points = points 
        self.min = np.min(z_list) 
        self.max = np.max(z_list)
        
    def contains(self, p:float): 
        return self.min <= p <= self.max
    
    def __eq__(self, other): 
        return (self.min, self.max) == (other.min, other.max)
class wedgePatch(): 
    
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

    def contains_p(self, point:float, layer:int): 
        
        sp = self.superpoints[layer] 
        return sp.contains(point)
    
    def __eq__(self, other): 
        if not isinstance(other, wedgePatch): 
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

class wedgeCover(): 
    
    def __init__(self, env:Environment, data:WedgeData): 
        self.n_patches = 0 
        self.patches = [] 
        self.env = env 
        self.data = data 
        self.fitting_lines = [] 
        self.superPoints = [] 
        
    def add_patch(self, curr_patch:wedgePatch): 
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

    def solve(self, lining:str = "SolveS", z0=0, n = 16, nlines:int=100, show = True): 
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

    def S_loop15(self, z0 = 0, stop = 1, n = 16):
        
        loops = self.n_patches - 1
        mins = []

        for i in range(5):
            y = 5*(i+1)
            last_patch = self.patches[loops].superpoints
            row_data = last_patch[i].points
            row_list = np.array([row_data[x].z for x in range(len(row_data))])
            amin = (row_list[n-1]-z0)/(y/100)
            mins.append(amin)

        min_index = np.argmin(np.array(mins))
        patch_ingredients = []
        min_value = (last_patch[min_index].points[n-1].z-z0)/(5*(min_index+1)/100)

        for i in range(5):
            y = 5*(i+1)
            row_data = self.data.array
            row_list = np.array([row_data[i][x].z for x in range(len(row_data[i]))])
            closest_index = np.argmin(np.abs((row_list-z0)/(y/100) - min_value))
            stop_index = np.argmin(np.abs(row_list - (stop*(100-z0)*y/25 + z0)))

            if stop_index != len(row_list)-1:
                stop_index += 1
            stop_value = row_list[stop_index]

            if closest_index > stop_index - n + 1:
                if (stop_index - n + 1) < 0:
                    stop_index = n - 1
                patch_ingredients.append(wedgeSuperPoint(row_data[i][stop_index-n+1:stop_index+1]))
            
            else:
                if closest_index == 0:
                    closest_index = 1
                patch_ingredients.append(wedgeSuperPoint(row_data[i][closest_index-1:closest_index + n - 1]))

        new_patch = wedgePatch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_loop15(z0, stop, n = n)

    def S_rloop15(self, z0 = 0, stop = -1, n = 16):

        loops = self.n_patches - 1
        mins = []
        stop_count = 0

        for i in range(5):
        	y = 5*(i+1)
        	last_patch = self.patches[loops].superpoints
        	row_data = last_patch[i].points
        	row_list = np.array([row_data[x].z for x in range(len(row_data))])
        	amin = (row_list[0]-z0)/(y/100)
        	mins.append(amin)

        min_index = np.argmax(np.array(mins))
        patch_ingredients = []
        min_value = (last_patch[min_index].points[0].z-z0)/((min_index+1)/20)

        for i in range(5):
            y= 5*(i+1)
            row_data = self.data.array
            row_list = np.array([row_data[i][x].z for x in range(len(row_data[i]))])
            closest_index = np.argmin(np.abs((row_list-z0)/((i+1)/20) - min_value))
            stop_index = np.argmin(np.abs(row_list - (((stop*(z0+100)*y)/25+z0))))

            if stop_index != 0:
                stop_index -= 1
            stop_value = row_list[stop_index]

            if closest_index < stop_index + n - 1:
                if (stop_index + n) > len(self.data.array[i]):
                    stop_index = int(len(self.data.array[i]) - stop_index-1)
                patch_ingredients.append(wedgeSuperPoint(row_data[i][stop_index:stop_index+n]))
            else:
                patch_ingredients.append(wedgeSuperPoint(row_data[i][closest_index-n+2:closest_index+2]))


        new_patch = wedgePatch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_rloop15(z0, stop, n = n)         

    def solveS(self, z0 = 0, stop = 1, n = 16):

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            y = 5*(row + 1)
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            start_index = np.argmin(np.abs(row_list - (((-z0-100)*y)/25+z0)))
            if start_index != 0:
                start_index -= 1
            init_patch.append(wedgeSuperPoint(row_data[row][start_index:start_index+n]))

        #add to patch
        self.add_patch(wedgePatch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_loop15(z0=z0, stop=stop, n=n)
        return

    def solveS_reverse(self, z0 = 0, stop = -1, n = 16):

        init_patch = []

        #pick first 16 points in each layer
        for row in range(5):
            y = 5*(row+1)
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            start_index = np.argmin(np.abs((row_list - ((100-z0)*y/25 + z0))))

            if start_index != len(self.data.array[row])-1:
                start_index += 1
            init_patch.append(wedgeSuperPoint(row_data[row][start_index-n+1:start_index+1]))

        #add to patch
        self.add_patch(wedgePatch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_rloop15(z0 = z0, stop = -1, n = n)
        return

    def solveS_center2(self, center = 0, z0 = 0, stop = 'none', n = 16):
        init_patch = []

        #pick center 16 points based on 
        for row in range(5):
            y = 5*(row+1)
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            center_index = np.argmin(np.abs(row_list - ((y*(center-z0)/25)+z0)))
            if (center_index-int(n/2)) < 0:
                center_index = int(n/2)
            elif (center_index+int(n/2)) > len(self.data.array[row]):
                center_index = len(self.data.array[row]) - int(n/2)
            init_patch.append(wedgeSuperPoint(row_data[row][center_index-int(n/2):center_index+int(n/2)]))
        #add to patch
        self.add_patch(wedgePatch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, stop = 0, n = n)
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
                self.S_rloop15(z0 = z0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
            #starts right of center
            if center > 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, n = n)
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
                self.S_rloop15(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            #run main algorithm
            n_patch_start = self.n_patches
            self.S_loop15(z0 = z0, n = n)
            self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
            init_patch.append(wedgeSuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))
        #add to patch
        self.add_patch(wedgePatch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_loop15(z0 = z0, stop = 0, n = n)
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
            self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
            init_patch.append(wedgeSuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))

        self.add_patch(wedgePatch(self.env, tuple(init_patch)))

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
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
                self.S_reverse_relaxed_gap(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            #run main algorithm
            n_patch_start = self.n_patches
            self.S_relaxed_gap(z0 = z0, n = n)
            self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
            init_patch.append(wedgeSuperPoint(self.data.array[row, center_index-int(n/2):center_index+int(n/2)]))

        self.add_patch(wedgePatch(self.env, tuple(init_patch)))

        if stop == 'center':
            # starts left of center
            if center < 0:
                n_patch_start = self.n_patches
                self.S_relaxed_gap(z0 = z0, stop = 0, n = n)
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
                self.add_patch(wedgePatch(self.env, tuple(init_patch)))
                self.S_reverse_relaxed_gap(z0 = z0, stop = 0, n = n)
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:

            n_patch_start = self.n_patches
            self.S_relaxed_gap(z0 = z0, n = n)
            del self.patches[self.n_patches-1]
            self.add_patch(wedgePatch(self.env, tuple(init_patch)))
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
                patch_ingredients.append(wedgeSuperPoint(self.data.array[i][stop_index-n+1:stop_index+1]))

            
            else:
                patch_ingredients.append(wedgeSuperPoint(self.data.array[i][closest_index:closest_index+n]))

        new_patch = wedgePatch(self.env, tuple(patch_ingredients))

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

            if closest_index < stop_index + n:
                if (stop_index + n) > len(self.data.array[i]):
                    stop_index = int(len(self.data.array[i]) - n)
                patch_ingredients.append(wedgeSuperPoint(self.data.array[i][stop_index:stop_index+n]))
            
            else:
                patch_ingredients.append(wedgeSuperPoint(self.data.array[i][closest_index-n+1:closest_index+1]))

        new_patch = wedgePatch(self.env, tuple(patch_ingredients))

        if np.array_equal(new_patch.superpoints, last_patch):
            return

        else:
            self.add_patch(new_patch)
            return self.S_reverse_relaxed_gap(z0=z0, stop=stop, n = n)
             
    def plot(self, data=True, lines=True, patches=True): 
        
        if self.n_patches == 0: 
            raise Exception("You must run the solve method first to generate the patches. ")

        if patches == False: 
            if lines == False and data == False: 
                pass 
            elif lines == False and data == True: 
                self.data.plot(show = False) 
                plt.show() 
            elif lines == True and data == False: 
                for line in self.fitting_lines: 
                    line.plot('r') 
                plt.show() 
            elif lines == True and data == True: 
                self.data.plot(show = False) 
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
                    self.data.plot(show = False) 
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
                    self.data.plot(show = False) 
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
                    

