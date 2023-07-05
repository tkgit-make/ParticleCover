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
        self.points = (1 / slope) * np.array([0] + self.env.radii) + start 

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
        
        if not(-env.bottom_layer_lim <= start <= env.bottom_layer_lim): 
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

        ycoor = 25
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
    
    def __init__(self, points:list):
        
        if len(points) not in [16, 32]:
            raise Exception("This superpoint does not have 16 or 32 points. ")
        
        layers = set([p.layer_num for p in points])
        if len(layers) != 1: 
            raise Exception("These points are not all on the same layer. ")
        
        self.points = points 
        self.min = min(points, key=attrgetter("z")) 
        self.max = max(points, key=attrgetter("z")) 
        self.layer = layers[0]
        
    def contains(self, p:Point): 
        return self.min <= p.z <= self.max 
    
    def __eq__(self, other): 
        return (self.min, self.max) == (other.min, other.max) and self.layers == other.layers

class Patch(): 
    
    def __init__(self, env:Environment, superpoints:tuple): 
        self.env = env 
        
        if len(superpoints) != env.layers: 
            raise Exception("The patch layers does not match environment layers. ")
        
        self.superpoints = superpoints
        
    def contains(self, line:Line): 
        
        for i in range(len(self.superpoints)): 
            if not self.superpoints[i].min <= line.points[i] <= self.superpoints[i].max: 
                return False 
            
        return True

    def contains_p(self, p:Point): 
        relevant_superpoint = self.superpoints[p.layer_num - 1]
        return relevant_superpoint.contains(p)
    
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

class Cover(): 
    
    def __init__(self, data:DataSet): 
        self.n_patches = 0 
        self.patches = [] 
        self.env = data.env 
        self.data = data 
        self.fitting_lines = [] 
        self.superPoints = [] 
        
        
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
        """Loop for creating patches left to right

        Args:
            z0 (num, optional): Places to generate patch. Defaults to 0.
            stop (num, optional): stopping location, normalized to 1m. Defaults to 1.
            n (int, optional): points per layer per patch. Defaults to 16.

        Returns:
            function: reruns loop if it hasn't reached the end of the dataset
        """

        #count how many times this has been run
        loops = self.n_patches - 1
        #create list for points closest to starting line and patch ingredients
        mins = []
        patch_ingredients = []
        #creates count for terminating patch making. loop stops when all layers are beyond line from (z0, 0) to (100, 25)
        term = 0

        #loops through layers
        for i in range(5):
            y = 5*(i+1)
            #reads last patch made
            last_patch = self.patches[loops].superpoints
            #create compatible arrays from data structure
            row_data = last_patch[i].points
            row_list = np.array([row_data[x].z for x in range(len(row_data))])
            #rescales point for layer and add to mins list
            amin = (row_list[n-1]-z0)/(y/100)
            mins.append(amin)

        #find which layer of the next n points from last patch stops first and find rescaled value of that point
        min_index = np.argmin(np.array(mins))
        min_value = (last_patch[min_index].points[n-1].z-z0)/(5*(min_index+1)/100)
        
        #loops through layers again
        for i in range(5):
            y = 5*(i+1)
            #create compatible array from data structure
            row_data = self.data.array
            row_list = np.array([row_data[i][x].z for x in range(len(row_data[i]))])
            #finds point closest to line from (z0, 0) to leftmost rescaled point
            closest_index = np.argmin(np.abs((row_list-z0)/(y/100) - min_value))
            #find where the stopping index is based on the line from (z0, 0) to (100*stop, 25)
            stop_index = np.argmin(np.abs(row_list - (stop*(100-z0)*y/25 + z0)))

            #add one to stop index in case it is left of the line from (z0, 0) to (100*stop, 25)
            #this makes sure there is full coverage
            ########if stop_index != len(row_list)-1: (old conditional that I'm not sure if we need)
            stop_index += 1

            #checks to see if patch will go past stop index, if so, add one to term variable
            if closest_index + n -1 > stop_index:
                term += 1


            #if there is not enough points left, pick last n points
            if closest_index + n - 1 > len(row_list):
                patch_ingredients.append(SuperPoint(row_data[i][len(row_list)-n:]))
            
            #if there are enough points left, pick point closest to slope and next n-1 points
            else:
                #makes sure there won't be an error of negative indices
                if closest_index == 0:
                    closest_index = 1
                #closest_index - 1 insures point is to left of line ie ensuring patches overlap
                patch_ingredients.append(SuperPoint(row_data[i][closest_index-1:closest_index + n - 1]))

        #creates new patch
        new_patch = Patch(self.env, tuple(patch_ingredients))

        #if all layers have points beyond stop index, add patch and stop
        if term == 5:
            self.add_patch(new_patch)
            return

        #if new patches are still being created, add patch to cover instance and repeat loop
        else:
            self.add_patch(new_patch)
            return self.S_loop15(z0, stop, n = n)

    def S_rloop15(self, z0 = 0, stop = -1, n = 16):
        """Loop for creating patches right to left

        Args:
            z0 (num, optional): Places to generate patch. Defaults to 0.
            stop (num, optional): stopping location, normalized to 1m. Defaults to -1.
            n (int, optional): points per layer per patch. Defaults to 16.

        Returns:
            function: reruns loop if it hasn't reached the end of the dataset
        """

        #count how many times this has been run
        loops = self.n_patches - 1
        #create list for points closest to starting line and patch ingredients
        mins = []
        patch_ingredients = []
        #creates count for terminating patch making. loop stops when all layers are beyond line from (z0, 0) to (-100, 25)
        term = 0

        #loops through layers
        for i in range(5):
            y = 5*(i+1)
            #reads last patch made
            last_patch = self.patches[loops].superpoints
            #create compatible arrays from data structure
            row_data = last_patch[i].points
            row_list = np.array([row_data[x].z for x in range(len(row_data))])
            #rescales point for layer and add to mins list
            amin = (row_list[0]-z0)/(y/100)
            mins.append(amin)

        #find which layer of the next n points from last patch stops first and find rescaled value of that point
        min_index = np.argmax(np.array(mins))
        min_value = (last_patch[min_index].points[0].z-z0)/((min_index+1)/20)

        #loops through layers again
        for i in range(5):
            y= 5*(i+1)
            #create compatible array from data structure
            row_data = self.data.array
            row_list = np.array([row_data[i][x].z for x in range(len(row_data[i]))])
            #finds point closest to line from (z0, 0) to leftmost rescaled point
            closest_index = np.argmin(np.abs((row_list-z0)/((i+1)/20) - min_value))
            #find where the stopping index is based on the line from (z0, 0) to (-100*stop, 25)
            stop_index = np.argmin(np.abs(row_list - (((stop*(z0+100)*y)/25+z0))))

            #subtract one from stop index in case it is right of the line from (z0, 0) to (-100*stop, 25)
            #this makes sure there is full coverage
            ########if stop_index != 0: (old conditional that I'm not sure if we need)
            stop_index -= 1

            #checks to see if patch will go past stop index, if so, add one to term variable
            if closest_index - n + 1 < stop_index:
                term += 1

            #if there aren't enough points left, pick leftmost n points
            if closest_index + 2 < n:
                patch_ingredients.append(SuperPoint(row_data[i][:n]))

            #if there are enough points left, pick point closest to slope and n-1 points to the left
            else:
                #makes sure there won't be an error of indices beyond length of list
                if closest_index == len(row_list) - 1:
                    closest_index -=1
                #closest_index + 2 insures point is to right of line ie ensures patches overlap
                patch_ingredients.append(SuperPoint(row_data[i][closest_index - n + 2:closest_index + 2]))

        #creates new patch
        new_patch = Patch(self.env, tuple(patch_ingredients))

        #if all layers have points beyond stop index, add patch and stop
        if term == 5:
            self.add_patch(new_patch)
            return

        #if new patches are still being created, add patch to cover instance and repeat loop
        else:
            self.add_patch(new_patch)
            return self.S_rloop15(z0, stop, n = n)         

    def solveS(self, z0 = 0, stop = 1, n = 16):
        """Creates patches left to right

        Args:
            z0 (num, optional): Places to generate patch. Defaults to 0.
            stop (num, optional): stopping location, normalized to 1m. Defaults to 1.
            n (int, optional): points per layer per patch. Defaults to 16.

        Returns:
            function: runs loop to make patches
        """
        #create list for inital patch
        init_patch = []

        #loops through each layer and picks n points closest to (z0, 0) and (-100, 25)
        for row in range(5):
            y = 5*(row + 1)
            #create compatible arrays from data structure
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            #picks picks n points closest to (z0, 0) and (-100, 25) 
            start_index = np.argmin(np.abs(row_list - (((-z0-100)*y)/25+z0)))
            #subtract one from stop index in case it is right of the line from (z0, 0) to (-100, 25)
            if start_index != 0:
                start_index -= 1
            #add superpoint
            init_patch.append(SuperPoint(row_data[row][start_index:start_index+n]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_loop15(z0=z0, stop=stop, n=n)
        return

    def solveS_reverse(self, z0 = 0, stop = -1, n = 16):
        """Creates patches right to left

        Args:
            z0 (num, optional): Places to generate patch. Defaults to 0.
            stop (num, optional): stopping location, normalized to 1m. Defaults to -1.
            n (int, optional): points per layer per patch. Defaults to 16.

        Returns:
            function: runs loop to make patches
        """
        #create list for inital patch
        init_patch = []

        #loops through layers and picks picks n points closest to (z0, 0) and (100, 25) 
        for row in range(5):
            y = 5*(row+1)
            #create compatible array from data structure
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            #picks picks n points closest to (z0, 0) and (100, 25) 
            start_index = np.argmin(np.abs((row_list - ((100-z0)*y/25 + z0))))
            #add one to stop index in case it is left of the line from (z0, 0) to (100, 25)
            if start_index != len(self.data.array[row])-1:
                start_index += 1
            #add superpoint
            init_patch.append(SuperPoint(row_data[row][start_index-n+1:start_index+1]))

        #add to patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #run main algorithm
        self.S_rloop15(z0 = z0, stop = -1, n = n)
        return

    def solveS_center2(self, center = 0, z0 = 0, stop = 'none', n = 16):
        """generate patches starting from center or specified value

        Args:
            center (num, optional): picks where patch making starts from -100 to 100. Defaults to 0.
            z0 (num, optional): Places to generate patch. Defaults to 0.
            stop (str, optional): 'none' or 'center', if center, patch making stops at z = 0. Defaults to 'none'.
            n (int, optional): points per layer per patch. Defaults to 16.
        """
        #create list for inital patch
        init_patch = []

        #loops through layers and picks picks 16 points closest to (z0, 0) and (0, center) 
        for row in range(5):
            y = 5*(row+1)
            #create compatible array from data structure
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            #picks n/2 points left and right of point closest to line from (0, 0) to (center, 25)
            center_index = np.argmin(np.abs(row_list - ((y*(center-z0)/25)+z0)))
            #conditionals make sure no negative indices indices past length of array
            if (center_index-int(n/2)) < 0:
                center_index = int(n/2)
            elif (center_index+int(n/2)) > len(self.data.array[row]):
                center_index = len(self.data.array[row]) - int(n/2)
            init_patch.append(SuperPoint(row_data[row][center_index-int(n/2):center_index+int(n/2)]))
        #add initial patch
        self.add_patch(Patch(self.env, tuple(init_patch)))

        #for solveQ loops when it needs to stop at line from (0, 0) to (center, 25)
        if stop == 'center':
            # starts left of center
            if center < 0:
                #create index for deleting a patch
                n_patch_start = self.n_patches
                #generates patches right of starting, stopping at center
                self.S_loop15(z0 = z0, stop = 0, n = n)
                #add initial patch again
                self.add_patch(Patch(self.env, tuple(init_patch)))
                #generate point left of starting
                self.S_rloop15(z0 = z0, n = n)
                #delete one of the inital patches so no duplices
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
            #starts right of center
            if center > 0:
                #create index for deleting a patch
                n_patch_start = self.n_patches
                #generates patches right of starting
                self.S_loop15(z0 = z0, n = n)
                #add initial patch again
                self.add_patch(Patch(self.env, tuple(init_patch)))
                #generates patches left of starting, stopping at center
                self.S_rloop15(z0 = z0, stop = 0, n = n)
                #delete one of the initial patches so no duplicates
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:
            #create index for deleting a patch
            n_patch_start = self.n_patches
            #generates patches right of starting
            self.S_loop15(z0 = z0, n = n)
            #add initial patch again
            self.add_patch(Patch(self.env, tuple(init_patch)))
            #generates patches left of starting
            self.S_rloop15(z0 = z0, n = n)
            #delete one of the initial patches so no duplicates
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 1
            return
        
    def solveQ(self, z0 = 0, n = 16):
        """solves starting at Q1 and Q3 (-50 and 50)

        Args:
            z0 (num, optional): Places to generate patch. Defaults to 0.
            n (int, optional): points per layer per patch. Defaults to 16.
        """
        #solves center starting at -50 and 50, ending at center
        self.solveS_center2(center = -50, stop = 'center', z0 = z0, n = n) 
        self.solveS_center2(center = 50, stop = 'center', z0 = z0, n = n)
        return
             
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
                    plt.savefig(f"temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == False and data == True: 
                name = 0
                for patch in self.patches: 
                    self.data.plot(show = False) 
                    patch.plot("b")
                    plt.savefig(f"temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == True and data == False: 
                name = 0
                for patch in self.patches: 
                    for line in self.fitting_lines: 
                        if patch.contains(line): 
                            line.plot("r")
                    patch.plot("b")
                
                    plt.savefig(f"temp_image_dir/{str(name).zfill(2)}.png")
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
                    plt.savefig(f"temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
   
                
                        
            image_files = glob.glob("temp_image_dir/*.png")
            
            image_files.sort()

            # index to keep track of current image
            current_image_idx = 0
                    
            while True:
                # read and display the current image
                image = imread(image_files[current_image_idx])
                imshow("Slideshow", image)

                # wait for a key press and check its value
                key = waitKey(0) & 0xFF

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
                
            destroyAllWindows()