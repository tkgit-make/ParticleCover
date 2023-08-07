import numpy as np 
import matplotlib.pyplot as plt 
from src.coverers.data_structs import * 
from src.coverers.line import *
import math
import cv2 
import os 
import glob
from time import time 
from src.coverers.parallelogram import *
from src.debug import * 
        
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
        
    def contains(self, p): 
        try:
             p = p.z
        except:
            p = p
        return self.min <= p <= self.max
    
    def __eq__(self, other): 
        return (self.min, self.max) == (other.min, other.max)

class wedgePatch(): 
    
    # Should be hashable (nvm we can't make it hashable) 
    
    def __init__(self, env:Environment, superpoints:tuple, apexZ0): 
        self.env = env 
        self.end_layer = -1
        self.apexZ0 = apexZ0
        
        if len(superpoints) != env.num_layers: 
            raise Exception("The patch layers does not match environment layers. ")
        
        self.superpoints = superpoints
        # first superpoint in array should be the 1st layer 
        
        self.getParallelograms()
        
    def contains(self, line:Line): 
        
        for i in range(len(self.superpoints)): 
            if not self.superpoints[i].contains(line.points[i+1]): 
                return False 
            
        return True

    def contains_p(self, point:float, layer:int): 
        
        sp = self.superpoints[layer] 
        return sp.contains(point)

    def add_end(self, layer):
        self.end_layer = layer
    
    def __eq__(self, other): 
        if not isinstance(other, wedgePatch): 
            return NotImplemented 
        
        for i in range(len(self.superpoints)): 
            if self.superpoints[i] != other.superpoints[i]: 
                return False 
            
        return True 
    
    def straightLineProjector(self, z_top, z_j, j): 
        radii_leverArm = self.env.radii_leverArm[j-1]
        return z_top - (z_top - z_j) * radii_leverArm 

    def getParallelograms(self): 
        
        parallelograms = [] 
        
        # min and max of the top superpoint of a patch
        top_layer_zmin = max(self.superpoints[-1].min, -self.env.top_layer_lim)
        top_layer_zmax = min(self.superpoints[-1].max, self.env.top_layer_lim)
        
        if top_layer_zmin > top_layer_zmax: 
            top_layer_zmin = self.env.top_layer_lim + 1
            top_layer_zmax = top_layer_zmin

        for j, superpoint in enumerate(self.superpoints[:-1], start=1): 
            z_j_min = superpoint.min 
            z_j_max = superpoint.max 
            
            a = self.straightLineProjector(top_layer_zmax, z_j_min, j)
            b = self.straightLineProjector(top_layer_zmax, z_j_max, j)
            
            pSlope = self.env.parallelogramSlopes[j-1]
            
            Parallelogram = parallelogram(j, top_layer_zmin, top_layer_zmax, a, b, pSlope)
            parallelograms.append(Parallelogram)
            
        self.parallelograms = parallelograms
    
    def plot(self, color='g'): 
        heights = self.env.radii
        for i in range(self.env.num_layers): 
            sp = self.superpoints[i] 
            
            max_height = self.env.radii[-1]
            
            plt.plot([sp.min, sp.max], [heights[i], heights[i]], c=color)
        
        plt.xticks(np.arange(-self.env.top_layer_lim, self.env.top_layer_lim, 10))
        plt.yticks(np.arange(0, max_height + 1, self.env.num_layers))
        # plt.show()

class wedgeCover(): 
    
    def __init__(self, env:Environment, data:DataSet): 
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

    def solve(self, lining:str = "makePatches_Projective", apexZ0=0, ppl = 16, nlines:int=100, show = True):
        if show == True:
            fitting_lines = []
            if (type(apexZ0) == int) or (type(apexZ0) == float):
                lGen = LineGenerator(self.env, apexZ0)
                fitting_lines = fitting_lines + lGen.generateEvenGrid(nlines)
            else:
                for s in apexZ0:
                    lGen = LineGenerator(self.env, s)
                    fitting_lines = fitting_lines + lGen.generateEvenGrid(nlines)
            self.fitting_lines = fitting_lines
        
        if lining == 'makePatches_Projective':
            try:
                for s in apexZ0:
                    self.makePatches_Projective(apexZ0=s, ppl = ppl, leftRight = True)
            except:
                self.makePatches_Projective(apexZ0=apexZ0, ppl = ppl, leftRight = True)
            return
        
        elif (lining == 'makePatches_Projective_center') or (lining == 'c'):
            try:
                for s in apexZ0:
                    self.makePatches_Projective_center(apexZ0=s, ppl = ppl)
            except:
                self.makePatches_Projective_center(apexZ0=apexZ0, ppl = ppl)
            return
        
        elif (lining == 'makePatches_Projective_quartile') or (lining == 'q'):
            try:
                for s in apexZ0:
                    self.makePatches_Projective_quartile(apexZ0=s, ppl = ppl)
            except:
                self.makePatches_Projective_quartile(apexZ0=apexZ0, ppl = ppl)
            return
        
        elif lining == 'lr':
            try:
                for s in apexZ0:
                    self.makePatches_Projective(apexZ0=s, ppl = ppl, leftRight = True)
            except:
                self.makePatches_Projective(apexZ0=apexZ0, ppl = ppl, leftRight = True)
            return
        
        elif (lining == 'rl') or (lining == 'makePatches_Projective_Rightleft'):
            try:
                for s in apexZ0:
                    self.makePatches_Projective(apexZ0=s, ppl = ppl, leftRight = False)
            except:
                self.makePatches_Projective(apexZ0=apexZ0, ppl = ppl, leftRight = False)
            return
        else:
            raise("Please choose valid solving method")

    def makePatches_Projective_Loop(self, apexZ0 = 0, stop = 1, ppl = 16, leftRight = True):
            """Loop for creating patches left to right or right to left depending on argument leftRight

            Args:
                apexZ0 (num, optional): Places to generate patch. Defaults to 0.
                stop (num, optional): stopping location, normalized to 1m. Defaults to 1.
                ppl (int, optional): points per patch per layer. Defaults to 16.
                leftRight(Bool): If set to true, make patches from left to right, if false, then make from right to left

            Returns:
                function: reruns loop if it hasn't reached the end of the dataset
            """
            #count how many times this has been run
            loops = self.n_patches - 1
            #reads last patch made
            last_patch = self.patches[loops].superpoints
            #create list for points closest to starting line and patch ingredients
            lambdaZ_list = []
            patch_ingredients = []
            #creates count for terminating patch making. loop stops when all layers are beyond line from (z0, 0) to (100, 25)
            term = 0

            #loops through layers
            for i in range(self.env.num_layers):
                y = self.env.radii[i]
                #create compatible arrays from data structure
                row_data = last_patch[i].points
                row_list = np.array([row_data[x].z for x in range(len(row_data))])
                #rescales point for layer and add to mins list
                if leftRight == True:
                    lambdaZ = (row_list[ppl-1]-apexZ0)/y 
                else:
                    lambdaZ = (row_list[0]-apexZ0)/(y)
                lambdaZ_list.append(lambdaZ)

            #find which layer of the next n points from last patch stops first and find rescaled value of that point
            if leftRight == True:
                end_index = np.argmin(lambdaZ_list)
                min_lambdaZ = min(lambdaZ_list)
            else:
                end_index = np.argmax(lambdaZ_list)
                min_lambdaZ = max(lambdaZ_list)

            self.patches[loops].add_end(end_index+1)

            r_max = self.env.radii[-1]
            z_max = self.env.top_layer_lim

            #row_data[layer] gives spacepoints in layer
            row_data = self.data.array
            #loops through layers again
            for i in range(self.env.num_layers):
                y = self.env.radii[i]
                row_list = np.array([row_data[i][x].z for x in range(len(row_data[i]))])
                #finds point closest to line from (z0, 0) to leftmost rescaled point
                closest_index = np.argmin(np.abs((row_list-apexZ0)/(y) - min_lambdaZ))
                #find where the stopping index is based on the line from (z0, 0) to z-edge of outermost layer
                if leftRight == True:
                    stop_index = np.argmin(np.abs(row_list - (stop*(z_max-apexZ0)*y/r_max + apexZ0)))

                    #add one to stop index in case it is left of the line from (z0, 0) to (100*stop, 25)
                    #this makes sure there is full coverage
                    if stop_index != len(row_list)-1:
                        stop_index += 1

                    #checks to see if patch will go past stop index, if so, add one to term variable
                    if closest_index + ppl - 1 > stop_index:
                        term += 1


                    #if there is not enough points left, pick last n points
                    if closest_index + ppl - 1 > len(row_list):
                        patch_ingredients.append(wedgeSuperPoint(row_data[i][len(row_list)-ppl:]))
                    
                    #if there are enough points left, pick point closest to slope and next n-1 points
                    else:
                        #makes sure there won't be an error of negative indices
                        if closest_index == 0:
                            closest_index = 1
                        #closest_index - 1 insures point is to left of line ie ensuring patches overlap
                        patch_ingredients.append(wedgeSuperPoint(row_data[i][closest_index-1:closest_index + ppl - 1]))
                    #print(term)
                else:
                    stop_index = np.argmin(np.abs(row_list - (stop*(z_max+apexZ0)*y/r_max+apexZ0)))
                    #for the extremely specific condition where two z's are equal and it is the edgepoint
                    try:
                        if row_list[closest_index] == row_list[closest_index+1]:
                            closest_index = closest_index + 1 
                    except:
                        pass

                    if stop_index != 0: 
                        stop_index -= 1
                    #checks to see if patch will go past stop index, if so, add one to term variable
                    if closest_index - ppl + 2 <= stop_index:
                        term += 1

                    #if there aren't enough points left, pick leftmost n points
                    if closest_index + 2 < ppl:
                        patch_ingredients.append(wedgeSuperPoint(row_data[i][:ppl]))

                    #if there are enough points left, pick point closest to slope and n-1 points to the left
                    else:
                        #makes sure there won't be an error of indices beyond length of list
                        if closest_index == len(row_list) - 1:
                            closest_index -=1
                        #closest_index + 2 ensures point is to right of line ie ensures patches overlap
                        patch_ingredients.append(wedgeSuperPoint(row_data[i][closest_index - ppl + 2:closest_index + 2]))
                    #print(term)
            #add superpoints to patch
            new_patch = wedgePatch(self.env, tuple(patch_ingredients), apexZ0=apexZ0)
            #add patch to cover
            self.add_patch(new_patch)
            
            #if all layers have points beyond stop index, stop
            if term == 5:
                return
            #if new patches are still being created, repeat loop
            else:
                return self.makePatches_Projective_Loop(apexZ0, stop, ppl = ppl, leftRight = leftRight)
            
    def makePatches_Projective(self, apexZ0 = 0, stop = 1, ppl = 16, leftRight = True):
            """Creates patches left to right or right to left depending on argument leftRight

            Args:
                apexZ0 (num, optional): Places to generate patch. Defaults to 0.
                stop (num, optional): stopping location, normalized to 1m. Defaults to 1.
                n (int, optional): points per patch per layer. Defaults to 16.
                leftRight(Bool): If set to true, make patches from left to right, if false, then make from right to left

            Returns:
                function: runs loop to make patches
            """
            #create list for inital patch
            if (leftRight == False) & (stop == 1):
                stop = -1
            init_patch = []

            #row_data[layer] contains spacepoints for each layer
            row_data = self.data.array
            #loops through each layer and picks n points closest to (z0, 0) and (-100, 25)
            for row in range(self.env.num_layers):
                y = self.env.radii[row]
                #create compatible arrays from data structure
                row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
                #picks picks n points closest to line from (z0, 0) to (-100, 25) (top left point)
                r_max = self.env.radii[-1]
                z_max = self.env.top_layer_lim

                if leftRight == True:
                    start_index = np.argmin(np.abs(row_list - (((-apexZ0-z_max)*y)/r_max+apexZ0)))
                    #subtract one from stop index in case it is right of the line from (z0, 0) to (-100, 25)
                    if start_index != 0:
                        start_index -= 1
                    #add superpoint to patch
                    init_patch.append(wedgeSuperPoint(row_data[row][start_index:start_index+ppl]))
                else:
                    start_index = np.argmin(np.abs((row_list - ((z_max-apexZ0)*y/r_max + apexZ0))))
                    #add one to stop index in case it is left of the line from (z0, 0) to (100, 25)
                    if start_index != len(row_list)-1:
                        start_index += 1
                    #add superpoint
                    init_patch.append(wedgeSuperPoint(row_data[row][start_index-ppl+1:start_index+1]))

                

            #add patch to cover
            self.add_patch(wedgePatch(self.env, tuple(init_patch), apexZ0=apexZ0))

            #run main algorithm
            self.makePatches_Projective_Loop(apexZ0=apexZ0, stop=stop, ppl = ppl, leftRight = leftRight)
            return

    def makePatches_Projective_center(self, center = 0, apexZ0 = 0, stop = 'none', ppl = 16):
        """generate patches starting from center or specified value

        Args:
            center (num, optional): picks where patch making starts from -100 to 100. Defaults to 0.
            apexZ0 (num, optional): Places to generate patch. Defaults to 0.
            stop (str, optional): 'none' or 'center', if center, patch making stops at z = 0. Defaults to 'none'.
            n (int, optional): points per layer per patch. Defaults to 16.
        """
        #create list for inital patch
        init_patch = []
        r_max = self.env.radii[-1]
        #loops through layers and picks picks 16 points closest to (z0, 0) and (0, center) 
        for row in range(self.env.num_layers):
            y = self.env.radii[row]
            #create compatible array from data structure
            row_data = self.data.array
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            #picks n/2 points left and right of point closest to line from (0, 0) to (center, 25)
            center_index = np.argmin(np.abs(row_list - ((y*(center-apexZ0)/r_max)+apexZ0)))
            #conditionals make sure no negative indices indices past length of array
            if (center_index-int(ppl/2)) < 0:
                center_index = int(ppl/2)
                #init_patch.append(wedgeSuperPoint(row_data[row][center_index-int(n/2):center_index+int(n/2)])) #DONT COMMENT IN BUT IS AN ALTERNATIVE
            elif (center_index+int(ppl/2)) > len(self.data.array[row]):
                center_index = len(self.data.array[row]) - int(ppl/2)
                #init_patch.append(wedgeSuperPoint(row_data[row][len(self.data.array-16:len(self.data.array)])) #DONT COMMENT IN BUT IS AN ALTERNATIVE
            if self.data.array[row][center_index].z >= 0:
                init_patch.append(wedgeSuperPoint(row_data[row][center_index-int(ppl/2):center_index+int(ppl/2)]))
            else:
                init_patch.append(wedgeSuperPoint(row_data[row][center_index-int(ppl/2)+1:center_index+int(ppl/2)+1]))            
            
            #init_patch.append(wedgeSuperPoint(row_data[row][center_index-int(ppl/2):center_index+int(ppl/2)]))
        #add initial patch
        self.add_patch(wedgePatch(self.env, tuple(init_patch), apexZ0=apexZ0))

        #for solveQ loops when it needs to stop at line from (0, 0) to (center, 25)
        if stop == 'center':
            # starts left of center
            if center < 0:
                #create index for deleting a patch
                n_patch_start = self.n_patches
                #generates patches right of starting, stopping at center
                self.makePatches_Projective_Loop(apexZ0, stop = 0, ppl = ppl, leftRight = True)
                #add initial patch again
                self.add_patch(wedgePatch(self.env, tuple(init_patch), apexZ0 = apexZ0))
                #generate point left of starting
                self.makePatches_Projective_Loop(apexZ0, stop = -1, ppl = ppl, leftRight = False)
                #delete one of the inital patches so no duplices
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
            #starts right of center
            if center > 0:
                #create index for deleting a patch
                n_patch_start = self.n_patches
                #generates patches right of starting
                self.makePatches_Projective_Loop(apexZ0, ppl = ppl, leftRight = True)
                #add initial patch again
                self.add_patch(wedgePatch(self.env, tuple(init_patch), apexZ0=apexZ0))
                #generates patches left of starting, stopping at center
                self.makePatches_Projective_Loop(apexZ0, stop = 0, ppl = ppl, leftRight = False)
                #delete one of the initial patches so no duplicates
                del self.patches[n_patch_start-1]
                self.n_patches = self.n_patches - 1
                return
        else:
            #create index for deleting a patch
            n_patch_start = self.n_patches
            #generates patches right of starting
            self.makePatches_Projective_Loop(apexZ0, ppl = ppl, leftRight = True)
            #add initial patch again
            self.add_patch(wedgePatch(self.env, tuple(init_patch), apexZ0=apexZ0))
            #generates patches left of starting
            self.makePatches_Projective_Loop(apexZ0, ppl = ppl, stop = -1, leftRight = False)
            #delete one of the initial patches so no duplicates
            del self.patches[n_patch_start-1]
            self.n_patches = self.n_patches - 1
            return          
        
    def makePatches_Projective_quartile(self, apexZ0 = 0, ppl = 16):
        """solves starting at Q1 and Q3 (-50 and 50)

        Args:
            apexZ0 (num, optional): Places to generate patch. Defaults to 0.
            n (int, optional): points per layer per patch. Defaults to 16.
        """
        #solves center starting at -50 and 50, ending at center
        quartile_value = self.env.top_layer_lim / 2
        self.makePatches_Projective_center(center = -quartile_value, stop = 'center', apexZ0 = apexZ0, ppl = ppl) 
        self.makePatches_Projective_center(center = quartile_value, stop = 'center', apexZ0 = apexZ0, ppl = ppl)
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
