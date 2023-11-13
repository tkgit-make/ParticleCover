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
            if (np.size(points) != 32) and (np.size(points) != 31):
                raise Exception("This patch does not have 16 or 32/31 points in each layer")
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
        self.left_end_layer = -1
        self.right_end_layer = -1
        self.left_end_lambdaZ = None
        self.right_end_lambdaZ = None
        self.apexZ0 = apexZ0

        self.shadow_fromTopToInnermost_topL_jL = None
        self.shadow_fromTopToInnermost_topL_jR = None
        self.shadow_fromTopToInnermost_topR_jL = None
        self.shadow_fromTopToInnermost_topR_jR = None

        if len(superpoints) != env.num_layers: 
            raise Exception("The patch layers does not match environment layers. ")
        
        self.superpoints = superpoints
        # first superpoint in array should be the 1st layer 
        
        self.getParallelograms()
        self.getParallelograms_v1()
        self.get_acceptanceCorners()
        self.get_end_layer()
        
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
        #return the z0 that is projected by z_top and z_j of layer j
        return z_top - (z_top - z_j) * radii_leverArm 

    def straightLineProjectorFromLayer1(self, z_1, z_j, j):
        radii_leverArm = self.env.radii[j-1] / (self.env.radii[j-1] - self.env.radii[0])
        #return z_0 projected by z_j and z_1
        return z_j - (z_j - z_1) * radii_leverArm
    
    def straightLineProjectorFromLayerIJtoK(self, z_i, z_j, i, j, k):
        # i,j,k = 0 implies beam axis
        radius_i = 0
        radius_j = 0
        radius_k = 0
        if (i == 0):
            radius_i = 0
        else:
            radius_i = self.env.radii[i-1]
        if (j == 0):
            radius_j = 0
        else:
            radius_j = self.env.radii[j-1]            
        if (k == 0):
            radius_k = 0
        else:
            radius_k = self.env.radii[k-1]            
        radii_leverArm = (radius_k - radius_i) / (radius_j - radius_i)
        #return z_k projected by z_j and z_i
        return z_i + (z_j - z_i) * radii_leverArm
    
    def getParallelograms(self): 

        parallelograms = [] 
        
        # min and max of the top superpoint of a patch
        z1_min = max(self.superpoints[0].min, -self.env.trapezoid_edges[0])
        z1_max = min(self.superpoints[0].max, self.env.trapezoid_edges[0])
        
        if z1_min > z1_max: 
            z1_min = self.env.trapezoid_edges[0] + 1
            z1_max = z1_min

        for j, superpoint in enumerate(self.superpoints[1:], start=2): 
            z_j_min = superpoint.min 
            z_j_max = superpoint.max 
            
            a = self.straightLineProjectorFromLayerIJtoK(z1_min, z_j_max, 1, j, self.env.num_layers)
            b = self.straightLineProjectorFromLayerIJtoK(z1_max, z_j_max, 1, j, self.env.num_layers)
            c = self.straightLineProjectorFromLayerIJtoK(z1_min, z_j_min, 1, j, self.env.num_layers)
            d = self.straightLineProjectorFromLayerIJtoK(z1_max, z_j_min, 1, j, self.env.num_layers)
            if j != self.env.num_layers:
                pSlope = self.env.parallelogramSlopes[j-1]
            else:
                pSlope = np.Inf
            Parallelogram = parallelogram(j, z1_min, z1_max, a, b, c, d, pSlope)
            parallelograms.append(Parallelogram)
        
        self.parallelograms = parallelograms

    def getShadows(self, zTopMin, zTopMax):
        
        # AVK add shadows from topLayer coordinate to innermost layer
        zTop_min = max(zTopMin, -self.env.trapezoid_edges[self.env.num_layers-1])
        zTop_max = min(zTopMax, self.env.trapezoid_edges[self.env.num_layers-1])
        topL_jL = []
        topL_jR = []
        topR_jL = []
        topR_jR = []
        for j, superpoint in enumerate(self.superpoints[:-1], start=1):
            z_j_min = superpoint.min 
            z_j_max = superpoint.max 
            topL_jL.append(self.straightLineProjectorFromLayerIJtoK(zTop_min, z_j_min, self.env.num_layers, j, 1))
            topL_jR.append(self.straightLineProjectorFromLayerIJtoK(zTop_min, z_j_max, self.env.num_layers, j, 1))
            topR_jL.append(self.straightLineProjectorFromLayerIJtoK(zTop_max, z_j_min, self.env.num_layers, j, 1))
            topR_jR.append(self.straightLineProjectorFromLayerIJtoK(zTop_max, z_j_max, self.env.num_layers, j, 1))
        self.shadow_fromTopToInnermost_topL_jL = max(topL_jL)
        self.shadow_fromTopToInnermost_topL_jR = min(topL_jR)
        self.shadow_fromTopToInnermost_topR_jL = max(topR_jL)
        self.shadow_fromTopToInnermost_topR_jR = min(topR_jR)
        
    def getParallelograms_v1(self): 

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
            Parallelogram = parallelogram_v1(j, top_layer_zmin, top_layer_zmax, a, b, pSlope)
            parallelograms.append(Parallelogram)
        
        self.parallelograms_v1 = parallelograms

    def get_acceptanceCorners(self):

        self.squareAcceptance = True # Ashutosh
        self.flatTop = True # Ashutosh
        self.flatBottom = True # Ashutosh
        self.triangleAcceptance = False # Ashutosh

        #corner list in z_top
        a_corner_list = [pgram.shadow_bottomL_jR for pgram in self.parallelograms]
        b_corner_list = [pgram.shadow_bottomR_jR for pgram in self.parallelograms]
        c_corner_list = [pgram.shadow_bottomL_jL for pgram in self.parallelograms]
        d_corner_list = [pgram.shadow_bottomR_jL for pgram in self.parallelograms]

        self.a_corner = (self.parallelograms[0].z1_min, min(a_corner_list))
        self.b_corner = (self.parallelograms[0].z1_max, min(b_corner_list))
        self.c_corner = (self.parallelograms[0].z1_min, max(c_corner_list))
        self.d_corner = (self.parallelograms[0].z1_max, max(d_corner_list))

        # is layer5 the most restrictive acceptance? 
        if min(a_corner_list) != a_corner_list[self.env.num_layers-2]:
            self.squareAcceptance = False
            self.flatTop = False
        if min(b_corner_list) != b_corner_list[self.env.num_layers-2]:
            self.squareAcceptance = False
            self.flatTop = False
        if max(c_corner_list) != c_corner_list[self.env.num_layers-2]:
            self.squareAcceptance = False
            self.flatBottom = False
        if max(d_corner_list) != d_corner_list[self.env.num_layers-2]:
            self.squareAcceptance = False
            self.flatBottom = False

        # is the acceptance a triangle shape?
        if (self.c_corner[1] > self.a_corner[1]):
            self.triangleAcceptance = True
            self.c_corner = (self.c_corner[0],self.b_corner[1])
            self.a_corner = (self.a_corner[0],self.b_corner[1])
        if (self.b_corner[1] < self.d_corner[1]):
            self.triangleAcceptance = True
            self.b_corner = (self.b_corner[0],self.c_corner[1])
            self.d_corner = (self.d_corner[0],self.c_corner[1])

    def get_acceptanceCorners_v0(self):

        self.squareAcceptance = True # Ashutosh
        
        a_corner_list = [pgram.shadow_topR_jL for pgram in self.parallelograms]
        b_corner_list = [pgram.shadow_topR_jR for pgram in self.parallelograms]
        c_corner_list = [pgram.shadow_topL_jL for pgram in self.parallelograms]
        d_corner_list = [pgram.shadow_topL_jR for pgram in self.parallelograms]
        #print ("a_corner_list: ", a_corner_list)
        #print ("b_corner_list: ", b_corner_list)
        #print ("c_corner_list: ", c_corner_list)
        #print ("d_corner_list: ", d_corner_list)
        
        self.a_corner = (self.parallelograms[-1].top_layer_zmax, max(a_corner_list))
        self.b_corner = (self.parallelograms[-1].top_layer_zmax, min(b_corner_list))
        self.c_corner = (self.parallelograms[-1].top_layer_zmin, max(c_corner_list))
        self.d_corner = (self.parallelograms[-1].top_layer_zmin, min(d_corner_list))

        # is layer1 the most restrictive acceptance? 
        if np.argmax(a_corner_list) != 0:
            self.squareAcceptance = False
        if np.argmin(b_corner_list) != 0:
            self.squareAcceptance = False
        if np.argmax(c_corner_list) != 0:
            self.squareAcceptance = False
        if np.argmin(d_corner_list) != 0:
            self.squareAcceptance = False

        if self.b_corner <= self.a_corner:

            a_corner_list = []
            b_corner_list = []

            for layer in np.arange(self.env.num_layers-2)+1:
                #compute a_corners by intersecting b_line of layer 2, 3, 4 with a_line of layer 1
                a_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmax, 
                    self.parallelograms[layer].shadow_topR_jR),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmax, self.parallelograms[0].shadow_topR_jL),
                    self.parallelograms[0].pSlope))

                #compute b_corners by intersecting b_line of layer 2, 3, 4 with b_line of layer 1
                b_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmax, 
                    self.parallelograms[layer].shadow_topR_jR),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmax, self.parallelograms[0].shadow_topR_jR),
                    self.parallelograms[0].pSlope))
            
            #take most restrictive a_corner and b_corner to be the one with the lowest z_top value
            self.a_corner = a_corner_list[np.argmin([intersect[0] for intersect in a_corner_list])]
            self.b_corner = b_corner_list[np.argmin([intersect[0] for intersect in b_corner_list])]
            
            #fixes for triangles
            if self.b_corner[0] < self.d_corner[0]:
                self.b_corner = self.d_corner
            if self.b_corner[1] > self.env.beam_axis_lim:
                self.b_corner = ( self.a_corner[0], self.env.beam_axis_lim)

            #print(self.a_corner)

        if self.c_corner[1] >= self.d_corner[1]:
            
            c_corner_list = []
            d_corner_list = []
            
            for layer in np.arange(self.env.num_layers-2)+1:
                #compute d_corners by intersecting c_line of layer 2, 3, 4 with d_line of layer 1
                d_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmin, 
                    self.parallelograms[layer].shadow_topL_jL),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmin, self.parallelograms[0].shadow_topL_jR),
                    self.parallelograms[0].pSlope))
                
                #compute c_corners by intersecting c_line of layer 2, 3, 4 with c_line of layer 1
                c_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmin, 
                    self.parallelograms[layer].shadow_topL_jL),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmin, self.parallelograms[0].shadow_topL_jL),
                    self.parallelograms[0].pSlope))
            
            #take most restrictive d_corner and c_corner to be the one with the highest z_top value
            self.d_corner = d_corner_list[np.argmax([intersect[0] for intersect in d_corner_list])]
            self.c_corner = c_corner_list[np.argmax([intersect[0] for intersect in c_corner_list])]

            #fixes for triangles
            if self.c_corner[0] > self.a_corner[0]:
                self.c_corner = self.a_corner
    
    def get_acceptanceCorners_v0(self):

        self.squareAcceptance = True # Ashutosh
        
        a_corner_list = [pgram.shadow_topR_jL for pgram in self.parallelograms]
        b_corner_list = [pgram.shadow_topR_jR for pgram in self.parallelograms]
        c_corner_list = [pgram.shadow_topL_jL for pgram in self.parallelograms]
        d_corner_list = [pgram.shadow_topL_jR for pgram in self.parallelograms]
        #print ("a_corner_list: ", a_corner_list)
        #print ("b_corner_list: ", b_corner_list)
        #print ("c_corner_list: ", c_corner_list)
        #print ("d_corner_list: ", d_corner_list)
        
        self.a_corner = (self.parallelograms[-1].top_layer_zmax, max(a_corner_list))
        self.b_corner = (self.parallelograms[-1].top_layer_zmax, min(b_corner_list))
        self.c_corner = (self.parallelograms[-1].top_layer_zmin, max(c_corner_list))
        self.d_corner = (self.parallelograms[-1].top_layer_zmin, min(d_corner_list))

        # is layer1 the most restrictive acceptance? 
        if np.argmax(a_corner_list) != 0:
            self.squareAcceptance = False
        if np.argmin(b_corner_list) != 0:
            self.squareAcceptance = False
        if np.argmax(c_corner_list) != 0:
            self.squareAcceptance = False
        if np.argmin(d_corner_list) != 0:
            self.squareAcceptance = False

        if self.b_corner <= self.a_corner:

            a_corner_list = []
            b_corner_list = []

            for layer in np.arange(self.env.num_layers-2)+1:
                #compute a_corners by intersecting b_line of layer 2, 3, 4 with a_line of layer 1
                a_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmax, 
                    self.parallelograms[layer].shadow_topR_jR),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmax, self.parallelograms[0].shadow_topR_jL),
                    self.parallelograms[0].pSlope))

                #compute b_corners by intersecting b_line of layer 2, 3, 4 with b_line of layer 1
                b_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmax, 
                    self.parallelograms[layer].shadow_topR_jR),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmax, self.parallelograms[0].shadow_topR_jR),
                    self.parallelograms[0].pSlope))
            
            #take most restrictive a_corner and b_corner to be the one with the lowest z_top value
            self.a_corner = a_corner_list[np.argmin([intersect[0] for intersect in a_corner_list])]
            self.b_corner = b_corner_list[np.argmin([intersect[0] for intersect in b_corner_list])]
            
            #fixes for triangles
            if self.b_corner[0] < self.d_corner[0]:
                self.b_corner = self.d_corner
            if self.b_corner[1] > self.env.beam_axis_lim:
                self.b_corner = ( self.a_corner[0], self.env.beam_axis_lim)

            #print(self.a_corner)

        if self.c_corner[1] >= self.d_corner[1]:
            
            c_corner_list = []
            d_corner_list = []
            
            for layer in np.arange(self.env.num_layers-2)+1:
                #compute d_corners by intersecting c_line of layer 2, 3, 4 with d_line of layer 1
                d_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmin, 
                    self.parallelograms[layer].shadow_topL_jL),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmin, self.parallelograms[0].shadow_topL_jR),
                    self.parallelograms[0].pSlope))
                
                #compute c_corners by intersecting c_line of layer 2, 3, 4 with c_line of layer 1
                c_corner_list.append(calc_line_intersection(
                    (self.parallelograms[layer].top_layer_zmin, 
                    self.parallelograms[layer].shadow_topL_jL),
                    self.parallelograms[layer].pSlope,
                    (self.parallelograms[0].top_layer_zmin, self.parallelograms[0].shadow_topL_jL),
                    self.parallelograms[0].pSlope))
            
            #take most restrictive d_corner and c_corner to be the one with the highest z_top value
            self.d_corner = d_corner_list[np.argmax([intersect[0] for intersect in d_corner_list])]
            self.c_corner = c_corner_list[np.argmax([intersect[0] for intersect in c_corner_list])]

            #fixes for triangles
            if self.c_corner[0] > self.a_corner[0]:
                self.c_corner = self.a_corner

    def plot(self, color='g'): 
        heights = self.env.radii
        for i in range(self.env.num_layers): 
            sp = self.superpoints[i] 
            
            max_height = self.env.radii[-1]
            
            plt.plot([sp.min, sp.max], [heights[i], heights[i]], c=color)
        
        plt.xticks(np.arange(-self.env.top_layer_lim, self.env.top_layer_lim, 10))
        plt.yticks(np.arange(0, max_height + 1, self.env.num_layers))
        # plt.show()
    
    def get_end_layer(self):
        lambdaZ_left_list = []
        lambdaZ_right_list = []
        for layer in range(self.env.num_layers):
            lambdaZ_left_list.append((self.superpoints[layer].min-self.apexZ0)/self.env.radii[layer])
            lambdaZ_right_list.append((self.superpoints[layer].max-self.apexZ0)/self.env.radii[layer])
        self.left_end_layer = np.argmax(lambdaZ_left_list)
        self.right_end_layer = np.argmin(lambdaZ_right_list)
        self.left_end_lambdaZ = max(lambdaZ_left_list)
        self.right_end_lambdaZ = min(lambdaZ_right_list)

class wedgeCover(): 
    
    def __init__(self, env:Environment, data:DataSet): 
        self.n_patches = 0 
        self.patches = [] 
        self.env = env 
        self.data = data 
        self.fitting_lines = [] 
        self.superPoints = [] 
        self.all_patches = []
        self.real_patch_list = []
        
    def add_patch(self, curr_patch:wedgePatch): 
        if self.n_patches == 0: 
            self.patches.append(curr_patch) 
            self.all_patches.append(curr_patch)
            self.real_patch_list.append(True)
            self.n_patches += 1 
        else:
            prev_patch = self.patches[-1] 
            prev_sp = prev_patch.superpoints 
            curr_sp = curr_patch.superpoints 
            
            for l in range(len(prev_sp)): 
                if (prev_sp[l].min != curr_sp[l].min) or (prev_sp[l].max != curr_sp[l].max): 
                    self.patches.append(curr_patch) 
                    self.all_patches.append(curr_patch)
                    self.real_patch_list.append(True)
                    self.n_patches += 1 
                    break
    
    def delete_patch(self, index):
        
        del self.patches[index]
        self.real_patch_list[index] = False

    def solve(self, lining:str = "makePatches_Projective", apexZ0=0, ppl = 16, nlines:int=100, leftRight:bool =True, show = True):

        # AVK remove identicalness of z-values of adjacent hits
        for row in range(self.env.num_layers):
            foundIdentical = False
            firstTime = True
            while (foundIdentical or firstTime):
                foundIdentical = False
                for x in range(len(self.data.array[row])-1):
                    if (self.data.array[row][x].z == self.data.array[row][x+1].z):
                        self.data.array[row][x+1].z += 0.00001
                        foundIdentical = True # search again to make sure there are no other identical ones 
                firstTime = False
                if foundIdentical:
                    self.data.array[row].sort(key=lambda point: point.z)

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
        
        if lining == 'makePatches_Projective_Leftright':
            try:
                for s in apexZ0:
                    self.makePatches_Projective(apexZ0=s, ppl = ppl, leftRight = True)
            except:
                self.makePatches_Projective(apexZ0=apexZ0, ppl = ppl, leftRight = True)
            return

        if lining == 'makePatches_ShadowQuilt_fromEdges':
            try:
                for s in apexZ0:
                    self.makePatches_ShadowQuilt_fromEdges(apexZ0=s, ppl = ppl, leftRight = leftRight)
            except:
                self.makePatches_ShadowQuilt_fromEdges(apexZ0=apexZ0, ppl = ppl, leftRight = leftRight)
            return
        if lining == 'makePatches_ShadowQuilt_fromCenter':
            try:
                for s in apexZ0:
                    self.makePatches_ShadowQuilt_fromCenter(apexZ0=s, ppl = ppl, leftRight = leftRight)
            except:
                self.makePatches_ShadowQuilt_fromCenter(apexZ0=apexZ0, ppl = ppl, leftRight = leftRight)
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

    def get_index_from_z(self, layer, z_value, alignment = 'closest'):
        layer_data = np.array([self.data.array[layer][x].z for x in range(len(self.data.array[layer]))])
        index = np.argmin(np.abs((layer_data - z_value)))
        if alignment == 'closest':
            return index
        
        if alignment == 'above':
            if layer_data[index] > z_value:
                return index
            else:
                return index + 1
            
        if alignment == 'below':
            if layer_data[index] < z_value:
                return index
            else:
                return index - 1

    def makePatches_ShadowQuilt_fromEdges_v0(self, apexZ0 = 0, stop = 1, ppl = 16, leftRight = True):
        """This method uses the geometry of shadows to generate patches based on superpoints
            the outer layer and the z0 of the collision point. 

        Args:
            apexZ0 (num, optional): Collision point on the z axis for the first patches Defaults to 0.
            stop (num, optional): Where to stop, normalized to 1. Defaults to 1.
            ppl (int, optional): Points per patch per layer. Defaults to 16.
            leftRight (bool, optional): If False, goes from right to left instead of left to right. Defaults to True.
        """

        '''
        #First, make list of contiguous superpoints from the outermost layer
        z_top_superpoint_edge_index = []
        top_layer_points = self.data.array[self.env.num_layers-1]
        top_row_list = np.array([top_layer_points[x].z for x in range(len(top_layer_points))])
        top_start_index = np.argmin(np.abs(top_row_list + self.env.top_layer_lim+self.env.boundaryPoint_offset))
        top_end_index = np.argmin(np.abs(top_row_list - (self.env.top_layer_lim+self.env.boundaryPoint_offset)))
        if leftRight == False:
            temp = top_start_index
            top_start_index = top_end_index
            top_end_index = temp

        num_points_z_top = abs(top_start_index-top_end_index)

        for sp in range(int(num_points_z_top/(ppl-1))):
            if leftRight == True:
                z_top_superpoint_edge_index.append((top_start_index+((ppl-1)*sp), top_start_index+((ppl-1)*sp)+(ppl-1)))
            else:
                z_top_superpoint_edge_index.append((top_start_index-((ppl-1)*sp)-(ppl-1), top_start_index-((ppl-1)*sp)))
        if num_points_z_top%(ppl-1) !=0:
            if leftRight == True:
                z_top_superpoint_edge_index.append((top_end_index-(ppl-1),top_end_index))
            else:
                z_top_superpoint_edge_index.append((top_end_index,top_end_index+(ppl-1)))

        
        print("total points in top layer: ", len(top_row_list))
        print("index list: ", z_top_superpoint_edge_index)
        print("start: ", top_start_index, "value: ", top_row_list[top_start_index])
        print("end: ", top_end_index, "value: ",  top_row_list[top_end_index])
        '''
        #self.makePatches_Projective(leftRight=False)

        initial_apexZ0 = self.env.beam_axis_lim
        apexZ0 = initial_apexZ0
        
        first_row_count = 0
        c_corner = np.inf
        bottom_layer_min = 0
        highest_c_corner = -self.env.top_layer_lim
        while (c_corner > -self.env.beam_axis_lim) & (bottom_layer_min > -self.env.trapezoid_edges[0]):
            #make top row in acceptance space by going to lower z0 from previous patch
            #patch is pushed up against right side trapezoid boundary in real space
            self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                         z_top = self.env.top_layer_lim + self.env.boundaryPoint_offset, 
                                         leftRight=False)
            #pick a value as next apexZ0 value
            apexZ0 = self.patches[-1].a_corner[1]
            c_corner = self.patches[-1].c_corner[1]
            first_row_count += 1
            highest_c_corner = max(highest_c_corner, self.patches[-1].c_corner[0])
            bottom_layer_min = self.patches[-1].superpoints[0].min
            #print("top row c_corner: ", c_corner)

        """
        #add top row rightmost patch again in order to make patches going down
        self.add_patch(self.patches[0])
        #make patches going down
        self.makePatches_Projective_Loop(leftRight=False, ppl=ppl, apexZ0 = self.patches[0].b_corner, stop = -1)
        #delete original first patch
        del self.patches[first_row_count]
        self.n_patches -= 1
        
        #go through each row and make rest of patches based on y value of rightmost column
        for column in range(first_row_count+1):
            apexZ0 = self.patches[first_row_count+column].a_corner
            while apexZ0 >= -self.env.beam_axis_lim: 
                self.makePatch_alignedToLine(apexZ0 = apexZ0, z_top = self.patches[first_row_count+column].superpoints[self.env.num_layers-1].max, leftRight=False, ppl = ppl)
                apexZ0 = self.patches[-1].a_corner   
    
        """

        initial_apexZ0 = -self.env.beam_axis_lim
        apexZ0 = initial_apexZ0
        last_row_count = 0
        b_corner = -np.inf
        bottom_layer_max = 0
        lowest_b_corner = self.env.top_layer_lim
        while (b_corner < self.env.beam_axis_lim) & (bottom_layer_max < self.env.trapezoid_edges[0]): 
            #make bottom row in acceptance space by going to higher z0 from previous patch
            #patch is pushed up against left side trapezoid boundary in real space
            self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                         z_top = -self.env.top_layer_lim - self.env.boundaryPoint_offset, 
                                         leftRight=True)
            #pick a value as next apexZ0 value
            apexZ0 = self.patches[-1].d_corner[1]
            b_corner = self.patches[-1].b_corner[1]
            last_row_count += 1
            bottom_layer_max = self.patches[-1].superpoints[0].max
            lowest_b_corner = min(lowest_b_corner, self.patches[-1].b_corner[0])
            #print("bottom row b_corner: ", b_corner)

        total_num_rows = 2
        z_top_eff_list = []
        for p in np.arange(first_row_count, first_row_count + last_row_count):
            #gets right end lambdas since we're going from bottom row and apexz0 value of each patch
            patch_end_lambdaZ = self.patches[p].right_end_lambdaZ
            patch_apexZ0 = self.patches[p].apexZ0

            #use limiting lambda and apexZ0 to compute effective z_top
            z_top_eff = self.env.radii[-1]*patch_end_lambdaZ + patch_apexZ0
            z_top_eff_list.append(z_top_eff)

        #pass z_top_min as next z_top value
        z_top_min =  min(z_top_eff_list)
        z_top_min_for_row = z_top_min 

        z_top_eff_list = []
        for p in np.arange(first_row_count):
            #gets right end lambdas since we're going from bottom row and apexz0 value of each patch
            patch_end_lambdaZ = self.patches[p].left_end_lambdaZ
            patch_apexZ0 = self.patches[p].apexZ0

            #use limiting lambda and apexZ0 to compute effective z_top
            z_top_eff = self.env.radii[-1]*patch_end_lambdaZ + patch_apexZ0
            z_top_eff_list.append(z_top_eff)

        #pass z_top_min as next z_top value
        z_top_max =  max(z_top_eff_list)
        z_top_max_for_row = z_top_max

        #while z_top_min_for_row < z_top_max_for_row:
        z_top_min = lowest_b_corner
        z_top_max = highest_c_corner
        for _ in range(1):


            initial_apexZ0 = -self.env.beam_axis_lim
            apexZ0 = initial_apexZ0
            previous_row_z_top =[]
            b_corner = -np.inf
            bottom_layer_max = 0
            while (b_corner < self.env.beam_axis_lim) & (bottom_layer_max < self.env.trapezoid_edges[0]):
                #building on top of the bottom row in acceptance space
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                             z_top = z_top_min, leftRight=True)
                patch_end_lambdaZ = self.patches[-1].right_end_lambdaZ
                previous_row_z_top.append(self.env.radii[-1]*patch_end_lambdaZ + apexZ0)
                #pick a value as next apexZ0 value
                plt.axhline(z_top_min, color = 'r')
                apexZ0 = self.patches[-1].d_corner[1]
                plt.axvline(apexZ0, color = 'k')
                plt.axvline(self.patches[-1].c_corner[1], color = 'b')
                plt.axhline(self.patches[-1].c_corner[0], color = 'b')

                
                last_b_corner = b_corner
                b_corner = self.patches[-1].b_corner[1]
                #bottom_layer_max = self.patches[-1].superpoints[0].max
                print(self.patches[-1].b_corner)
                z_top_min = self.patches[-1].d_corner[0]
                
                '''
                print("row b_corner: ", self.patches[-1].a_corner[1], b_corner, self.patches[-1].c_corner[1], apexZ0)
                print("z_top_min: ",z_top_min)
                if last_b_corner == b_corner:
                    for p in range(3):
                        for i in range(self.env.num_layers):
                            print(f"layer {i+1}: ", self.patches[-(p+1)].superpoints[i].min, self.patches[-(p+1)].superpoints[i].max)
                    #del self.patches[-1]
                    break
                '''
                
                #first_row_count += 1
            total_num_rows += 1
            z_top_min_for_row =  min(previous_row_z_top)
            z_top_min = z_top_min_for_row

            #if z_top_min_for_row > z_top_max_for_row:
            #    #print(total_num_rows)
            #    break

            initial_apexZ0 = self.env.beam_axis_lim
            apexZ0 = initial_apexZ0
            c_corner = np.inf
            bottom_layer_min = 0
            previous_row_z_top =[]
            print(lowest_b_corner)
            while (c_corner > -self.env.beam_axis_lim)  & (bottom_layer_min > -self.env.trapezoid_edges[0]):
                #make first row using last ppl points
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                             z_top = z_top_max, leftRight=False)
                patch_end_lambdaZ = self.patches[-1].left_end_lambdaZ
                previous_row_z_top.append(self.env.radii[-1]*patch_end_lambdaZ + apexZ0)
                #pick a value as next apexZ0 value
                apexZ0 = self.patches[-1].a_corner[1]
                c_corner = self.patches[-1].c_corner[1]
                #first_row_count += 1
                #print("row c_corner: ", c_corner)
                #bottom_layer_min = self.patches[-1].superpoints[0].min
                z_top_max = self.patches[-1].a_corner[0]

            z_top_max_for_row =  max(previous_row_z_top)
            z_top_max = z_top_max_for_row
            total_num_rows += 1
            #print(total_num_rows)

        '''
        #add top row rightmost patch again in order to make patches going down
        self.add_patch(self.patches[0])
        #make patches going down
        self.makePatches_Projective_Loop(leftRight=True, ppl=ppl, apexZ0 = self.patches[0].c_corner, stop = 1)
        #delete original first patch
        del self.patches[first_row_count]
        self.n_patches -= 1
        

        
        #go through each row and make rest of patches based on y value of rightmost column
        for column in range(first_row_count+1):
            apexZ0 = self.patches[first_row_count+column].a_corner
            while apexZ0 >= -self.env.beam_axis_lim: 
                self.makePatch_alignedToLine(apexZ0 = apexZ0, z_top = self.patches[first_row_count+column].superpoints[self.env.num_layers-1].max, leftRight=True, ppl = ppl)
                apexZ0 = self.patches[-1].a_corner       
        '''
    def makePatches_ShadowQuilt_fromEdges_v1(self, apexZ0 = 0, stop = 1, ppl = 16, leftRight = True):
        """This method uses the geometry of shadows to generate patches based on superpoints
            the outer layer and the z0 of the collision point. 

        Args:
            apexZ0 (num, optional): Collision point on the z axis for the first patches Defaults to 0.
            stop (num, optional): Where to stop, normalized to 1. Defaults to 1.
            ppl (int, optional): Points per patch per layer. Defaults to 16.
            leftRight (bool, optional): If False, goes from right to left instead of left to right. Defaults to True.
        """

        '''
        #First, make list of contiguous superpoints from the outermost layer
        z_top_superpoint_edge_index = []
        top_layer_points = self.data.array[self.env.num_layers-1]
        top_row_list = np.array([top_layer_points[x].z for x in range(len(top_layer_points))])
        top_start_index = np.argmin(np.abs(top_row_list + self.env.top_layer_lim+self.env.boundaryPoint_offset))
        top_end_index = np.argmin(np.abs(top_row_list - (self.env.top_layer_lim+self.env.boundaryPoint_offset)))
        if leftRight == False:
            temp = top_start_index
            top_start_index = top_end_index
            top_end_index = temp

        num_points_z_top = abs(top_start_index-top_end_index)

        for sp in range(int(num_points_z_top/(ppl-1))):
            if leftRight == True:
                z_top_superpoint_edge_index.append((top_start_index+((ppl-1)*sp), top_start_index+((ppl-1)*sp)+(ppl-1)))
            else:
                z_top_superpoint_edge_index.append((top_start_index-((ppl-1)*sp)-(ppl-1), top_start_index-((ppl-1)*sp)))
        if num_points_z_top%(ppl-1) !=0:
            if leftRight == True:
                z_top_superpoint_edge_index.append((top_end_index-(ppl-1),top_end_index))
            else:
                z_top_superpoint_edge_index.append((top_end_index,top_end_index+(ppl-1)))

        
        print("total points in top layer: ", len(top_row_list))
        print("index list: ", z_top_superpoint_edge_index)
        print("start: ", top_start_index, "value: ", top_row_list[top_start_index])
        print("end: ", top_end_index, "value: ",  top_row_list[top_end_index])
        '''
        #self.makePatches_Projective(leftRight=False)
        #initialize the method to make patches floor by floor, topFloor down and bottomFloor up
        z_top_max = self.env.beam_axis_lim
        z_top_min = -self.env.top_layer_lim

        initial_apexZ0 = self.env.beam_axis_lim
        initial_apexZ0 = self.env.trapezoid_edges[0] # switched from z0 to z1
        print ("initial_apexZ0: ", initial_apexZ0)
        apexZ0 = initial_apexZ0        
        first_row_count = 0
        c_corner = np.inf
        bottom_layer_min = 0
        while (c_corner > -self.env.trapezoid_edges[0]) & (bottom_layer_min > -self.env.trapezoid_edges[0]):
            #make top row in acceptance space by going to lower z1 from previous patch
            #patch is pushed up against right side trapezoid boundary in real space
            self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                         z_top = self.env.top_layer_lim + self.env.boundaryPoint_offset, 
                                         leftRight=False)
            print("after seedPatch, N_patches: ", len(self.patches))
            seed_apexZ0 = apexZ0
            #pick "a" value as next apexZ0 value
            apexZ0 = self.patches[-1].a_corner[1]
            c_corner = self.patches[-1].c_corner[1]
            print ("a: ", self.patches[-1].a_corner, "b: ", self.patches[-1].b_corner, "c: ", self.patches[-1].c_corner, "d: ", self.patches[-1].d_corner)
            first_row_count += 1
            bottom_layer_min = self.patches[-1].superpoints[0].min
            z_top_min = max(z_top_min, self.patches[-1].parallelograms[0].top_layer_zmin)
            #z_top_min = self.patches[-1].parallelograms[0].top_layer_zmin
            #print("before squareAcceptance check, z_top_min: ", z_top_min)
            if self.patches[-1].squareAcceptance == False:
                complementary_apexZ0 = self.patches[-1].parallelograms[0].shadow_topL_jL
                #print("complementary_apexZ0 ", complementary_apexZ0)
                # check that outermost superpoint is not too long
                self.makePatch_alignedToLine(apexZ0 = seed_apexZ0, ppl = ppl,
                                             z_top = self.env.top_layer_lim + self.env.boundaryPoint_offset,
                                             leftRight=False, double_middleLayers_ppl = True)
                if self.patches[-1].squareAcceptance == False:
                    # outermost superpoint is too long, cut it short
                    check_c_list = []
                    for j, superpoint in enumerate(self.patches[-1].superpoints[1:-1], start=2):
                        check_c_list.append(self.patches[-1].straightLineProjectorFromLayer1(self.patches[-1].superpoints[0].min, superpoint.min, j))
                        #print("layer ", j, " self.patches[-1].superpoints[0].min, superpoint.min ", self.patches[-1].superpoints[0].min, superpoint.min)
                    #print("check_c_list ", check_c_list)
                    # smallest c value corresponds to shortest middleLayer
                    shortestSP = np.argmin(check_c_list) + 1
                    #print("shortestSP layer ", shortestSP+1)
                    # project SP1 and shortestSP to outermost layer
                    z_top_min = max(z_top_min, self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].superpoints[0].min,
                                                                                     self.patches[-1].superpoints[shortestSP].min,
                                                                                     1, shortestSP+1, self.env.num_layers))
                    complementary_apexZ0 = self.patches[-1].straightLineProjector(z_top_min, self.patches[-1].superpoints[0].min, 1)
                    #print("shortened, complementary_apexZ0 ", complementary_apexZ0, " updated z_top_min ", z_top_min)
                # delete the check patch
                #del self.patches[-1]
                self.delete_patch(-1)
                self.n_patches -= 1
                # now make actual complementary patch 
                self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_min, leftRight=True)
                c_corner = self.patches[-1].c_corner[1]
                #print("complementary's c corner: ",self.patches[-1].c_corner, " d corner: ", self.patches[-1].d_corner)
                pass
            #print("maybe after complementaryPatch, N_patches: ", len(self.patches))
        print("after topFloor, N_patches: ", len(self.patches))
        
        initial_apexZ0 = -self.env.beam_axis_lim
        initial_apexZ0 = -self.env.trapezoid_edges[0] # switched from z0 to z1    
        apexZ0 = initial_apexZ0
        #last_row_count = 0
        b_corner = -np.inf
        bottom_layer_max = 0
        while ((b_corner < self.env.trapezoid_edges[0]) & (bottom_layer_max < self.env.trapezoid_edges[0])): 
            #make bottom row in acceptance space by going to higher z0 from previous patch
            #patch is pushed up against left side trapezoid boundary in real space
            self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                         z_top = -self.env.top_layer_lim - self.env.boundaryPoint_offset, 
                                         leftRight=True)
            print("after creating seedPatch, N_patches: ", len(self.patches))
            seed_apexZ0 = apexZ0
            z_top_max = min(z_top_max, self.patches[-1].parallelograms[0].top_layer_zmax)    
           #pick a value as next apexZ0 value
            apexZ0 = self.patches[-1].d_corner[1]
            b_corner = self.patches[-1].b_corner[1]
            #last_row_count += 1
            bottom_layer_max = self.patches[-1].superpoints[0].max
            if self.patches[-1].squareAcceptance == False:
                complementary_apexZ0 = self.patches[-1].parallelograms[0].shadow_topR_jR
                print("complementary_apexZ0: ", complementary_apexZ0, "z_top_max: ", z_top_max)                               
                # check that outermost superpoint is not too long                                                                                        
                for layer in range(self.env.num_layers):                                                                                                
                    print("layer ", layer+1," superpoint min: ", self.patches[-1].superpoints[layer].min, " max: ", self.patches[-1].superpoints[layer].max)                                                                                                             
                self.makePatch_alignedToLine(apexZ0 = seed_apexZ0, ppl = ppl,
                                             z_top = -self.env.top_layer_lim - self.env.boundaryPoint_offset,
                                             leftRight=True, double_middleLayers_ppl = True)
                #for layer in range(self.env.num_layers):   
                    #print("layer ", layer+1," superpoint min: ", self.patches[-1].superpoints[layer].min, " max: ", self.patches[-1].superpoints[layer].max)                                                                                                         
                print("after creating checkPatch, N_patches: ", len(self.patches))    
                if self.patches[-1].squareAcceptance == False:
                    # outermost superpoint is too long, cut it short                                                                                     
                    check_b_list = []
                    for j, superpoint in enumerate(self.patches[-1].superpoints[1:-1], start=2):
                        check_b_list.append(self.patches[-1].straightLineProjectorFromLayer1(self.patches[-1].superpoints[0].max, superpoint.max, j))
                        #print("layer ", j, " self.patches[-1].superpoints[0].min, superpoint.min ", self.patches[-1].superpoints[0].max, superpoint.max) 
                    #print("check_b_list ", check_b_list) 
                    # smallest c value corresponds to shortest middleLayer         
                    shortestSP = np.argmax(check_b_list) + 1
                    #print("shortestSP layer ", shortestSP+1)                                                                                            
                    # project SP1 and shortestSP to outermost layer         
                    z_top_max = min(z_top_max, self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].superpoints[0].max,
                                                                                     self.patches[-1].superpoints[shortestSP].max,
                                                                                     1, shortestSP+1, self.env.num_layers))
                    complementary_apexZ0 = self.patches[-1].straightLineProjector(z_top_max, self.patches[-1].superpoints[0].max, 1)
                    print("shortened: updated complementary_apexZ0 ", complementary_apexZ0, " updated z_top_max ", z_top_max, len(self.patches))  
                # delete the check patch                               
                #del self.patches[-1]
                self.delete_patch(-1)
                self.n_patches -= 1
                print("after deleting checkPatch, N_patches: ", len(self.patches))
                # now make actual complementary patch  
                self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_max, leftRight=False)
                b_corner = self.patches[-1].b_corner[1]
                print("after creating complementary Patch, N_patches: ", len(self.patches), "b_corner: ", self.patches[-1].b_corner)
        print("after groundFloor, N_patches: ", len(self.patches))
        '''    
        #startCommentAshutosh    
        total_num_rows = 2
        z_top_eff_list = []
        for p in np.arange(first_row_count, first_row_count + last_row_count):
            #gets right end lambdas since we're going from bottom row and apexz0 value of each patch
            patch_end_lambdaZ = self.patches[p].right_end_lambdaZ
            patch_apexZ0 = self.patches[p].apexZ0

            #use limiting lambda and apexZ0 to compute effective z_top
            z_top_eff = self.env.radii[-1]*patch_end_lambdaZ + patch_apexZ0
            z_top_eff_list.append(z_top_eff)

        #pass z_top_min as next z_top value
        z_top_min =  min(z_top_eff_list)
        z_top_min_for_row = z_top_min 

        z_top_eff_list = []
        for p in np.arange(first_row_count):
            #gets right end lambdas since we're going from bottom row and apexz0 value of each patch
            patch_end_lambdaZ = self.patches[p].left_end_lambdaZ
            patch_apexZ0 = self.patches[p].apexZ0

            #use limiting lambda and apexZ0 to compute effective z_top
            z_top_eff = self.env.radii[-1]*patch_end_lambdaZ + patch_apexZ0
            z_top_eff_list.append(z_top_eff)

        #pass z_top_min as next z_top value
        z_top_max =  max(z_top_eff_list)
        z_top_max_for_row = z_top_max

        #while z_top_min_for_row < z_top_max_for_row:
        z_top_min = lowest_b_corner
        z_top_max = highest_c_corner
        
        #endCommentAshutosh 
        '''                               
#        for _ in range(1):
        #print ("after top and bottom floor, z_top_max: " , z_top_max, " z_top_min: ", z_top_min, "N_patches: ", len(self.patches))
        while (z_top_max < z_top_min):
            initial_apexZ0 = -self.env.trapezoid_edges[0]
            apexZ0 = initial_apexZ0
            #previous_row_z_top =[]
            b_corner = -np.inf
            bottom_layer_max = 0
            while (b_corner < self.env.trapezoid_edges[0]) & (bottom_layer_max < self.env.trapezoid_edges[0]):
                #building on top of the bottom row in acceptance space
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, z_top = z_top_max, leftRight=True)
                #print ("z_top_max: ", z_top_max, "apexZ0: ", apexZ0, "N_patches: ", len(self.patches))
                #for layer in range(self.env.num_layers):
                    #print("layer ", layer+1," superpoint min: ", self.patches[-1].superpoints[layer].min, " max: ", self.patches[-1].superpoints[layer].max)
#                patch_end_lambdaZ = self.patches[-1].right_end_lambdaZ
#                previous_row_z_top.append(self.env.radii[-1]*patch_end_lambdaZ + apexZ0)
                #pick a value as next apexZ0 value
#                plt.axhline(z_top_min, color = 'r')
                apexZ0 = self.patches[-1].d_corner[1] + self.env.parallelogramSlopes[0]*(z_top_max - self.patches[-1].d_corner[0])
#                plt.axvline(apexZ0, color = 'k')
#                plt.axvline(self.patches[-1].c_corner[1], color = 'b')
#                plt.axhline(self.patches[-1].c_corner[0], color = 'b')               
                last_b_corner = b_corner
                b_corner = self.patches[-1].b_corner[1]
                bottom_layer_max = self.patches[-1].superpoints[0].max
                z_top_max_upperFloor = self.patches[-1].parallelograms[0].top_layer_zmax
                #print("seedPatch b_corner: ", self.patches[-1].b_corner, "z_top_min: ", z_top_min)
                if (self.patches[-1].squareAcceptance == False) and (self.patches[-1].b_corner[0] <= z_top_min):
                    #print ("noSqAc, z_top_max_upperFloor: ", z_top_max_upperFloor, "b corner: ", self.patches[-1].parallelograms[0].shadow_topR_jR,len(self.patches))
                    self.makePatch_alignedToLine(apexZ0 = self.patches[-1].parallelograms[0].shadow_topR_jR, ppl = ppl, z_top = z_top_max_upperFloor, leftRight=False)
                    #print ("noSqAc b_corner: ", self.patches[-1].b_corner, len(self.patches))
                    b_corner = self.patches[-1].b_corner[1]
                    pass
            z_top_max = z_top_max_upperFloor
            #print("climbing one Floor, N_patches: ", len(self.patches))

            '''
                print("row b_corner: ", self.patches[-1].a_corner[1], b_corner, self.patches[-1].c_corner[1], apexZ0)
                print("z_top_min: ",z_top_min)
                if last_b_corner == b_corner:
                    for p in range(3):
                        for i in range(self.env.num_layers):
                            print(f"layer {i+1}: ", self.patches[-(p+1)].superpoints[i].min, self.patches[-(p+1)].superpoints[i].max)
                    #del self.patches[-1]
                    break
                
                #first_row_count += 1
            total_num_rows += 1
            z_top_min_for_row =  min(previous_row_z_top)
            z_top_min = z_top_min_for_row

            '''
            if (z_top_max < z_top_min):
                initial_apexZ0 = self.env.trapezoid_edges[0]
                apexZ0 = initial_apexZ0
                c_corner = np.inf
                bottom_layer_min = 0
                previous_row_z_top =[]
                while (c_corner > -self.env.trapezoid_edges[0])  & (bottom_layer_min > -self.env.trapezoid_edges[0]):
                    #make first row using last ppl points
                    self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, z_top = z_top_min, leftRight=False)
                    #patch_end_lambdaZ = self.patches[-1].left_end_lambdaZ
                    #previous_row_z_top.append(self.env.radii[-1]*patch_end_lambdaZ + apexZ0)
                    #pick a value as next apexZ0 value
                    apexZ0 = self.patches[-1].a_corner[1]
                    c_corner = self.patches[-1].c_corner[1]
                    #first_row_count += 1
                    bottom_layer_min = self.patches[-1].superpoints[0].min
                    z_top_min_roof = self.patches[-1].parallelograms[0].top_layer_zmin
                    if self.patches[-1].squareAcceptance == False:
                        self.makePatch_alignedToLine(apexZ0 = self.patches[-1].parallelograms[0].shadow_topL_jL, ppl = ppl, z_top = z_top_min_roof, leftRight=True)
                        c_corner = self.patches[-1].c_corner[1]
                        pass
                    
                z_top_min = z_top_min_roof
            #z_top_max_for_row =  max(previous_row_z_top)
            #z_top_max = z_top_max_for_row
            #total_num_rows += 1
            #print(total_num_rows)
         
#        print(z_top_min,z_top_max)
    def makePatches_ShadowQuilt_fromEdges(self, apexZ0 = 0, stop = 1, ppl = 16, leftRight = True):
        """This method uses the geometry of shadows to generate patches based on superpoints
            the outer layer and the z0 of the collision point. 

        Args:
            apexZ0 (num, optional): starting point on the z1 axis for the first patches Defaults to 0.
            stop (num, optional): Where to stop, normalized to 1. Defaults to 1.
            ppl (int, optional): Points per patch per layer. Defaults to 16.
            leftRight (bool, optional): If False, goes from right to left instead of left to right. Defaults to True.
        """
        #self.makePatches_Projective(leftRight=False)
        #initialize the method to make patches floor by floor, topFloor down and bottomFloor up

        fix42 = True
        apexZ0 = self.env.trapezoid_edges[0] # switched from z0 to z1
        while (apexZ0 > -self.env.trapezoid_edges[0]):
            z_top_min = -self.env.top_layer_lim
            complementary_apexZ0 = 0
            first_row_count = 0
            c_corner = np.inf
            z_top_max = self.env.top_layer_lim + self.env.boundaryPoint_offset
            if (len(self.patches) > 0):
                z_top_max = min(z_top_max,
                            self.patches[-1].straightLineProjectorFromLayerIJtoK(-self.env.beam_axis_lim,apexZ0,0,1,self.env.num_layers))
            nPatchesInColumn = 0
            projectionOfCcornerToBeam = 0
            while (c_corner > -self.env.trapezoid_edges[self.env.num_layers-1]) and (nPatchesInColumn<100000000) and (projectionOfCcornerToBeam < self.env.beam_axis_lim):
                nPatchesInColumn += 1
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, z_top = z_top_max, leftRight=False)
                print ('top layer from ', self.patches[-1].superpoints[self.env.num_layers-1].max, ' to ', self.patches[-1].superpoints[self.env.num_layers-1].min, ' z_top_max: ', z_top_max)
                print('original:',self.patches[-1].a_corner, 'for patch',len(self.patches))
                print('original:',self.patches[-1].b_corner)
                print('original:',self.patches[-1].c_corner)
                print('original:',self.patches[-1].d_corner)
                for j, sp in enumerate(self.patches[-1].superpoints[1:-1], start=2):
                    print (j,'superpoint:',sp.min,sp.max,'shadowTop from L1Max:',self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].superpoints[0].max,sp.min,1,j,self.env.num_layers),self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].superpoints[0].max,sp.max,1,j,self.env.num_layers),'from L1Min:',self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].superpoints[0].min,sp.min,1,j,self.env.num_layers),self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].superpoints[0].min,sp.max,1,j,self.env.num_layers))
                #original_c = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].c_corner[1], 'below')
                #original_d = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].d_corner[1], 'below')
                original_c = self.patches[-1].c_corner[1]
                original_d = self.patches[-1].d_corner[1]
                c_corner = original_c
                repeat_patch = False
                repeat_original = False
                if len(self.patches) > 2:
                    repeat_original = (self.patches[-1].superpoints[self.env.num_layers-1] == self.patches[-3].superpoints[self.env.num_layers-1]) and (self.patches[-1].superpoints[0] == self.patches[-3].superpoints[0]) and (self.patches[-1].superpoints[1] == self.patches[-3].superpoints[1]) and (self.patches[-1].superpoints[2] == self.patches[-3].superpoints[2]) and (self.patches[-1].superpoints[3] == self.patches[-3].superpoints[3])
                #if (self.patches[-1].triangleAcceptance == True):
                    #original_c = self.patches[-1].b_corner[1]
                seed_apexZ0 = apexZ0
                projectionOfCcornerToBeam = self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].c_corner[1],self.patches[-1].c_corner[0],self.env.num_layers,1,0)
                squarePatch_alternate1 = ((self.patches[-1].a_corner[1] > z_top_max) and (self.patches[-1].b_corner[1] > z_top_max) and self.patches[-1].flatBottom)
                squarePatch_alternate2 = ((self.patches[-1].a_corner[1] > z_top_max) and self.patches[-1].flatBottom)
                notChoppedPatch = (self.patches[-1].squareAcceptance) or squarePatch_alternate1 or squarePatch_alternate2
                madeComplementaryPatch = False
                nPatchesAtOriginal = len(self.patches)
                print('squareAcceptance: ', self.patches[-1].squareAcceptance, 'triangleAcceptance: ', self.patches[-1].triangleAcceptance, ' projectionOfCcornerToBeam: ', projectionOfCcornerToBeam,'notChoppedPatch',notChoppedPatch)
                if (not notChoppedPatch) and (self.patches[-1].c_corner[1] > -self.env.trapezoid_edges[self.env.num_layers-1]) and (projectionOfCcornerToBeam < self.env.beam_axis_lim):
                    complementary_apexZ0 = self.patches[-1].superpoints[0].min
                    if (self.patches[-1].triangleAcceptance == True) and not(repeat_original):
                        z_top_min = self.patches[-1].d_corner[1]
                    else:
                        print('z_top_min before:', z_top_min, 'superpoints[self.env.num_layers-1].min:', self.patches[-1].superpoints[self.env.num_layers-1].min)
                        z_top_min = max(-self.env.top_layer_lim, self.patches[-1].superpoints[self.env.num_layers-1].min)
                    self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_min, leftRight=True)
                    madeComplementaryPatch = True
                    print('complementary: ', self.patches[-1].a_corner, ' for z_top_min:', z_top_min)
                    print('complementary: ', self.patches[-1].b_corner,'for patch',len(self.patches))
                    print('complementary: ', self.patches[-1].c_corner)
                    print('complementary: ', self.patches[-1].d_corner)
                    #complementary_a = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].a_corner[1], 'above')
                    #complementary_b = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].b_corner[1], 'above')
                    complementary_a = self.patches[-1].a_corner[1]
                    complementary_b = self.patches[-1].b_corner[1]
                    white_space_height = max(original_c - complementary_a, original_d - complementary_b)
                    previous_white_space_height = -1
                    counter = 0
                    counterUpshift = 0
                    current_z_top_index = -1
                    previous_z_top_min = -999
                    #while ((white_space_height > 0) or (abs(white_space_height) > 10)) and (counter < 5):
                    #while (counterUpshift < 100) and (white_space_height != 0) and ((counter < 15) or (white_space_height > 0)) and (self.patches[-1].c_corner[1] > -self.env.trapezoid_edges[self.env.num_layers-1]):
                    
                    #while not((white_space_height < 0) and (previous_white_space_height >= 0)) and ((self.patches[-1].c_corner[1] > -self.env.trapezoid_edges[self.env.num_layers-1]) or (white_space_height > 0)) and (current_z_top_index < (len(self.data.array[self.env.num_layers-1])-1)) and (self.patches[-2].triangleAcceptance == False) :
                    while not((white_space_height <= 0) and (previous_white_space_height >= 0)) and (abs(white_space_height)>0.000001) and ((self.patches[-1].c_corner[1] > -self.env.trapezoid_edges[self.env.num_layers-1]) or (white_space_height > 0)) and (current_z_top_index < (len(self.data.array[self.env.num_layers-1])-1)) and not(repeat_patch) and not(repeat_original):
                        print()
                        if (len(self.patches) > 2):
                            print('original c:', original_c, ' ', self.patches[-2].c_corner[1], '|| original d:', original_d, ' ', self.patches[-2].d_corner[1])
                        print('complementary_a:', complementary_a, ' ', self.patches[-1].a_corner[1], ' || complementary_b:', complementary_b, ' ', self.patches[-1].b_corner[1])
                        current_z_top_index = self.get_index_from_z(self.env.num_layers-1, z_top_min)
                        print('current white_space_height: ', white_space_height)
                        print('counter: ',counter, ' counterUpshift: ', counterUpshift)
                        print('orig_ztop: ', current_z_top_index, 'orig_z_top_min: ', z_top_min)
                        current_z_i_index = tuple(self.get_index_from_z(
                            layer,
                            self.patches[-1].straightLineProjectorFromLayerIJtoK(complementary_apexZ0,z_top_min,1,self.env.num_layers,layer+1)
                            ) for layer in range(self.env.num_layers))
                        if (z_top_min == previous_z_top_min):
                            current_z_top_index += 1
                            new_z_i_index = tuple(oldIndex+1 for oldIndex in current_z_i_index)
                        previous_z_top_min = z_top_min
                        if (white_space_height < 0):
                            counter +=1
                            current_z_top_index -= 1
                            new_z_i_index = tuple(oldIndex-1 for oldIndex in current_z_i_index)
                        else:
                            counterUpshift += 1
                            current_z_top_index += 1
                            new_z_i_index = tuple(oldIndex+1 for oldIndex in current_z_i_index)
                        current_z_top_index = min(current_z_top_index,len(self.data.array[self.env.num_layers-1])-1)
                        new_z_i_index = tuple(min(z_i_index,len(self.data.array[layer])-1) for layer, z_i_index in enumerate(new_z_i_index)) 
                        new_z_i_index = tuple(max(z_i_index,0) for layer, z_i_index in enumerate(new_z_i_index))
                        new_z_i = tuple(self.data.array[layer][new_z_i_index[layer]].z for layer in range(self.env.num_layers))
                        new_z_i_atTop = tuple(self.patches[-1].straightLineProjectorFromLayerIJtoK(complementary_apexZ0,new_z_i[layer],1,layer+1,self.env.num_layers) for layer in range(1,self.env.num_layers))
                        layerWithSmallestShift = 1 + np.argmin(np.abs(np.array(new_z_i_atTop)-previous_z_top_min))
                        for layer in range(self.env.num_layers-1):
                            print (layer+1, ' new_z_i_atTop: ', new_z_i_atTop[layer], ' shift_i_ztop: ', new_z_i_atTop[layer]-previous_z_top_min,
                                ' layerWithSmallestShift: ', layerWithSmallestShift)
                        z_top_min = self.data.array[self.env.num_layers-1][current_z_top_index].z
                        z_top_min = new_z_i_atTop[layerWithSmallestShift-1] # AVK try smallest shift
                        if abs(z_top_min-previous_z_top_min) < 0.000001:
                            z_top_min = self.data.array[self.env.num_layers-1][current_z_top_index].z
                        if abs(z_top_min-previous_z_top_min) < 0.000001:
                            z_top_min = self.data.array[self.env.num_layers-2][current_z_top_index].z
                        if abs(z_top_min-previous_z_top_min) < 0.000001:
                            z_top_min = self.data.array[self.env.num_layers-3][current_z_top_index].z
                        if ((z_top_min-previous_z_top_min)*(white_space_height)) < 0:
                            z_top_min = new_z_i_atTop[self.env.num_layers-2]
                        print('new_def_z_top_min_diff:',z_top_min-self.data.array[self.env.num_layers-1][current_z_top_index].z)
                        print('new_ztop_index: ', current_z_top_index, ' new_z_i_index: ', new_z_i_index, ' new_z_top_min: ', z_top_min, ' shift_ztop:', z_top_min-previous_z_top_min)
                        nPatchesAtComplementary = len(self.patches)
                        if (nPatchesAtComplementary > nPatchesAtOriginal):
                            print('deleted complementary: ', self.patches[-1].a_corner, 'for patch',len(self.patches))
                            print('deleted complementary: ', self.patches[-1].b_corner)
                            print('deleted complementary: ', self.patches[-1].c_corner)
                            print('deleted complementary: ', self.patches[-1].d_corner)
                            #del self.patches[-1]
                            self.delete_patch(-1)
                            self.n_patches -= 1
                        self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_min, leftRight=True)
                        #complementary_a = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].a_corner[1], 'above')
                        #complementary_b = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].b_corner[1], 'above')
                        complementary_a = self.patches[-1].a_corner[1]
                        complementary_b = self.patches[-1].b_corner[1]
                        previous_white_space_height = white_space_height
                        white_space_height = max(original_c - complementary_a, original_d - complementary_b)
                        print('complementary_a:', complementary_a, ' ', self.patches[-1].a_corner[1], ' || complementary_b:', complementary_b, ' ', self.patches[-1].b_corner[1], ' new z_top_min: ', z_top_min)
                        print('new white_space_height: ', white_space_height)
                        print('adjusted complementary: ', self.patches[-1].a_corner, ' for z_top_min:', z_top_min)
                        print('adjusted complementary: ', self.patches[-1].b_corner, 'for patch',len(self.patches))
                        print('adjusted complementary: ', self.patches[-1].c_corner)
                        print('adjusted complementary: ', self.patches[-1].d_corner)
                    if (self.n_patches > 3) and fix42:
                        if (self.patches[-1].superpoints[self.env.num_layers-1] == self.patches[-3].superpoints[self.env.num_layers-1]) and (self.patches[-1].superpoints[0] == self.patches[-3].superpoints[0]) and (self.patches[-1].superpoints[1] == self.patches[-3].superpoints[1]) and (self.patches[-1].superpoints[2] == self.patches[-3].superpoints[2]) and (self.patches[-1].superpoints[3] == self.patches[-3].superpoints[3]):
                            repeat_patch = True
                            print (self.patches[-1].superpoints[self.env.num_layers-1].min, self.patches[-1].superpoints[self.env.num_layers-1].max, ' repeat_patch: ', repeat_patch)
                            #del self.patches[-1]
                            self.delete_patch[-1]
                            self.n_patches -= 1
                            current_z_top_index -= 1
                            z_top_min = self.data.array[self.env.num_layers-1][current_z_top_index].z
                            z_top_min = new_z_i_atTop[layerWithSmallestShift-1] # AVK try smallest shift    
                            self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_min, leftRight=True)

                c_corner = self.patches[-1].c_corner[1]
                projectionOfCcornerToBeam = self.patches[-1].straightLineProjectorFromLayerIJtoK(c_corner,self.patches[-1].c_corner[0],self.env.num_layers,1,0)

                saved_apexZ0 = self.patches[-1].c_corner[0]

                if madeComplementaryPatch:
                    self.patches[-1].getShadows(z_top_min,z_top_max)
                    self.patches[-2].getShadows(z_top_min,z_top_max)
                    #abAlignment = self.patches[-1].a_corner[1]-self.patches[-2].b_corner[1]
                    #cdAlignment = self.patches[-1].c_corner[1]-self.patches[-2].d_corner[1]
                    #print('alignment:',self.patches[-1].a_corner[1]-self.patches[-2].b_corner[1], self.patches[-1].c_corner[1]-self.patches[-2].d_corner[1])
                    #print('z_top_min:', z_top_min, 'z_top_max:', z_top_max, 'shadow_fromTopToInnermost_topL_jL:', self.patches[-1].shadow_fromTopToInnermost_topL_jL, 'shadow_fromTopToInnermost_topL_jR:', self.patches[-1].shadow_fromTopToInnermost_topL_jR, 'shadow_fromTopToInnermost_topR_jL:', self.patches[-1].shadow_fromTopToInnermost_topR_jL, 'shadow_fromTopToInnermost_topR_jR:', self.patches[-1].shadow_fromTopToInnermost_topR_jR)
                    original_topR_jL = self.patches[-2].shadow_fromTopToInnermost_topR_jL
                    originalPartialTop = (original_topR_jL > complementary_apexZ0) and (original_topR_jL < apexZ0) and (abs(self.patches[-2].straightLineProjectorFromLayerIJtoK(original_topR_jL,z_top_max,1,self.env.num_layers,0))<20*self.env.beam_axis_lim)
                    original_topL_jL = self.patches[-2].shadow_fromTopToInnermost_topL_jL
                    originalPartialBottom = (original_topL_jL > complementary_apexZ0) and (original_topL_jL < apexZ0) and (abs(self.patches[-2].straightLineProjectorFromLayerIJtoK(original_topL_jL,z_top_min,1,self.env.num_layers,0))<20*self.env.beam_axis_lim)
                    complementary_topR_jR = self.patches[-1].shadow_fromTopToInnermost_topR_jR
                    complementaryPartialTop = (complementary_topR_jR > complementary_apexZ0) and (complementary_topR_jR < apexZ0) and (abs(self.patches[-1].straightLineProjectorFromLayerIJtoK(complementary_topR_jR,z_top_max,1,self.env.num_layers,0))<20*self.env.beam_axis_lim)
                    complementary_topL_jR = self.patches[-1].shadow_fromTopToInnermost_topL_jR
                    complementaryPartialBottom = (complementary_topL_jR > complementary_apexZ0) and (complementary_topL_jR < apexZ0) and (abs(self.patches[-1].straightLineProjectorFromLayerIJtoK(complementary_topL_jR,z_top_min,1,self.env.num_layers,0))<20*self.env.beam_axis_lim)

                    horizontalShiftTop = original_topR_jL - complementary_topR_jR
                    horizontalShiftBottom = original_topL_jL - complementary_topL_jR

                    complementary_topR_jL = self.patches[-1].shadow_fromTopToInnermost_topR_jL
                    complementary_topL_jL = self.patches[-1].shadow_fromTopToInnermost_topL_jL
                    original_topR_jR = self.patches[-2].shadow_fromTopToInnermost_topR_jR
                    original_topL_jR = self.patches[-2].shadow_fromTopToInnermost_topL_jR

                    originalSaved_topR_jR = original_topR_jR
                    originalSaved_topL_jR = original_topL_jR
                    originalSaved_topR_jL = original_topR_jL
                    originalSaved_topL_jL = original_topL_jL
                    
                    complementarySaved_topR_jL = complementary_topR_jL
                    complementarySaved_topL_jL = complementary_topL_jL
                    
                    horizontalOverlapTop = max(complementary_topR_jL - original_topR_jL, complementary_topR_jR - original_topR_jR)
                    horizontalOverlapBottom = max(complementary_topL_jL - original_topL_jL, complementary_topL_jR - original_topL_jR)
                    horizontalOverlapTop = -1
                    horizontalOverlapBottom = -1
                    newGapTop = -0.000001
                    newGapBottom = -0.000001
                    
                    makeHorizontallyShiftedPatch = False
                    shifted_Align = apexZ0
                    doShiftedPatch = True
                    # decide whether to shift original patch or complementary patch
                    # compute z0 associated with b-corner of original patch and c-corner of complementary patch
                    z0_original_bCorner = self.patches[-2].straightLineProjectorFromLayerIJtoK(apexZ0,z_top_max,1,self.env.num_layers,0)
                    z0_complementary_cCorner = self.patches[-1].straightLineProjectorFromLayerIJtoK(complementary_apexZ0,z_top_min,1,self.env.num_layers,0)
                    shiftOriginal = True
                    if (z0_original_bCorner < 0):
                        shiftOriginal = False
                        shifted_Align = complementary_apexZ0
                    if (z0_complementary_cCorner > 0):
                        shiftOriginal = True
                        shifted_Align = apexZ0

                    if (horizontalShiftTop > 0 or horizontalShiftBottom > 0):
                        print('originalPartialTop:',originalPartialTop,'complementaryPartialTop:',complementaryPartialTop,'originalPartialBottom:',originalPartialBottom,'complementaryPartialBottom:',complementaryPartialBottom, original_topR_jL, original_topL_jL, complementary_topR_jR, complementary_topL_jR,'horizontalOverlapTop:',horizontalOverlapTop,'horizontalOverlapBottom:',horizontalOverlapBottom)
                    while ((horizontalShiftTop > 0 and originalPartialTop and complementaryPartialTop) or (horizontalShiftBottom > 0 and originalPartialBottom and complementaryPartialBottom)) and doShiftedPatch and (horizontalOverlapTop <= 0) and (horizontalOverlapBottom <= 0) and (newGapTop<0 or newGapBottom<0):
                        print('horizontalShifts:',horizontalShiftTop,horizontalShiftBottom, 'shifted_Align:',shifted_Align)
                        newZtop = z_top_max
                        if shiftOriginal:
                            shifted_Align -= max(horizontalShiftTop,horizontalShiftBottom) #+ min(newGapTop,newGapBottom)
                        else :
                            shifted_Align += max(horizontalShiftTop,horizontalShiftBottom)
                            newZtop = z_top_min
                        if (makeHorizontallyShiftedPatch):  
                            #del self.patches[-1]
                            self.delete_patch(-1)                                                                                          
                            self.n_patches -= 1                        
                        self.makePatch_alignedToLine(apexZ0 = shifted_Align, ppl = ppl, z_top = newZtop, leftRight = (not shiftOriginal))
                        self.patches[-1].getShadows(z_top_min,z_top_max)
                        if shiftOriginal:
                            original_topR_jL = self.patches[-1].shadow_fromTopToInnermost_topR_jL
                            original_topL_jL = self.patches[-1].shadow_fromTopToInnermost_topL_jL
                            original_topR_jR = self.patches[-1].shadow_fromTopToInnermost_topR_jR
                            original_topL_jR = self.patches[-1].shadow_fromTopToInnermost_topL_jR
                        else :
                            complementary_topR_jR = self.patches[-1].shadow_fromTopToInnermost_topR_jR
                            complementary_topL_jR = self.patches[-1].shadow_fromTopToInnermost_topL_jR
                            complementary_topR_jL = self.patches[-1].shadow_fromTopToInnermost_topR_jL
                            complementary_topL_jL = self.patches[-1].shadow_fromTopToInnermost_topL_jL

                        horizontalShiftTop = original_topR_jL - complementary_topR_jR
                        horizontalShiftBottom = original_topL_jL - complementary_topL_jR
                        if (shiftOriginal and self.patches[-1].straightLineProjectorFromLayerIJtoK(original_topR_jR,z_top_max,1,self.env.num_layers,0)<self.env.beam_axis_lim):
                            horizontalOverlapTop = max(complementary_topR_jL - original_topR_jL, complementary_topR_jR - original_topR_jR)
                            horizontalOverlapBottom = max(complementary_topL_jL - original_topL_jL, complementary_topL_jR - original_topL_jR)
                            print('horizontalOverlapTop:',horizontalOverlapTop,'horizontalOverlapBottom:',horizontalOverlapBottom)

                        """    
                        if shiftOriginal:
                            newGapTop = original_topR_jR - originalSaved_topR_jL
                            newGapBottom = original_topL_jR - originalSaved_topL_jL
                        else:
                            newGapTop = complementary_topR_jL - complementarySaved_topR_jR
                            newGapBottom = complementary_topL_jL - complementarySaved_topL_jR
                        print('newGapTop:',newGapTop,'newGapBottom:',newGapBottom)
                        """
                        """
                        if (makeHorizontallyShiftedPatch):
                            if (horizontalOverlapTop > 0) or (horizontalOverlapBottom > 0):
                                del self.patches[-1]
                                self.n_patches -= 1        
                            else:
                                if (len(self.patches) > 1):
                                    del self.patches[-2]
                                    self.n_patches -= 1
                        """
                        print('original_topR_jL:',original_topR_jL,'complementary_topR_jR',complementary_topR_jR,'original_topL_jL',original_topL_jL,'complementary_topL_jR',complementary_topL_jR,'shiftOriginal',shiftOriginal)
                        makeHorizontallyShiftedPatch = True
                        print('updated_horizontalShifts:',horizontalShiftTop,horizontalShiftBottom, 'shifted_Align:',shifted_Align)
                    if (makeHorizontallyShiftedPatch):
                        if ((self.patches[-1].straightLineProjectorFromLayerIJtoK(shifted_Align,newZtop,1,self.env.num_layers,0) > self.env.beam_axis_lim)) and shiftOriginal:
                            if (len(self.patches) > 2):
                                #del self.patches[-3]
                                self.delete_patch(-3)
                                self.n_patches -= 1
                        #if ((self.patches[-1].straightLineProjectorFromLayerIJtoK(shifted_Align,newZtop,1,self.env.num_layers,0) < -self.env.beam_axis_lim)) and (not shiftOriginal):
                            #del self.patches[-2]
                            #self.n_patches -= 1
                            
                z_top_max = c_corner
                print('+++++++++++++++++++++++ c_corner: ', c_corner)

            apexZ0 = self.patches[-1].c_corner[0]
            apexZ0 = saved_apexZ0
            print('=======================================================  z1_Align: ', apexZ0)
        #for i in range(3):
            #del self.patches[-1]                                                                                                  
        #for i in range(34):
            #del self.patches[-1]                                                                                                  
        #for i in range(8):
            #del self.patches[-4]
        
    def makePatches_ShadowQuilt_fromCenter(self, apexZ0 = 0, stop = 1, ppl = 16, leftRight = True):
        """This method uses the geometry of shadows to generate patches based on superpoints
            the outer layer and the z0 of the collision point. 

        Args:
            apexZ0 (num, optional): Collision point on the z axis for the first patches Defaults to 0.
            stop (num, optional): Where to stop, normalized to 1. Defaults to 1.
            ppl (int, optional): Points per patch per layer. Defaults to 16.
            leftRight (bool, optional): If False, goes from right to left instead of left to right. Defaults to True.
        """

        apexZ0 = 0 # switched from z0 to z1, start at center of z1
#      while (apexZ0 > -self.env.trapezoid_edges[0]):
        z_top_min = -self.env.top_layer_lim
        complementary_apexZ0 = 0
        first_row_count = 0
        b_corner = -np.inf
        z_top_max = -self.env.top_layer_lim - self.env.boundaryPoint_offset
        nPatchesInColumn = 0
        projectionOfCcornerToBeam = 0
        while (b_corner < self.env.trapezoid_edges[self.env.num_layers-1]) and (nPatchesInColumn<10):
            nPatchesInColumn += 1
            self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, z_top = z_top_min, leftRight=True)
            print ('top layer from ', self.patches[-1].superpoints[self.env.num_layers-1].max, ' to ', self.patches[-1].superpoints[self.env.num_layers-1].min, ' z_top_min: ', z_top_min)
            print(self.patches[-1].a_corner)
            print(self.patches[-1].b_corner)
            print(self.patches[-1].c_corner)
            print(self.patches[-1].d_corner)
            #original_c = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].c_corner[1], 'below')
            #original_d = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].d_corner[1], 'below')
            original_a = self.patches[-1].a_corner[1]
            original_b = self.patches[-1].b_corner[1]
            a_corner = original_a
            seed_apexZ0 = apexZ0
            projectionOfCcornerToBeam = self.patches[-1].straightLineProjectorFromLayerIJtoK(self.patches[-1].c_corner[1],self.patches[-1].c_corner[0],self.env.num_layers,1,0)
            print('squareAcceptance: ', self.patches[-1].squareAcceptance, 'triangleAcceptance: ', self.patches[-1].triangleAcceptance, ' projectionOfCcornerToBeam: ', projectionOfCcornerToBeam)
            if (self.patches[-1].squareAcceptance == False) and (self.patches[-1].b_corner[1] < self.env.trapezoid_edges[self.env.num_layers-1]):
                complementary_apexZ0 = self.patches[-1].superpoints[0].max
                if (self.patches[-1].triangleAcceptance == True):
                    z_top_max = self.patches[-1].b_corner[1]
                else:     
                    z_top_max = max(z_top_max, self.patches[-1].superpoints[self.env.num_layers-1].max)
                self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_max, leftRight=False)
                print('complementary: ', self.patches[-1].a_corner, ' for z_top_max:', z_top_max)
                print('complementary: ', self.patches[-1].b_corner)
                print('complementary: ', self.patches[-1].c_corner)
                print('complementary: ', self.patches[-1].d_corner)
                #complementary_a = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].a_corner[1], 'above')
                #complementary_b = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].b_corner[1], 'above')
                complementary_c = self.patches[-1].c_corner[1]
                complementary_d = self.patches[-1].d_corner[1]
                white_space_height = max(original_a - complementary_c, original_b - complementary_d)
                previous_white_space_height = -1
                counter = 0
                counterUpshift = 0
                current_z_top_index = -1
                #while ((white_space_height > 0) or (abs(white_space_height) > 10)) and (counter < 5):
                #while (counterUpshift < 100) and (white_space_height != 0) and ((counter < 15) or (white_space_height > 0)) and (self.patches[-1].c_corner[1] > -self.env.trapezoid_edges[self.env.num_layers-1]):
                
                #while not((white_space_height < 0) and (previous_white_space_height >= 0)) and ((self.patches[-1].c_corner[1] > -self.env.trapezoid_edges[self.env.num_layers-1]) or (white_space_height > 0)) and (current_z_top_index < (len(self.data.array[self.env.num_layers-1])-1)) and (self.patches[-2].triangleAcceptance == False) :
                while not((white_space_height <= 0) and (previous_white_space_height >= 0)) and ((self.patches[-1].b_corner[1] < self.env.trapezoid_edges[self.env.num_layers-1]) or (white_space_height > 0)) and (current_z_top_index < (len(self.data.array[self.env.num_layers-1])-1)) :
                    print()
                    if (len(self.patches) > 2):
                        print('original_a:', original_a, ' ', self.patches[-2].a_corner[1], '|| original_b:', original_b, ' ', self.patches[-2].b_corner[1])
                    print('complementary_c:', complementary_c, ' ', self.patches[-1].c_corner[1], ' || complementary_d:', complementary_d, ' ', self.patches[-1].d_corner[1])
                    current_z_top_index = self.get_index_from_z(self.env.num_layers-1, z_top_max)
                    print('current white_space_height: ', white_space_height)
                    print('counter: ',counter)
                    print('orig ztop: ', current_z_top_index)
                    if (white_space_height < 0):
                        current_z_top_index += 1
                        counterUpshift += 1
                    else:
                        counter +=1
                        current_z_top_index -= 1
                    print('new ztop: ', current_z_top_index, ' arrayLength: ', len(self.data.array[self.env.num_layers-1]))
                    current_z_top_index = min(current_z_top_index,len(self.data.array[self.env.num_layers-1])-1)
                    z_top_max = self.data.array[self.env.num_layers-1][current_z_top_index].z
                    #del self.patches[-1]
                    self.delete_patch(-1)
                    self.n_patches -= 1
                    self.makePatch_alignedToLine(apexZ0 = complementary_apexZ0, ppl = ppl, z_top = z_top_max, leftRight=False)
                    #complementary_a = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].a_corner[1], 'above')
                    #complementary_b = self.get_index_from_z(self.env.num_layers-1, self.patches[-1].b_corner[1], 'above')
                    complementary_c = self.patches[-1].c_corner[1]
                    complementary_d = self.patches[-1].d_corner[1]
                    previous_white_space_height = white_space_height
                    white_space_height = max(original_a - complementary_c, original_b - complementary_d)
                    print('complementary_c:', complementary_c, ' ', self.patches[-1].c_corner[1], ' || complementary_d:', complementary_d, ' ', self.patches[-1].d_corner[1])
                    print('new white_space_height: ', white_space_height)
            b_corner = self.patches[-1].b_corner[1]
            z_top_min = b_corner
            print('b_corner: ', b_corner)
#            projectionOfCcornerToBeam = self.patches[-1].straightLineProjectorFromLayerIJtoK(c_corner,self.patches[-1].c_corner[0],self.env.num_layers,1,0)

        #apexZ0 = self.patches[-1].c_corner[0]
        #print('z1_Align: ', apexZ0)
    
        #for i in range(4):
            #del self.patches[-1]                                                                                                  
        #for i in range(34):
            #del self.patches[-1]                                                                                                  
        #for i in range(8):
            #del self.patches[-4]                                                                                                  
        
    def makePatches_ShadowQuilt_fromCenter_v2(self, apexZ0 = 0, stop = 1, z_top = 0, ppl = 16, leftRight = True):
        """This method uses the geometry of shadows to generate patches based on superpoints
            the outer layer and the z0 of the collision point. 

        Args:
            apexZ0 (num, optional): Collision point on the z axis for the first patches Defaults to 0.
            stop (num, optional): Where to stop, normalized to 1. Defaults to 1.
            ppl (int, optional): Points per patch per layer. Defaults to 16.
            leftRight (bool, optional): If False, goes from right to left instead of left to right. Defaults to True.
        """

        first_row_count = 0
        self.makePatches_Projective_center()
        initial_apexZ0 = self.env.beam_axis_lim
        first_row_center  = self.patches
        first_row_center = tuple(first_row_center)
        #print(self.n_patches)
        
        for patch in first_row_center:
            #plt.axvline(patch.a_corner[1], color = 'k')
            #plt.axhline(min(patch.a_corner[0], self.env.top_layer_lim), color = 'r')
            #print(patch.a_corner[0])
            c_corner = patch.c_corner[1]
            bottom_layer_min = 0
            apexZ0 = patch.a_corner[1]
            while (c_corner >= -self.env.beam_axis_lim)  & (bottom_layer_min > -self.env.trapezoid_edges[0]):
                row_z_top = min(patch.right_end_lambdaZ*self.env.radii[-1] + patch.apexZ0, self.env.top_layer_lim)
                row_z_top = min(patch.a_corner[0], self.env.top_layer_lim)
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                             z_top = row_z_top, leftRight=False)
                c_corner = self.patches[-1].c_corner[1]
                bottom_layer_min = self.patches[-1].superpoints[0].min
                apexZ0 = self.patches[-1].a_corner[1]
                #print("going left: ", self.patches[-1].a_corner[1], self.patches[-1].b_corner[1], self.patches[-1].c_corner[1], self.patches[-1].d_corner[1])
                #print("z_top: ", row_z_top)
                for i in range(5):
                    #print(self.patches[-1].superpoints[i].min)
                    pass
                if first_row_count == 0:
                    #plt.axvline(patch.a_corner[1], color = 'k')
                    #plt.axhline(patch.a_corner[0], color = 'r')
                    pass
            
            b_corner = patch.b_corner[1]
            apexZ0 = patch.d_corner[1]
            bottom_layer_max = 0            
            while (b_corner <= self.env.beam_axis_lim) & (bottom_layer_max < self.env.trapezoid_edges[0]):
                row_z_top = max(patch.left_end_lambdaZ*self.env.radii[-1] + patch.apexZ0, -self.env.top_layer_lim)
                row_z_top = max(patch.d_corner[0], -self.env.top_layer_lim)
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                             z_top = row_z_top, leftRight=True)
                b_corner = self.patches[-1].b_corner[1]
                bottom_layer_max = self.patches[-1].superpoints[0].max
                apexZ0 = self.patches[-1].d_corner[-1]
                #print("going right: ", self.patches[-1].a_corner[1], self.patches[-1].b_corner[1], self.patches[-1].c_corner[1], self.patches[-1].d_corner[1])
                #print("z_top: ", row_z_top)
                for i in range(5):
                    #print(self.patches[-1].superpoints[i].max)
                    pass
                if first_row_count == 1:    
                    #plt.axvline(patch.d_corner[1], color = 'k')
                    #plt.axhline(patch.left_end_lambdaZ*self.env.radii[-1] + patch.apexZ0, color = 'r')
                    pass
                #pick a value as next apexZ0 value
                #apexZ0 = self.patches[-1].a_corner[1]
            first_row_count += 1
        #print(self.n_patches)
        
        
    def makePatches_ShadowQuilt_fromCenter_v1(self, apexZ0 = 0, stop = 1, z_top = 0, ppl = 16, leftRight = True):
        """This method uses the geometry of shadows to generate patches based on superpoints
            the outer layer and the z0 of the collision point. 

        Args:
            apexZ0 (num, optional): Collision point on the z axis for the first patches Defaults to 0.
            stop (num, optional): Where to stop, normalized to 1. Defaults to 1.
            ppl (int, optional): Points per patch per layer. Defaults to 16.
            leftRight (bool, optional): If False, goes from right to left instead of left to right. Defaults to True.
        """

        first_row_count = 0
        self.makePatches_Projective_center()
        initial_apexZ0 = self.env.beam_axis_lim
        first_row_center  = self.patches
        first_row_center = tuple(first_row_center)
        #print(self.n_patches)
        
        for patch in first_row_center:
            #plt.axvline(patch.a_corner[1], color = 'k')
            #plt.axhline(min(patch.a_corner[0], self.env.top_layer_lim), color = 'r')
            #print(patch.a_corner[0])
            c_corner = patch.c_corner[1]
            bottom_layer_min = 0
            apexZ0 = patch.a_corner[1]
            while (c_corner >= -self.env.beam_axis_lim)  & (bottom_layer_min > -self.env.trapezoid_edges[0]):
                row_z_top = min(patch.right_end_lambdaZ*self.env.radii[-1] + patch.apexZ0, self.env.top_layer_lim)
                row_z_top = min(patch.a_corner[0], self.env.top_layer_lim)
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, z_top = row_z_top, leftRight=False)
                c_corner = self.patches[-1].c_corner[1]
                bottom_layer_min = self.patches[-1].superpoints[0].min
                apexZ0 = self.patches[-1].a_corner[1]
            #print("going left: ", self.patches[-1].a_corner[1], self.patches[-1].b_corner[1], self.patches[-1].c_corner[1], self.patches[-1].d_corner[1])
            #print("z_top: ", patch.parallelograms[0].top_layer_zmax, " z0: ", patch.parallelograms[0].shadow_topR_jL)
            #self.makePatch_alignedToLine(apexZ0 = patch.parallelograms[0].shadow_topR_jL, ppl = ppl,z_top = patch.parallelograms[0].top_layer_zmax, leftRight=False)
                self.makePatch_alignedToLine(apexZ0 = self.patches[-1].parallelograms[0].shadow_topL_jL, ppl = ppl,z_top = self.patches[-1].parallelograms[0].top_layer_zmin, leftRight=True)
                
                #print("z_top: ", row_z_top)
            #    for i in range(5):
                    #print(self.patches[-1].superpoints[i].min)
            #        pass
            #    if first_row_count == 0:
                    #plt.axvline(patch.a_corner[1], color = 'k')
                    #plt.axhline(patch.a_corner[0], color = 'r')
            #        pass
            
            b_corner = patch.b_corner[1]
            apexZ0 = patch.d_corner[1]
            bottom_layer_max = 0            
            while (b_corner <= self.env.beam_axis_lim) & (bottom_layer_max < self.env.trapezoid_edges[0]):
                row_z_top = max(patch.left_end_lambdaZ*self.env.radii[-1] + patch.apexZ0, -self.env.top_layer_lim)
                row_z_top = max(patch.d_corner[0], -self.env.top_layer_lim)
                self.makePatch_alignedToLine(apexZ0 = apexZ0, ppl = ppl, 
                                             z_top = row_z_top, leftRight=True)
                b_corner = self.patches[-1].b_corner[1]
                bottom_layer_max = self.patches[-1].superpoints[0].max
                apexZ0 = self.patches[-1].d_corner[-1]
                #print("going right: ", self.patches[-1].a_corner[1], self.patches[-1].b_corner[1], self.patches[-1].c_corner[1], self.patches[-1].d_corner[1])
                #print("z_top: ", row_z_top)
                self.makePatch_alignedToLine(apexZ0 = self.patches[-1].parallelograms[0].shadow_topR_jR, ppl = ppl,z_top = self.patches[-1].parallelograms[0].top_layer_zmax, leftRight=False)
                for i in range(5):
                    #print(self.patches[-1].superpoints[i].max)
                    pass
                if first_row_count == 1:    
                    #plt.axvline(patch.d_corner[1], color = 'k')
                    #plt.axhline(patch.left_end_lambdaZ*self.env.radii[-1] + patch.apexZ0, color = 'r')
                    pass
                #pick a value as next apexZ0 value
                #apexZ0 = self.patches[-1].a_corner[1]
            first_row_count += 1
        #print(self.n_patches)
        
    def makePatch_alignedToLine(self, apexZ0 = 0, z_top = -50, ppl = 16, leftRight = True, double_middleLayers_ppl = False):

        init_patch = []
        original_ppl = ppl
        alignmentAccuracy = 0.00001 # 0.1 micron

        #row_data[layer] contains spacepoints for each layer
        row_data = self.data.array
        #loops through each layer and picks n points closest to (z0, 0) and (-100, 25)
        #for row in range(self.env.num_layers-1,-1,-1):
        for row in range(self.env.num_layers):
            y = self.env.radii[row]
            #create compatible arrays from data structure
            row_list = np.array([row_data[row][x].z for x in range(len(row_data[row]))])
            #picks picks n points closest to line from (z0, 0) to (-100, 25) (top left point)
            r_max = self.env.radii[-1]
            #start_index = np.argmin(np.abs((row_list - ((z_top-apexZ0)*y/r_max + apexZ0))))
            #start_value = row_list[start_index] - ((z_top-apexZ0)*y/r_max + apexZ0)
            projectionToRow = (z_top-apexZ0)*(y-self.env.radii[0])/(r_max-self.env.radii[0]) + apexZ0 # apexZ0 is actually defined at layer1
            start_index = np.argmin(np.abs((row_list - projectionToRow)))
            start_value = row_list[start_index] - projectionToRow

            left_bound = np.argmin(np.abs((row_list + self.env.trapezoid_edges[row] + self.env.boundaryPoint_offset)))
            right_bound = np.argmin(np.abs((row_list - self.env.trapezoid_edges[row] - self.env.boundaryPoint_offset)))

            if (double_middleLayers_ppl == True) & (row != 0) & (row!=self.env.num_layers-1):
                ppl = original_ppl * 2 - 1
            else:
                ppl = original_ppl

            if leftRight == True:
                #subtract one from stop index in case it is right of the line from (z0, 0) to (-100, 25)
                if start_index != 0:
                    if start_value > alignmentAccuracy:
                        start_index -= 1
                    pass
                #add superpoint to patch
                if start_index + ppl > right_bound + 1:
                    init_patch.append(wedgeSuperPoint(row_data[row][right_bound+1-ppl:right_bound+1]))
                else:
                    init_patch.append(wedgeSuperPoint(row_data[row][start_index:start_index+ppl]))

            else:
                #add one to stop index in case it is left of the line from (z0, 0) to (100, 25)
                if start_index != len(row_list)-1:
                    print('row',row+1,'start_index',start_index,'start_value',start_value,'z:',row_list[start_index])
                    if start_value < -alignmentAccuracy:
                        start_index += 1
                        start_value = row_list[start_index] - projectionToRow
                        print('row',row+1,'updated start_index',start_index,'start_value',start_value,'z:',row_list[start_index])
                #add superpoint to patch 
                if start_index - ppl + 1 < left_bound:
                    init_patch.append(wedgeSuperPoint(row_data[row][left_bound:left_bound+ppl]))
                else:
                    init_patch.append(wedgeSuperPoint(row_data[row][start_index-ppl+1:start_index+1]))
                    
            #if (row == self.env.num_layers-1):
                # update z_top to the nearest point's coordinate in outermost layer so that other layers' points are accurately aligned
                #z_top = row_list[start_index]                
                
        #init_patch.reverse() # reverse the order of rows to make ascending order
        #add patch to cover
        self.add_patch(wedgePatch(self.env, tuple(init_patch), apexZ0=apexZ0))


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
        '''
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
        '''
        if leftRight == True:
            self.makePatch_alignedToLine(apexZ0=apexZ0, z_top=-self.env.top_layer_lim, ppl = ppl, leftRight=True)
        else:
            self.makePatch_alignedToLine(apexZ0=apexZ0, z_top=self.env.top_layer_lim, ppl = ppl, leftRight=False)
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
            if (self.data.array[row][center_index].z >= 0) or (center_index+16 == len(row_list)):
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
                #del self.patches[n_patch_start-1]
                self.delete_patch(n_patch_start-1)
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
                #del self.patches[n_patch_start-1]
                self.delete_patch(n_patch_start-1)
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
            #del self.patches[n_patch_start-1]
            self.delete_patch(n_patch_start-1)
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
                    plt.savefig(f"python/temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == False and data == True: 
                name = 0
                for patch in self.patches: 
                    self.data.plot(show = False) 
                    patch.plot("b")
                    plt.savefig(f"python/temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == True and data == False: 
                name = 0
                for patch in self.patches: 
                    for line in self.fitting_lines: 
                        if patch.contains(line): 
                            line.plot("r")
                    patch.plot("b")
                
                    plt.savefig(f"python/temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
                    
            elif lines == True and data == True: 
                name = 0
                for patch in self.patches: 
                    self.data.plot(show = False, show_lines=True) 
                    for line in self.fitting_lines: 
                        if patch.contains(line): 
                            line.plot("r")
                    patch.plot("b")
                    plt.savefig(f"python/temp_image_dir/{str(name).zfill(2)}.png")
                    plt.clf() 
                    name += 1 
   
                
                        
            image_files = glob.glob("python/temp_image_dir/*.png")
            
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

    def movie(self, z0_spacing=0.5, figSizeScale = 6.5):
        zInnerLayer = np.arange(-22, 22+z0_spacing, z0_spacing)
        
        if self.n_patches == 0: 
            raise Exception("You must run the solve method first to generate the patches. ")
        
        name = 0
        plt.figure(figsize = (1.7*self.env.beam_axis_lim/figSizeScale, self.env.top_layer_lim/figSizeScale))
        for i, patch in enumerate(self.all_patches): 
            for iz, zIn in enumerate(np.array(zInnerLayer)):
                
                list_of_intersections = []
                for j, patch in enumerate(self.all_patches[:i+1]): 
                    if (j == i) or (self.real_patch_list[j]):
                        list_of_segs = [pgram.crossSection(zIn) for pgram in patch.parallelograms]
                        overlap_of_superpoints = intersection(patch.env, list_of_segs, True) 
                        list_of_intersections.append(overlap_of_superpoints)

                
                plt.xlabel(r"$z_1$ (cm)", fontsize = 18)
                plt.ylabel(r"$z_{top}$ (cm)", fontsize = 18)
                plt.title("acceptance of cover", fontsize = 18)
                z1Lim = self.all_patches[-1].straightLineProjectorFromLayerIJtoK(-self.env.top_layer_lim ,self.env.beam_axis_lim,self.env.num_layers,0,1)
                plt.axline((z1Lim, -self.env.top_layer_lim), (self.env.trapezoid_edges[0], self.env.top_layer_lim),linewidth=1, color='black')
                plt.axline((-z1Lim, self.env.top_layer_lim), (-self.env.trapezoid_edges[0], -self.env.top_layer_lim),linewidth=1, color='black')

                colors = ["b", "r", "g", "c", "m", "y", "k", "chocolate", "indigo", "springgreen", "orange", "rosybrown", "tomato","olive", "deeppink"]
                    
                col = 0
                for line in list_of_intersections: 
                    #plt.xlim(-z0_luminousRegion,z0_luminousRegion)
                    plt.xlim(-self.env.trapezoid_edges[0],self.env.trapezoid_edges[0])
                    plt.ylim(-self.env.top_layer_lim, self.env.top_layer_lim)
                    plt.plot([zIn, zIn], [line.min_z5_accepted, line.max_z5_accepted], c=colors[col % len(colors)], alpha=0.3, linewidth=3)
                    col += 1

                #plt.show()
            plt.savefig(f"python/temp_image_dir/{str(name).zfill(2)}.png")
            plt.clf() 
            name += 1 

                
                        
        image_files = glob.glob("python/temp_image_dir/*.png")
        
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