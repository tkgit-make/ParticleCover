from src.coverers.data_structs import Environment
from src.debug import * 

class lineSegment(): 
    
    def __init__(self, min_z5_accepted:float, max_z5_accepted:float): 
        
        # the min and max z5 values that are accepted by 1 parallelogram
        
        if min_z5_accepted > max_z5_accepted: 
            raise Exception("The minimum z5 accepted is greater than the maximum z5 accepted.")
        
        self.min_z5_accepted = min_z5_accepted
        self.max_z5_accepted = max_z5_accepted
        
    
def intersection(env:Environment, lineSegments:list): 
    
    if len(lineSegments) != env.num_layers - 1: 
        raise Exception("The number of line segments does not match the number of layers in the environment minus 1")
    
    all_minimums = [lineSeg.min_z5_accepted for lineSeg in lineSegments]
    all_maximums = [lineSeg.max_z5_accepted for lineSeg in lineSegments]
    
    max_of_mins = max(all_minimums)
    min_of_maxes = min(all_maximums)
    
    if max_of_mins > min_of_maxes: 
        # this case returns a line segment of 0 length 
        # if the intersection is null
        return lineSegment(max_of_mins, max_of_mins)
    else: 
        return lineSegment(max_of_mins, min_of_maxes)


def insideFunction(sorted_lineSegments:list, total_measure:float=0.0): 
    
    if len(sorted_lineSegments) != 1: 

        first = sorted_lineSegments[0]
        second = sorted_lineSegments[1]
        
        if second.min_z5_accepted > first.max_z5_accepted: 
            total_measure += first.max_z5_accepted - first.min_z5_accepted 
        else: 
            union = lineSegment(first.min_z5_accepted, max(first.max_z5_accepted, second.max_z5_accepted))
            sorted_lineSegments[1] = union
            
        del sorted_lineSegments[0]
            
        return insideFunction(sorted_lineSegments, total_measure)
    
    else: 
        total_measure += sorted_lineSegments[0].max_z5_accepted - sorted_lineSegments[0].min_z5_accepted 
        return total_measure 

def unionOfLineSegments(lineSegments:list): 

    total_measure = 0.0 
    
    lineSegments.sort(key=lambda x : x.min_z5_accepted)
    total_measure = insideFunction(lineSegments)
    
    return total_measure 
            

class parallelogram(): 
    
    def __init__(self, layer_num, top_layer_zmin, top_layer_zmax, shadow_topR_jL, shadow_topR_jR, pSlope): 
        
        self.layer_num = layer_num
        self.pSlope = pSlope

        # shadow is the top base (at z5) of the parallelogram

        self.shadow_topR_jL = shadow_topR_jL    # a
        self.shadow_topR_jR = shadow_topR_jR    # b
        
        delta_ztop = top_layer_zmax - top_layer_zmin
        delta_z0 = delta_ztop * pSlope 
        
        
        self.shadow_topL_jL = shadow_topR_jL - delta_z0 # c
        self.shadow_topL_jR = shadow_topR_jR - delta_z0 # d
        
        self.top_layer_zmin = top_layer_zmin
        self.top_layer_zmax = top_layer_zmax
        
        if self.top_layer_zmax > 100 or self.top_layer_zmin < -100: 
            print(self.top_layer_zmin, self.top_layer_zmax)
        
        
    def crossSection(self, z0:float): 
        # vertical cross section of the parallelogram at particular z0
        # returns a line interval 
        
        if z0 <= self.shadow_topR_jL or z0 >= self.shadow_topL_jR: 
            # left of a or right of d
            return lineSegment(0.0, 0.0)
        
        if self.shadow_topR_jL <= z0 <= self.shadow_topR_jR: 
            segment_max = self.top_layer_zmax
        elif self.shadow_topR_jR < z0 < self.shadow_topL_jR: 
            segment_max = self.top_layer_zmax + (z0 - self.shadow_topR_jR)/self.pSlope

        if self.shadow_topL_jL <= z0 <= self.shadow_topL_jR: 
            segment_min = self.top_layer_zmin 
        elif self.shadow_topR_jL < z0 < self.shadow_topL_jL: 
            segment_min = self.top_layer_zmax + (z0 - self.shadow_topR_jL)/self.pSlope
        
        return lineSegment(segment_min, segment_max)
    
def calc_line_intersection(line1_point:tuple, line1_slope, line2_point:tuple, line2_slope):

    if line1_slope == line2_slope:
        raise Exception("The slope is the same, so it is either the same line or lines do not intersect")
    
    z_top = (line2_point[1]-line1_point[1] - (line2_slope*line2_point[0]) + (line1_slope*line1_point[0]))/(line1_slope-line2_slope)
    z_0 = line1_slope*(z_top-line1_point[0]) + line1_point[1]
    #z_top is the indepedent variable and z_0 is the dependent variable
    return tuple((z_top, z_0))