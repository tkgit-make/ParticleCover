import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import * 
import time 

def numCovers(events=1000): 
    # Runs a bunch of iterations by generating 1000 datasets and 
    # computing the cover. Then, it just looks at how many covers is
    # being generated for each dataset. The lower the distribution the better. 
    num_covers = [] 
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.solve(100) 
        num_covers.append(cover.n_patches)

    plt.hist(num_covers)
    plt.show() 
    
def acceptSlopePlot(): 
    
    pass 

def pointRepetitionFactor(events=10): 
    # for every event, we loop through all the points in the dataset and compute 
    # how many patches contain that point. The lower in general the better, since 
    # this is a metric of non-wastefulness 

    out = [] 
    
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.solve() 
        
        out2 = [] 
        
        for layer in range(env.layers): 
            for point in data.array[layer]: 
                num_in = 0
                for patch in cover.patches: 
                    if patch.contains_p(point, layer+1): 
                        num_in += 1
                        
                out2.append(num_in) 
                
        out.append(out2) 
        
    hi = 1
    for plot in out: 
        plt.hist(plot) 
        plt.savefig(f"{hi}.png")
        # plt.show() 
        plt.clf()
        hi += 1

pointRepetitionFactor()

# env = Environment()
# data = DataSet(env, n_points=150) 
# cover = Cover(env, data) 
# cover.solve(100) 

# num = 1
# for patch in cover.patches: 
#     patch.plot() 
#     plt.savefig(f"Muchang/images3/{num}.png")
#     plt.clf()
#     num += 1