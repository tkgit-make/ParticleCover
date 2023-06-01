import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import * 
import time 

# cluster_type - "LeftRight", "Center"
# lines_types - "LeftRight", "CenterGrid", "CenterSpread", "Randomized"


def numCovers(cluster_type:str, lines_type:str, events=1000): 
    # Runs a bunch of iterations by generating 1000 datasets and 
    # computing the cover. Then, it just looks at how many covers is
    # being generated for each dataset. The lower the distribution the better. 
    num_covers = [] 
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        
        if cluster_type == "LeftRight": 
            cover.cluster("LR") 
        elif cluster_type == "Center": 
            cover.cluster("C")
        else: 
            raise Exception("This is not a valid cluster type.")
        
        if lines_type == "LeftRight": 
            cover.solveGridLR(100) 
        elif lines_type == "Randomized": 
            cover.solveRandomized(100) 
        elif lines_type == "CenterSpread": 
            cover.solveCenterSpread(100)
        elif lines_type == "CenterGrid": 
            cover.solveCenterGrid(100) 
        else: 
            raise Exception("This is not a valid lines type.")
        
        num_covers.append(cover.n_patches)

    plt.hist(num_covers, 
                bins=np.arange(min(num_covers), max(num_covers)) - 0.5, 
                edgecolor='black', 
                rwidth=0.8
            )
    avg = np.mean(num_covers) 
    std = np.std(num_covers)
    print(f"({cluster_type}, {lines_type}) - {format(avg, '.2f')}, {format(std, '.2f')}")
    plt.title(f"Number of Patches per Cover ({cluster_type}, {lines_type})")
    plt.xlabel("Number of Patches")
    plt.ylabel("Number of Covers")
    plt.savefig(f"nPatches_({cluster_type}_{lines_type})")
    plt.show() 
    
def acceptSlopePlot(cluster_type:str, lines_type:str, events=100, lines=1000):
    
    percentage_accepted = [0 for _ in range(lines)] 
    
    
    for k in range(events): 
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        
        if cluster_type == "LeftRight": 
            cover.cluster("LR") 
        elif cluster_type == "Center": 
            cover.cluster("C")
        else: 
            raise Exception("This is not a valid cluster type.")
        
        if lines_type == "LeftRight": 
            cover.solveGridLR(100) 
        elif lines_type == "Randomized": 
            cover.solveRandomized(100) 
        elif lines_type == "CenterSpread": 
            cover.solveCenterSpread(100)
        elif lines_type == "CenterGrid": 
            cover.solveCenterGrid(100) 
        else: 
            raise Exception("This is not a valid lines type.")
        
        # for layer in cover.superPoints: 
        #     for sp in layer: 
        #         print(sp.min, sp.max)
        #     print("------")
        
        lg = LineGenerator(env, 0.0)
        test_lines = lg.generateGridLines(lines) 
        
        
        for i in range(len(test_lines)): 
            color = "r"
            for patch in cover.patches: 
                if patch.contains(test_lines[i]): 
                    color="g"
                    percentage_accepted[i] += 1 
                    break 
            
            if color == "r": 
                
                # print(i)
                pass

    percentage_accepted = [x / 100 for x in percentage_accepted]
    plt.plot(np.arange(1000), percentage_accepted, c="b")
    mean_accept = format(np.mean(percentage_accepted), ".3f")
    print(f"({cluster_type}, {lines_type}) - {mean_accept}")
    
    plt.title(f"Acceptance Rate ({cluster_type}, {lines_type})")
    plt.xlabel("Slope (Indexed from Left to Right)")
    plt.ylabel("Acceptance Probability")
    plt.savefig(f"Acceptance_Rate_({cluster_type}_{lines_type})")
    plt.show() 
            
def pointRepetitionFactor(cluster_type:str, lines_type:str, events=10): 
    # for every event, we loop through all the points in the dataset and compute 
    # how many patches contain that point. The lower in general the better, since 
    # this is a metric of non-wastefulness 

    out = [] 
    
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        
        if cluster_type == "LeftRight": 
            cover.cluster("LR") 
        elif cluster_type == "Center": 
            cover.cluster("C")
        else: 
            raise Exception("This is not a valid cluster type.")
        
        if lines_type == "LeftRight": 
            cover.solveGridLR(100) 
        elif lines_type == "Randomized": 
            cover.solveRandomized(100) 
        elif lines_type == "CenterSpread": 
            cover.solveCenterSpread(100)
        elif lines_type == "CenterGrid": 
            cover.solveCenterGrid(100) 
        else: 
            raise Exception("This is not a valid lines type.")
        
        
        out2 = [] 
        unaccept = []
        
        for layer in range(env.layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                
                for patch in cover.patches: 
                    if patch.contains_p(point, layer): 
                        num_in += 1
                
                # if num_in == 0: 
                #     plt.scatter(point, layer+1, s=2, c="r")
                #     unaccept.append((point, layer+1))
                # else: 
                #     plt.scatter(point, layer+1, s=2, c="b")
                        
                out2.append(num_in) 
                
        out += out2
        # print(out2)
        
    # plt.scatter(*zip(*unaccept)) 
    
    # plt.show() 
    print(f"({cluster_type}, {lines_type}) - {format(np.mean(out), '.2f')}")
        
    plt.hist(out, bins=np.arange(11) - 0.5, 
             edgecolor='black', 
             rwidth=0.8
            ) 
    plt.xlabel("Number of Covering Patches")
    plt.ylabel("Number of Points")
    
    plt.title(f"Point Repetition Factor ({cluster_type}, {lines_type})")
    plt.savefig(f"Muchang/Point_Repetition_Factor_({cluster_type}_{lines_type})")
    plt.show() 
    
    
    
numCovers("LeftRight", "LeftRight")
numCovers("LeftRight", "Randomized")
numCovers("LeftRight", "CenterGrid")
numCovers("LeftRight", "CenterSpread")
numCovers("Center", "LeftRight")
numCovers("Center", "Randomized")
numCovers("Center", "CenterGrid",)
numCovers("Center", "CenterSpread")

# acceptSlopePlot("LeftRight")
# acceptSlopePlot("Randomized") 
# acceptSlopePlot("CenterSpread")
# acceptSlopePlot("CenterGrid")

# numCovers("CenterGrid")
# numCovers("LeftRight", "CenterGrid")
# numCovers("Center", "CenterGrid")




# pointRepetitionFactor(100)
# numCovers("C")
