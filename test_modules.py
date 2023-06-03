import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import * 
import time 

# clustering - "LeftRight", "Center"
# linings - "LeftRight", "CenterGrid", "CenterSpread", "Randomized"


def numCovers(clustering:str, lining:str, events=1000, savefig=False): 
    # Runs a bunch of iterations by generating 1000 datasets and 
    # computing the cover. Then, it just looks at how many covers is
    # being generated for each dataset. The lower the distribution the better. 
    num_covers = [] 
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, nlines=100)
        
        num_covers.append(cover.n_patches)

    plt.hist(num_covers, 
                bins=np.arange(min(num_covers), max(num_covers)) - 0.5, 
                edgecolor='black', 
                rwidth=0.8
            )
    avg = np.mean(num_covers) 
    std = np.std(num_covers)
    print(f"({clustering}, {lining}) - {format(avg, '.2f')}, {format(std, '.2f')}")
    plt.title(f"Number of Patches per Cover ({clustering}, {lining})")
    plt.xlabel("Number of Patches")
    plt.ylabel("Number of Covers")
    if savefig == True: 
        plt.savefig(f"nPatches_({clustering}_{lining})")
    plt.show() 
    
def acceptSlopePlot(clustering:str, lining:str, events=100, lines=1000, savefig=False):
    
    percentage_accepted = [0 for _ in range(lines)] 
    
    
    for k in range(events): 
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, nlines=100)
        
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
    print(f"({clustering}, {lining}) - {mean_accept}")
    
    plt.title(f"Acceptance Rate ({clustering}, {lining})")
    plt.xlabel("Slope (Indexed from Left to Right)")
    plt.ylabel("Acceptance Probability")
    if savefig == True: 
        plt.savefig(f"Acceptance_Rate_({clustering}_{lining})")
    plt.show() 
            
def pointRepetitionFactor(clustering:str, lining:str, events=10, savefig=False): 
    # for every event, we loop through all the points in the dataset and compute 
    # how many patches contain that point. The lower in general the better, since 
    # this is a metric of non-wastefulness 

    out = [] 
    
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, nlines=100)
        
        out2 = [] 
        
        for layer in range(env.layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                for patch in cover.patches: 
                    if patch.contains_p(point, layer): 
                        num_in += 1
                        
                out2.append(num_in) 
                
        out += out2

    print(f"({clustering}, {lining}) - {format(np.mean(out), '.2f')}")
        
    plt.hist(out, bins=np.arange(11) - 0.5, 
             edgecolor='black', 
             rwidth=0.8
            ) 
    plt.xlabel("Number of Covering Patches")
    plt.ylabel("Number of Points")
    plt.title(f"Point Repetition Factor ({clustering}, {lining})")
    if savefig == True: 
        plt.savefig(f"Muchang/Point_Repetition_Factor_({clustering}_{lining})")
    plt.show() 
    
def idealData(clustering:str, lining:str): 
    env = Environment()
    data = DataSet(env, n_points=150, equal_spacing=True) 
    cover = Cover(env, data) 
    cover.solve(clustering=clustering, lining=lining, nlines=100)
    
    print(f"Number of Patches: {cover.n_patches}")
    cover.plot() 


idealData(clustering="", lining="SlopeCenterStack1")
# idealData(clustering="LeftRight", lining="LeftRight")