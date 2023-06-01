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
        cover.michelle() 
        num_covers.append(cover.n_patches)

    plt.hist(num_covers, 
                bins=np.arange(10, 25) - 0.5, 
                edgecolor='black', 
                rwidth=0.8
            )
    plt.title("Number of Patches per Dataset")
    plt.show() 
    
def acceptSlopePlot(events=100, lines=1000):
    
    percentage_accepted = [0 for _ in range(lines)] 
    
    
    for k in range(events): 
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.michelle()  
        
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
            #test_lines[i].plot(color = color)

    percentage_accepted = [x / 100 for x in percentage_accepted]
    plt.plot(np.arange(1000), percentage_accepted, c="b")
    print(np.mean(percentage_accepted))
    plt.title("Acceptance vs Slope")
    plt.show() 
                

def pointRepetitionFactor(events=10): 
    # for every event, we loop through all the points in the dataset and compute 
    # how many patches contain that point. The lower in general the better, since 
    # this is a metric of non-wastefulness 

    out = [] 
    
    for _ in range(events): 
        
        env = Environment()
        data = DataSet(env, n_points=150) 
        cover = Cover(env, data) 
        cover.michelle() 
        
        out2 = [] 
        unaccept = []
        
        for layer in range(env.layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                
                for patch in cover.patches: 
                    if patch.contains_p(point, layer): 
                        num_in += 1
                
                if num_in == 0: 
                    plt.scatter(point, layer+1, s=2, c="r")
                    # unaccept.append((point, layer+1))
                else: 
                    plt.scatter(point, layer+1, s=2, c="b")
                        
                out2.append(num_in) 
                
        out += out2
        print(out2)
        
    #plt.scatter(*zip(*unaccept)) 
    
    #plt.show() 
        
    plt.hist(out, bins=np.arange(11) - 0.5, 
             edgecolor='black', 
             rwidth=0.8
            ) 
    
    #print(len(out)) 
        
    plt.show() 
    

env = Environment()
data = DataSet(env, n_points=150) 
cover = Cover(env, data) 
cover.michelle() 

for patch in cover.patches: 
    patch.plot()
    plt.show() 
    plt.clf() 

# acceptSlopePlot()
#pointRepetitionFactor(100)
# numCovers()
