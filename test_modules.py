import numpy as np 
import matplotlib.pyplot as plt 
from data import * 
from cover import * 
import time 

# clustering - "LeftRight", "Center"
# linings - "LeftRight", "CenterGrid", "CenterSpread", "Randomized"


def numCovers(clustering:str = "", lining:str = "solveS", events=1000, savefig=False, ideal=False): 
    # Runs a bunch of iterations by generating 1000 datasets and 
    # computing the cover. Then, it just looks at how many covers is
    # being generated for each dataset. The lower the distribution the better. 
    num_covers = [] 
    for _ in range(events): 
        
        env = Environment()
        if ideal == True:
            data = DataSet(env, n_points=150, equal_spacing = True)
        else:
            data = DataSet(env, n_points=150)

        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, nlines=100)
        
        num_covers.append(cover.n_patches)
    avg = np.mean(num_covers)
    std = np.std(num_covers)
    plt.hist(num_covers, 
                bins=np.arange(np.min(num_covers), np.max(num_covers)+2)-0.5,
                edgecolor='black', 
                rwidth=0.8, label = f"mean: {format(avg, '.2f')}, stdev: {format(std, '.2f')}"
            )
    print(f"({clustering}, {lining}) - {format(avg, '.2f')}, {format(std, '.2f')}")
    plt.title(f"Number of Patches per Cover ({clustering}, {lining})", fontsize = '20')
    plt.xlabel("Number of Patches", fontsize = '16')
    plt.ylabel("Number of Covers", fontsize = '16')
    plt.legend(fontsize = '20')
    if savefig == True: 
        if ideal == True:
            plt.savefig(f"Figures/nPatches_({clustering}_{lining}_ideal)")
        else:
            plt.savefig(f"Figures/nPatches_({clustering}_{lining})")
    plt.show() 

def fourTests(clustering:str = "", lining:str = "solveS", solve_at = 0, z0 = 0, events=100, lines=1000, savefig=False, ideal=False):
    mean_list = []
    num_covers = [] 
    for z in np.array(z0):
        percentage_accepted = [0 for _ in range(lines)] 
        
        for k in range(events): 
            env = Environment()
            if ideal == True:
                data = DataSet(env, n_points=150, equal_spacing = True)
            else:
                data = DataSet(env, n_points=150)
            cover = Cover(env, data)
            cover.solve(clustering=clustering, z0 = solve_at, lining=lining, show = False)
            num_covers.append(cover.n_patches)
            
            lg = LineGenerator(env, z)
            test_lines = lg.generateEvenGrid(lines)
            
            
            for i in range(len(test_lines)): 
                for patch in cover.patches: 
                    if patch.contains(test_lines[i]): 
                        percentage_accepted[i] += 1 
                        break 
                
        percentage_accepted = [x / events for x in percentage_accepted]
        mean_list.append(np.mean(percentage_accepted))
        mean_accept = format(np.mean(percentage_accepted), ".3f")
        mean_num = format(np.mean(num_covers), ".1f")
        std_num = format(np.std(num_covers), ".1f")
        print(f"({clustering}, {lining}) - {mean_accept}")
    plt.scatter(z0, mean_list, color = 'r')
    plt.plot(z0, mean_list, color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance Rate',  fontsize = 16)
    plt.ylim(-0.1, 1.1)
    plt.title(f'z0 Offset vs Acceptance Rate for {lining}', fontsize = 16)
    PRF = pointRepetitionFactor(lining = lining, ideal = ideal, z0 = solve_at, show = False)
    plt.legend([f"""Number of Patches: {mean_num}+-{std_num}\nPoint Reptition Factor: {PRF[0]}+-{PRF[1]}\nPatches generated at z0 = {solve_at}"""],
        loc = 8, fontsize = 12)
    if savefig == True: 
        if ideal == True:
            plt.savefig(f"Figures/Accept_vs_z0_({lining}_ideal)")
        else:
            plt.savefig(f"Figures/Accept_vs_z0_({lining}_10)")
    plt.show()
    
    
def acceptSlopePlot(clustering:str = "", lining:str = "solveS", events=100, lines=1000, savefig=False, ideal=False):
    
    percentage_accepted = [0 for _ in range(lines)] 
    
    
    for k in range(events): 
        env = Environment()
        if ideal == True:
            data = DataSet(env, n_points=150, equal_spacing = True)
        else:
            data = DataSet(env, n_points=150)
        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, nlines=100)
        
        lg = LineGenerator(env, 0.0)
        test_lines = lg.generateEvenGrid(lines)
        co_tan = []
        
        
        for i in range(len(test_lines)): 
            color = "r"
            co_tan.append(100/test_lines[i].slope)
            for patch in cover.patches: 
                if patch.contains(test_lines[i]): 
                    color="g"
                    percentage_accepted[i] += 1 
                    break 
            
            if color == "r": 
                
                # print(i)
                pass

    percentage_accepted = [x / events for x in percentage_accepted]
    mean_accept = format(np.mean(percentage_accepted), ".3f")
    plt.plot(co_tan, percentage_accepted, c="b", label = "Mean acceptance: "+mean_accept)
    print(f"({clustering}, {lining}) - {mean_accept}")
    
    plt.title(f"Acceptance Rate ({clustering}, {lining})", fontsize = '20')
    plt.xlabel("dZ/dr", fontsize = '16')
    plt.ylabel("Acceptance Probability", fontsize = '16')
    plt.legend(fontsize = '16')
    if savefig == True: 
        if ideal == True:
            plt.savefig(f"Figures/Acceptance_Rate_({clustering}_{lining}_ideal)")
        else:
            plt.savefig(f"Figures/Acceptance_Rate_({clustering}_{lining})")
    plt.show() 
            
def pointRepetitionFactor(clustering:str = "", lining:str = "solveS", z0 = 0, events=10, savefig=False, ideal=False, show = True): 
    # for every event, we loop through all the points in the dataset and compute 
    # how many patches contain that point. The lower in general the better, since 
    # this is a metric of non-wastefulness 

    out = [] 
    
    for _ in range(events): 
        
        env = Environment()
        if ideal == True:
            data = DataSet(env, n_points=150, equal_spacing = True)
        else:
            data = DataSet(env, n_points=150)
        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, z0 = z0, nlines=100, show = False)
        
        out2 = [] 
        
        for layer in range(env.layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                for patch in cover.patches: 
                    if patch.contains_p(point, layer): 
                        num_in += 1
                        
                out2.append(num_in) 
                
        out += out2

        # print(out2)
        
    # plt.scatter(*zip(*unaccept)) 
    
    # plt.show()
    if show == False:
        return (format(np.mean(out), '.2f'), format(np.std(out), '.2f'))

    print(f"({clustering}, {lining}) mean - {format(np.mean(out), '.2f')}")
    print(f"({clustering}, {lining}) stdev - {format(np.std(out), '.2f')}")


        
    plt.hist(out, bins=np.arange(11) - 0.5, 
             edgecolor='black', 
             label = f"mean: {format(np.mean(out), '.2f')}, stdev: {format(np.std(out), '.2f')}",
             rwidth=0.8
            ) 
    plt.xlabel("Number of Covering Patches", fontsize = '16')
    plt.ylabel("Number of Points", fontsize = '16')
    plt.title(f"Point Repetition Factor ({clustering}, {lining})", fontsize = '20')
    plt.legend(fontsize = '16')
    if savefig == True:
        if ideal == True:
            plt.savefig(f"Figures/Point_Repetition_Factor_({clustering}_{lining}_ideal)")
        else:
            plt.savefig(f"Figures/Point_Repetition_Factor_({clustering}_{lining})")
    plt.show()
    return (format(np.mean(out), '.2f'), format(np.std(out), '.2f'))
    
def idealData(clustering:str = "", lining:str = "solveS"): 
    env = Environment()
    data = DataSet(env, n_points=150, equal_spacing=True) 
    cover = Cover(env, data) 
    cover.solve(clustering=clustering, lining=lining, nlines=100)
    
    print(f"Figures/Number of Patches: {cover.n_patches}")
    cover.plot() 


def duplicates(lining:str = "solveS", z0 = 0, events=1000, ideal=False):
    dupes = []
    env = Environment()

    for times in range(events):
        data = DataSet(env, n_points=150)
        cover = Cover(env, data)
        cover.solve(lining=lining, nlines=100, z0 = z0, show = False)
        cover_array = np.zeros((cover.n_patches, 5, 16))
        d = 0
        for i, n in enumerate(cover.patches):
            for j, p in enumerate(n.superpoints):
                cover_array[i, j, :] = np.array(p.points)

        for i in range(len(cover_array)):
            index = len(cover_array)-i
            for j in range(index-1):
                if np.all(cover_array[i] == cover_array[i+j+1]):
                    d += 1
                    print(cover_array[i])
                    print(len(cover_array))
                    print(i)
                    print(i+j+1)
                    cover.plot()
        dupes.append(d)

    print(f'{lining} - mean: {np.mean(dupes)} std: {np.std(dupes)}')
