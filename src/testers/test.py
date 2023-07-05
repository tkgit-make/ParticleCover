import numpy as np 
import matplotlib.pyplot as plt 
import time, ast

from coverers.data_structs import * 
from coverers.cover import * 
from src.readers.reader import *
from src.coverers.wedgecover import *
from src.deprecated.reader import * 

def wedge_test(lining:str = "solveS", solve_at = 0, z0 = np.arange(-15, 15.5, 0.5), n = 16, wedges = [0, 128], lines=1000, v = 'v3', savefig=False):
    """Creates acceptance vs z0 plot
    
    Args:
        lining (str, optional): solving method, default is solveS
        solve_at (int, optional): the z values the patches are being made in accordance to
        z0 (num or list, optional): array of z0 values we are testing over, default is range of -15 to 15 with 0.5 spacing
        n (int, optional): ppl (point per patch per layer), default is 16
        wedges (list, optional): which wedges, enter in list of [starting wedge, ending wedge]
        lines (int, optional): how many line to test acceptance with at each z0 value
        savefig (bool, optional): True to save figure
        v (str, optional): version of data, ensure data file is in directory as "wedgeData_{v}_128.txt"
    """

    #create list for z values we're testing the acceptance of, number of covers, and PRF
    mean_list = np.zeros((len(z0), wedges[1]-wedges[0]))
    num_covers = []
    PRF = []

    #read wedgeData file and create environment
    all_data = readFile(f'wedgeData_{v}_128.txt', wedges[1])
    env = Environment()

    #loop through all events
    for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
        #convert to existing data format
        data = convertToDataset(all_data[k])
        #add the 0.1 cm points
        data.add()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(z0 = solve_at, lining=lining, n = n, show = False)
        #append number of covers in the patch
        num_covers.append(cover.n_patches)
        out = [] 

        #these loops calculate PRF
        for layer in range(env.layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                for patch in cover.patches: 
                    if patch.contains_p(point, layer): 
                        num_in += 1
                        
                out.append(num_in)
        PRF.append(out)

        #these loops calculate line acceptance
        for iz, z in enumerate(np.array(z0)):
            percentage_accepted = 0 
                
            lg = LineGenerator(env, z)
            test_lines = lg.generateEvenGrid(lines)
            
            for i in range(len(test_lines)): 
                for patch in cover.patches:
                    if patch.contains(test_lines[i]): 
                        percentage_accepted += 1 
                        break 

            
            percentage_accepted = percentage_accepted/lines
            mean_list[iz, ik] = mean_list[iz, ik] + percentage_accepted
    
    mean_accept = format(np.mean(mean_list), ".3f")
    mean_num = format(np.mean(num_covers), ".1f")
    std_num = format(np.std(num_covers), ".1f")

    #sets minimum for plot
    if type(solve_at) == float:
        ymin = 0
    elif type(solve_at) == int:
        ymin = 0
    else:
        ymin = 0.9

    #creates plots and saves them
    plt.scatter(z0, np.mean(mean_list, axis = 1), color = 'r')
    plt.plot(z0, np.mean(mean_list, axis = 1), color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance',  fontsize = 16)
    plt.ylim(ymin, 1.0)
    plt.title(f'{lining}', fontsize = 16)
    PRFm = format(np.mean(out), '.2f')
    PRFs = format(np.std(out), '.2f')
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\nPatches with " + r'$z_0$' + f" = {solve_at}\nppl = {n}, " + r'$N_{wedges}$ ' + f"= {wedges[1]}, {v} events"],
        loc = 8, fontsize = 12)
    if savefig == True:
        try:
            at = len(solve_at)
        except:
            at = 0
        plt.savefig(f"Figures/wedge_test({lining}_{v.replace('.','')}_{at}_n{n})")
    plt.show()

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

def fourTests(clustering:str = "", lining:str = "solveS", solve_at = 0, z0 = 0, n = 16, events=100, lines=1000, savefig=False, ideal=False):
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
            cover.solve(clustering=clustering, z0 = solve_at, lining=lining, n = n, show = False)
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
    plt.ylim(0, 1.0)
    plt.title(f'z0 Offset vs Acceptance Rate for {lining}', fontsize = 16)
    PRF = pointRepetitionFactor(lining = lining, ideal = ideal, z0 = solve_at, n = n, show = False)
    plt.legend([f"""Number of Patches: {mean_num}+-{std_num}\nPoint Reptition Factor: {PRF[0]}+-{PRF[1]}\nPatches generated at z0 = {solve_at}"""],
        loc = 8, fontsize = 12)
    if savefig == True: 
        if ideal == True:
            plt.savefig(f"Figures/Accept_vs_z0_({lining}_ideal)")
        else:
            plt.savefig(f"Figures/Accept_vs_z0_({lining}_10_n{n})")
    plt.show()    
    
def acceptSlopePlot(clustering:str = "", lining:str = "solveS", events=100, lines=1000, z0= 0, solve_at=0, savefig=False, ideal=False, show = True, custom=''):
    
    percentage_accepted = [0 for _ in range(lines)] 
    
    
    for k in range(events): 
        env = Environment()
        if ideal == True:
            data = DataSet(env, n_points=150, equal_spacing = True)
        else:
            data = DataSet(env, n_points=150)
        if '' not in custom:
            data.input_data(custom)
        cover = Cover(env, data) 
        cover.solve(clustering=clustering, lining=lining, z0 = solve_at, show = False)
        
        lg = LineGenerator(env, z0)
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
    
    if show == False:
        return np.mean(percentage_accepted)

    print(f"({clustering}, {lining}) - {mean_accept}")
    plt.plot(co_tan, percentage_accepted, c="b", label = "Mean acceptance: "+mean_accept)
    
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
            
def pointRepetitionFactor(clustering:str = "", lining:str = "solveS", z0 = 0, n = 16, events=10, savefig=False, ideal=False, show = True): 
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
        cover.solve(clustering=clustering, lining=lining, z0 = z0, n = n, nlines=100, show = False)
        
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

def wedge_test_old(lining:str = "solveS", solve_at = 0, z0 = np.arange(-15, 15.5, 0.5), n = 16, wedges = [0, 128], lines=1000, savefig=False, v = 'v2'):
    mean_list = np.zeros((len(z0), wedges[1]-wedges[0]))
    num_covers = []
    PRF = []
    file = f'wedgeData_{v}_128.txt'

    with open(file) as f:
        for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
            d = np.array(ast.literal_eval(f.readline()))
            env = Environment()
            data = DataSet(env = env, n_points = 150)
            data.input_data(d, add = True)
            cover = Cover(env, data)
            cover.solve(z0 = solve_at, lining=lining, n = n, show = False)
            #data.plot(True)
            num_covers.append(cover.n_patches)
            out = [] 

            for layer in range(env.layers): 
                for point in data.array[layer]: 
                    
                    num_in = 0
                    for patch in cover.patches: 
                        if patch.contains_p(point, layer): 
                            num_in += 1
                            
                    out.append(num_in)
            PRF.append(out)

            for iz, z in enumerate(np.array(z0)):
                percentage_accepted = 0 
                    
                lg = LineGenerator(env, z)
                test_lines = lg.generateEvenGrid(lines)
                
                for i in range(len(test_lines)): 
                    for patch in cover.patches:
                        if patch.contains(test_lines[i]): 
                            percentage_accepted += 1 
                            break 

                
                percentage_accepted = percentage_accepted/lines
                mean_list[iz, ik] = mean_list[iz, ik] + percentage_accepted
            #print(ik)
    mean_accept = format(np.mean(mean_list), ".3f")
    mean_num = format(np.mean(num_covers), ".1f")
    std_num = format(np.std(num_covers), ".1f")
    if type(solve_at) == float:
        ymin = 0
    elif type(solve_at) == int:
        ymin = 0
    else:
        ymin = 0.90

    plt.scatter(z0, np.mean(mean_list, axis = 1), color = 'r', s = 10)
    plt.plot(z0, np.mean(mean_list, axis = 1), color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance',  fontsize = 16)
    plt.ylim(ymin, 1.0)
    plt.title(f'{lining}', fontsize = 16)
    PRFm = format(np.mean(out), '.2f')
    PRFs = format(np.std(out), '.2f')
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\nPatches with " + r'$z_0$' + f" = {np.round(np.array(solve_at), 2)}\nppl = {n}, " + r'$N_{wedges}$ ' + f"= {wedges[1]}, {v} events"],
        loc = 8, fontsize = 12)
    if savefig == True:
        try:
            at = len(solve_at)
        except:
            at = solve_at
        plt.savefig(f"Figures/wedge_test({lining}_{v.replace('.','')}_{at}_n{n})")
    if np.mean(np.mean(mean_list, axis = 1)) > 0.999:
        plt.show()
        return np.mean(np.mean(mean_list, axis = 1))
    else:
        plt.clf()
        return np.mean(np.mean(mean_list, axis = 1))

def wedgeSlopePlot(lining:str = "solveS", events=128, lines=1000, z0 = 0, savefig=False, show = True, v = 'v2'):
    
    percentage_accepted = [0 for _ in range(lines)]

    
    file = open(f'wedgeData_{v}_128.txt')
    '''    for i in range(7):
        line = file.readline()'''
    for k in range(events):
        line = file.readline()
        env = Environment()
        data = DataSet(env, n_points=150)
        d = np.array(ast.literal_eval(line))
        data.input_data(d, add = True)
        cover = Cover(env, data) 
        cover.solve(lining=lining, z0 = z0, show = False)
        
        lg = LineGenerator(env, z0)
        test_lines = lg.generateEvenGrid(lines)
        co_tan = []
        
        
        for i in range(len(test_lines)):
            co_tan.append(100/test_lines[i].slope)
            for patch in cover.patches: 
                if patch.contains(test_lines[i]): 
                    percentage_accepted[i] += 1 
                    break


    percentage_accepted = [x / events for x in percentage_accepted]
    mean_accept = format(np.mean(percentage_accepted), ".3f")
    
    if show == False:
        return np.mean(percentage_accepted)

    print(f"({lining}) - {mean_accept}")
    plt.plot(co_tan, percentage_accepted, c="b", label = "Mean acceptance: "+mean_accept)
    
    plt.title(f"Acceptance Rate ({lining})", fontsize = '20')
    plt.xlabel("dZ/dr", fontsize = '16')
    plt.ylabel("Acceptance Probability", fontsize = '16')
    plt.legend(fontsize = '16')
    if savefig == True: 
        plt.savefig(f"Figures/Acceptance_Rate_({lining})")
    plt.show() 