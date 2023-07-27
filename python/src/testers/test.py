from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.readers.reader import *
from src.coverers.wedgecover import *
import numpy as np 
import matplotlib.pyplot as plt 

def wedge_test(lining:str = "makePatches_Projective", apexZ0 = 0, z0 = np.arange(-15, 15.5, 0.5), ppl = 16, z0_luminousRegion = 15., wedges = [0, 128], lines=1000, v = 'v3', z0_cutoff = 100., accept_cutoff = 10., uniform_N_points = False, savefig=False):
    """Creates acceptance vs z0 plot
    
    Args:
        lining (str, optional): solving method, default is solveS
        apexZ0 (int, optional): the z values the patches are being made in accordance to
        z0 (num or list, optional): array of z0 values we are testing over, default is range of -15 to 15 with 0.5 spacing
        n (int, optional): ppl (point per patch per layer), default is 16
        wedges (list, optional): which wedges, enter in list of [starting wedge, ending wedge]
        lines (int, optional): how many line to test acceptance with at each z0 value
        savefig (bool, optional): True to save figure
        uniform_N_points(False or int): number of points in each layer or False to be not uniform
        v (str, optional): version of data, ensure data file is in directory as "wedgeData_{v}_128.txt"
    """

    #create list for z values we're testing the acceptance of, number of covers, and PRF
    mean_list = np.zeros(( wedges[1]-wedges[0], len(z0)))
    num_covers = []
    PRF = []
    data_string = f'{v} events'
    if uniform_N_points != False:
        data_string = f'Uniform {uniform_N_points} points'
        wedges = [0, 1]

    #read wedgeData file and create environment
    all_data = readFile(f'data/wedgeData_{v}_128.txt', wedges[1])
    #loop through all events
    for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
        #convert to existing data format
        env, points = all_data[k] 
        env = Environment(top_layer_lim = z0_cutoff, beam_axis_lim=z0_luminousRegion)
        data = DataSet(env)
        if uniform_N_points == False:
            data.importData(points)
        else:
            data.generateUniform([uniform_N_points, uniform_N_points, uniform_N_points, uniform_N_points, uniform_N_points])
        #add the 0.1 cm points
        data.add()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, show = False)
        #append number of covers in the patch
        num_covers.append(cover.n_patches)
        out = [] 

        #these loops calculate PRF
        for layer in range(env.num_layers): 
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
            mean_list[ik, iz] = mean_list[ik, iz] + percentage_accepted
    
    mean_num = format(np.mean(num_covers), ".1f")
    std_num = format(np.std(num_covers), ".1f")

    #sets minimum for plot
    if type(apexZ0) == float:
        ymin = 0
    elif type(apexZ0) == int:
        ymin = 0
    else:
        ymin = 0.9

    #creates plots and saves them
    plt.scatter(z0, np.mean(mean_list, axis = 0), color = 'r', s = 10)
    plt.plot(z0, np.mean(mean_list, axis = 0), color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance',  fontsize = 16)
    plt.ylim(ymin, 1.0)
    plt.title(f'{lining}', fontsize = 16)
    mask = np.abs(z0) <= accept_cutoff
    PRFm = format(np.mean(out), '.2f')
    PRFs = format(np.std(out), '.2f')
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\n" + r'$apexZ_0$' + f" = {apexZ0}, ppl = {ppl}\n" + r'$N_{wedges}$ ' + f"= {wedges[1]}, {data_string}\nAverage Acceptance [-{accept_cutoff}, {accept_cutoff}]: {np.round(np.mean(mean_list[:, mask])*100, 2)}%"],
        loc = 8, fontsize = 12)
    if savefig == True:
        try:
            at = len(apexZ0)
        except:
            at = 0
        plt.savefig(f"Figures/wedge_test({lining}_{data_string}_{at}_ppl{ppl})")
    plt.show()

def unaccepted_lines(apexZ0:list = [-10, 0, 10], wedge_number = 0, line_origin:list = [-5, 5], accepted = False, unaccepted = True, v = 'v3', z0_cutoff = 100., uniform_points = False):
    filepath = f"data/wedgeData_{v}_128.txt"
    f = open(f'data/{v}_patches.txt')
    filedata = readFile(filepath, stop=128, performance=False)
    env, points = filedata[wedge_number]
    env = Environment(top_layer_lim = z0_cutoff)
    ds = DataSet(env)
    datastring = f"Wedge {wedge_number} Event {v}"

    if uniform_points != False:
        datastring = f'Uniform {uniform_points} Points'
        ds.generateUniform([uniform_points,uniform_points,uniform_points, uniform_points,uniform_points])
    else:
        ds.importData(points)
    ds.add()
    colors = ['b', 'orange', 'm', 'c', 'k']
    plt.figure(figsize = (10, 8))
    lines_to_plot = []
    fitting_lines = []
    for z in line_origin:
        lGen = LineGenerator(env, z)
        fitting_lines = fitting_lines + lGen.generateEvenGrid(100)
    for i in range(len(apexZ0)):
        cov = wedgeCover(env, ds)
        cov.makePatches_Projective_center(apexZ0 = apexZ0[i])
        for j, patch in enumerate(cov.patches):
            top = env.radii[-1]
            lambdaZ_list_left = []
            lambdaZ_list_right = []
            to_plot_l = []
            to_plot_r = []
            for layer in range(cov.env.num_layers):
                to_plot_l.append(min(patch.superpoints[layer].z_values))
                to_plot_r.append(max(patch.superpoints[layer].z_values))
                left = min(patch.superpoints[layer].z_values)
                right = max(patch.superpoints[layer].z_values)   
                lambdaZ_list_left.append((left - apexZ0[i])/env.radii[layer])
                lambdaZ_list_right.append((right - apexZ0[i])/env.radii[layer])
            min_lambdaZ = np.max(lambdaZ_list_left)
            max_lambdaZ = np.min(lambdaZ_list_right)
            to_plot_r.reverse()
            
            #patch.plot(color = colors[int(j%5)])
            if j ==0:
                label = r"$apexZ_0$ = " + f"{apexZ0[i]}"
            else:
                label = '_'
            flipped_radii = [25, 20, 15, 10, 5]
            new_radiis = env.radii+flipped_radii 
            plt.plot(to_plot_l + to_plot_r, new_radiis, color = colors[i], label = label, alpha = 0.5)
            plt.fill(to_plot_l + to_plot_r, new_radiis, color = colors[i], alpha = 0.2)
            #plt.plot([apexZ0[i],min_lambdaZ*top+apexZ0[i], max_lambdaZ*top+apexZ0[i], apexZ0[i]], [0, top, top, 0], color = colors[i], label = label, alpha = 0.5)
            #plt.fill_between([apexZ0[i],min_lambdaZ*top+apexZ0[i], max_lambdaZ*top+apexZ0[i], apexZ0[i]], [0, top, top, 0], color = colors[i], alpha = 0.2)
            for line_index, line in enumerate(fitting_lines):
                if patch.contains(line) == True:
                    lines_to_plot.append(line_index)
    lines_to_plot = np.unique(lines_to_plot)
    for l, line in enumerate(fitting_lines):
        if l in lines_to_plot:
            if accepted == True:
                line.plot('g')
            else:
                pass
        else:
            if unaccepted == True:
                line.plot('r')
            else:
                pass
    if accepted == True:
        plt.plot([0],color = 'g', label = 'Lines Accepted')
    if unaccepted == True:
        plt.plot([0],color = 'r', label = f'Lines Not Accepted from {line_origin}')
    
    plt.legend(loc = 'lower left', fontsize = 12)
    plt.xlim(-z0_cutoff, z0_cutoff)
    ds.plot(True, False)

    plt.title(datastring, fontsize = 20)
    plt.show()

def minimal_cover_binary_search(lining:str = "makePatches_Projective", accept = 0.999, start = 'odd', ppl = 16, wedges = 128, v = 'v3', savefig = False):
    if start == 'odd':
        apexZ0 = [0]
        real_solve = [0]
    else:
        apexZ0 = []
        real_solve = []     
    reached = False
    z0 = np.arange(-15, 15.5, 0.5)
    lines = 1000
    left_index = 0
    right_index = int(len(z0)-1)
    last_left = int(len(z0)/2)
    last_right = int(len(z0)/2 + 1)
    data_string = f"{v} events"

    while reached == False:
        mean_list = np.zeros((len(z0), wedges))
        num_covers = []
        PRF = []
        file = f'data/wedgeData_{v}_128.txt'
        all_data = readFile(file, wedges)

        for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
            env, points = all_data[k] 
            data = DataSet(env)
            data.importData(points)
            #add the 0.1 cm points
            data.add()
            #solve for cover
            cover = wedgeCover(env, data)
            cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, show = False)
            num_covers.append(cover.n_patches)
            out = [] 

            for layer in range(env.num_layers): 
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
        z0_means = np.mean(mean_list, axis = 1)
        if np.all(z0_means > accept):
            break
        if np.all(z0_means[left_index:last_left] > accept):
            real_solve = np.append(real_solve, [z0[left_index]])
            print('new left', z0[left_index])
            last_left = left_index
            left_index = -1
        if np.all(z0_means[last_right:right_index] > accept):
            real_solve = np.append(real_solve, [z0[right_index]])
            print('new right', z0[right_index])
            last_right = right_index
            right_index = len(z0)
        left_index +=1
        right_index -=1
        real_solve = np.unique(real_solve)
        apexZ0 = np.copy(real_solve)

        apexZ0 = np.append(apexZ0, z0[left_index])
        apexZ0 = np.append(apexZ0,z0[right_index])
        apexZ0 = np.unique(apexZ0)
        print(left_index)
        print(right_index)
        print(np.sort(real_solve))

    ymin = accept
    mean_num = format(np.mean(num_covers), ".1f")
    std_num = format(np.std(num_covers), ".1f")

    plt.scatter(z0, np.mean(mean_list, axis = 1), color = 'r', s = 10)
    plt.plot(z0, np.mean(mean_list, axis = 1), color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance',  fontsize = 16)
    plt.ylim(ymin, 1.0)
    plt.title(f'{lining}', fontsize = 16)
    PRFm = format(np.mean(out), '.2f')
    PRFs = format(np.std(out), '.2f')
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\n" + r'$apexZ_0$' + f" = {apexZ0}, ppl = {ppl}\n" + r'$N_{wedges}$ ' + f"= {wedges[1]}, {data_string}\nAverage Acceptance: {np.round(np.mean(mean_list)*100, 2)}%"],
        loc = 8, fontsize = 12)
    if savefig == True:
        plt.savefig(f"Figures/min_cover_binary_{start}_({lining}_{data_string}_ppl{16})")
    plt.show()

def patch_ending_layer(lining = 'makePatches_Projective_center', apexZ0 = [-10, 0, 10], wedges = 1280, z_5 = 50., v = 'v3'):
    filepath = f"data/wedgeData_{v}_128.txt"
    filedata = readFile(filepath, stop=wedges, performance=False)
    ends = []
    plt.figure(figsize=(10, 7))
    for i in range(wedges):
        env, points = filedata[i]
        env = Environment(top_layer_lim = z_5)
        ds = DataSet(env)
        ds.importData(points)
        ds.add()
        cov = wedgeCover(env, ds)
        cov.solve('makePatches_Projective_center',apexZ0 = apexZ0)
        for patch in cov.patches:
            ends.append(patch.end_layer)
    plt.hist(ends, bins = np.arange(0.5, 6.5, 1), rwidth=0.8,edgecolor='black',
             label = f"Method: {lining}\n" + r"z_5: " + f"{[-z_5, z_5]}\n" + r"$apexZ_0$: " + f"{apexZ0} \nData: {v} Events with {wedges} wedges"
             )
    plt.title('Patch Ending Layer', fontsize = 20)
    plt.ylabel('Count', fontsize = 15)
    plt.xlabel('Layer Number', fontsize = 15)
    plt.legend(loc = 'upper left', fontsize = 12)
    plt.show()

def z099(lining:str = "makePatches_Projective", accept = 0.999, start = 'odd', ppl = 16, wedges = [0, 128], v = 'v3', savefig = False):
    #PUT IN BINARY SEARCH INSTEAD
    #Binary search should take 4 steps, whole thing should take 10 steps
    if start == 'odd':
        apexZ0 = [0]
        real_solve = [0]
    else:
        apexZ0 = []
        real_solve = []     
    reached = False
    z0 = np.arange(-15, 15.5, 0.5)
    lines = 1000
    left_index = 0
    right_index = int(len(z0)-1)
    last_left = int(len(z0)/2)
    last_right = int(len(z0)/2 + 1)
    data_string = f"{v} events"
    rf = False
    lf = False

    while reached == False:
        mean_list = np.zeros((len(z0), wedges[1]-wedges[0]))
        num_covers = []
        PRF = []
        file = f'data/wedgeData_{v}_128.txt'
        all_data = readFile(file, wedges[1])

        for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
            env, points = all_data[k] 
            data = DataSet(env)
            data.importData(points)
            #add the 0.1 cm points
            data.add()
            #solve for cover
            cover = wedgeCover(env, data)
            cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, show = False)
            num_covers.append(cover.n_patches)
            out = [] 

            for layer in range(env.num_layers): 
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
        z0_means = np.mean(mean_list, axis = 1)
        if np.all(z0_means > accept):
            break
        if np.all(z0_means[left_index:last_left] > accept):
            real_solve = np.append(real_solve, [z0[left_index]])
            print('new left', z0[left_index])
            last_left = left_index
            left_index = -1
        if np.all(z0_means[last_right:right_index] > accept):
            real_solve = np.append(real_solve, [z0[right_index]])
            print('new right', z0[right_index])
            last_right = right_index
            right_index = len(z0)
        left_index +=1
        right_index -=1
        real_solve = np.unique(real_solve)
        apexZ0 = np.copy(real_solve)

        apexZ0 = np.append(apexZ0, z0[left_index])
        apexZ0 = np.append(apexZ0,z0[right_index])
        apexZ0 = np.unique(apexZ0)
        print(left_index)
        print(right_index)
        print(np.sort(real_solve))

    ymin = accept
    mean_num = format(np.mean(num_covers), ".1f")
    std_num = format(np.std(num_covers), ".1f")

    plt.scatter(z0, np.mean(mean_list, axis = 1), color = 'r', s = 10)
    plt.plot(z0, np.mean(mean_list, axis = 1), color = 'k')
    plt.xlabel('z0 offset [cm]', fontsize = 16)
    plt.ylabel('Acceptance',  fontsize = 16)
    plt.ylim(ymin, 1.0)
    plt.title(f'{lining}', fontsize = 16)
    PRFm = format(np.mean(out), '.2f')
    PRFs = format(np.std(out), '.2f')
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\n" + r'$apexZ_0$' + f" = {apexZ0}, ppl = {ppl}\n" + r'$N_{wedges}$ ' + f"= {wedges[1]}, {data_string}\nAverage Acceptance: {np.round(np.mean(mean_list)*100, 2)}%"],
        loc = 8, fontsize = 12)
    if savefig == True:
        try:
            at = len(apexZ0)
        except:
            at = apexZ0
        plt.savefig(f"Figures/minimal_z0_{start}_({lining}_{data_string}_ppl{16})")
    plt.show()

def numCovers(lining:str = "makePatches_Projective", events=1000, apexZ0 = 0, ppl = 16, savefig=False, v = 'v3'): 
    """Returns histogram of the number of patches

    Args:
        lining (str, optional): Solving method. Defaults to "makePatches_Projective".
        events (int, optional): number of wedges. Defaults to 1000.
        apexZ0 (num or list, optional): z0 locations patches are generated at. Defaults to 0.
        n (int, optional): points per patch per layer. Defaults to 16.
        savefig (bool, optional): Saves figure in Figures folder if True. Defaults to False.
        v (str, optional): Data version type. Defaults to 'v3'.
    """
    num_covers = [] 
    file = f'data/wedgeData_{v}_128.txt'
    all_data = readFile(file, events)
    for k in range(events): 
        env, points = all_data[k] 
        data = DataSet(env)
        data.importData(points)
        #add the 0.1 cm points
        data.add()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, show = False)
        
        num_covers.append(cover.n_patches)
    avg = np.mean(num_covers)
    std = np.std(num_covers)
    plt.hist(num_covers, 
                bins=np.arange(np.min(num_covers), np.max(num_covers)+2)-0.5,
                edgecolor='black', 
                rwidth=0.8, label = f"mean: {format(avg, '.2f')}, stdev: {format(std, '.2f')}"
            )
    print(f"({lining}) - {format(avg, '.2f')}, {format(std, '.2f')}")
    plt.title(f"Number of Patches per Cover ({lining})", fontsize = '20')
    plt.xlabel("Number of Patches", fontsize = '16')
    plt.ylabel("Number of Covers", fontsize = '16')
    plt.legend(fontsize = '20')
    if savefig == True: 
        plt.savefig(f"Figures/nPatches_({lining})")
    plt.show()    
    
def acceptSlopePlot(lining:str = "makePatches_Projective", events=128, lines=1000, z0= 0, ppl = 16, apexZ0=0, v = 'v3', savefig=False, show = True):
    """Generates plot of acceptance for a given z0 value

    Args:
        lining (str, optional): Solving method. Defaults to "makePatches_Projective".
        events (int, optional): number of wedges. Defaults to 128.
        lines (int, optional): Number of lines tested. Defaults to 1000.
        z0 (num, optional): Where lines are coming from on z axis. Defaults to 0.
        n (int, optional): points per patch per layer. Defaults to 16.
        apexZ0 (num or list, optional): z0 locations patches are generated at. Defaults to 0.
        v (str, optional): data version. Defaults to 'v3'.
        savefig (bool, optional): Saves figure in Figures folder if True. Defaults to False.
        show (bool, optional): Set to False to just return acceptance list without plot. Defaults to True.

    Returns:
        list: Acceptance list
    """
    
    percentage_accepted = [0 for _ in range(lines)] 
    file = f'data/wedgeData_{v}_128.txt'
    all_data = readFile(file, events)
    
    for k in range(events): 
        env, points = all_data[k] 
        data = DataSet(env)
        data.importData(points)
        #add the 0.1 cm points
        data.add()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(lining=lining, apexZ0 = apexZ0, ppl = ppl, show = False)
        
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
                pass

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
            
def pointRepetitionFactor(lining:str = "makePatches_Projective", events=128, apexZ0 = 0, ppl = 16, savefig=False, v = 'v3', show = True): 
    """for every event, we loop through all the points in the dataset and 
        compute how many patches contain that point. The lower in general 
        the better, since this is a metric of non-wastefulness 

    Args:
        lining (str, optional): Solving method. Defaults to "makePatches_Projective".
        events (int, optional): number of wedges. Defaults to 128.
        apexZ0 (num or list, optional): z0 locations patches are generated at. Defaults to 0.
        n (int, optional): points per patch per layer. Defaults to 16.
        savefig (bool, optional): Saves figure in Figures folder if True. Defaults to False.
        v (str, optional): data version. Defaults to 'v3'.
        show (bool, optional): Set to False to just return acceptance list without plot. Defaults to True.

    Returns:
        tuple: mean PRF and standard deviation
    """

    out = []
    file = f'data/wedgeData_{v}_128.txt'
    all_data = readFile(file, events)

    for k in range(events): 
        env, points = all_data[k] 
        data = DataSet(env)
        data.importData(points)
        #add the 0.1 cm points
        data.add()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(lining=lining, apexZ0 = apexZ0, ppl = ppl, show = False)

        
        out2 = [] 
        
        for layer in range(env.num_layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                for patch in cover.patches: 
                    if patch.contains_p(point, layer): 
                        num_in += 1
                        
                out2.append(num_in) 
                
        out += out2

    if show == False:
        return (format(np.mean(out), '.2f'), format(np.std(out), '.2f'))

    print(f"({lining}) mean - {format(np.mean(out), '.2f')}")
    print(f"({lining}) stdev - {format(np.std(out), '.2f')}")


        
    plt.hist(out, bins=np.arange(11) - 0.5, 
             edgecolor='black', 
             label = f"mean: {format(np.mean(out), '.2f')}, stdev: {format(np.std(out), '.2f')}",
             rwidth=0.8
            ) 
    plt.xlabel("Number of Covering Patches", fontsize = '16')
    plt.ylabel("Number of Points", fontsize = '16')
    plt.title(f"Point Repetition Factor ({lining})", fontsize = '20')
    plt.legend(fontsize = '16')
    if savefig == True:
        plt.savefig(f"Figures/Point_Repetition_Factor_({lining})")
    plt.show()
    return (format(np.mean(out), '.2f'), format(np.std(out), '.2f'))

def pointRepetitionFactorLayer(lining:str = "makePatches_Projective_center", wedges=128, apexZ0 = 0, z_5 = 100., ppl = 16, savefig=False, v = 'v3'): 
    """for every event, we loop through all the points in the dataset and 
        compute how many patches contain that point. The lower in general 
        the better, since this is a metric of non-wastefulness 

    Args:
        lining (str, optional): Solving method. Defaults to "makePatches_Projective".
        events (int, optional): number of wedges. Defaults to 128.
        apexZ0 (num or list, optional): z0 locations patches are generated at. Defaults to 0.
        n (int, optional): points per patch per layer. Defaults to 16.
        savefig (bool, optional): Saves figure in Figures folder if True. Defaults to False.
        v (str, optional): data version. Defaults to 'v3'.
        show (bool, optional): Set to False to just return acceptance list without plot. Defaults to True.

    Returns:
        tuple: mean PRF and standard deviation
    """

    file = f'data/wedgeData_{v}_128.txt'
    all_data = readFile(file, wedges)
    env = Environment(z_5)
    plt.figure(figsize=(27, 7))
    x_edges = np.array(env.radii)*(env.top_layer_lim-env.beam_axis_lim)/(env.radii[-1]) + env.beam_axis_lim
    out = [[],[],[],[],[]]
    for k in range(wedges): 
        env, points = all_data[k] 
        env = Environment(z_5)
        data = DataSet(env)
        data.importData(points)
        #add the 0.1 cm points
        data.add()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(lining=lining, apexZ0 = apexZ0, ppl = ppl, show = False)

        
        out2 = [[],[],[],[],[]]
        
        for layer in range(env.num_layers): 
            for point in data.array[layer]: 
                
                num_in = 0
                for patch in cover.patches:
                    if patch.contains_p(point, layer): 
                        num_in += 1        
                if (np.abs(point.z) <= x_edges[layer]) & (num_in != 0):
                    out2[layer].append(num_in)
            out[layer] = out[layer] + out2[layer]

    ylim = max([len(out[x]) for x in range(5)])
    for layer in range(env.num_layers):
        plt.subplot(1, 5, layer+1)
        plt.hist(out[layer],
                edgecolor='black', bins=np.arange(6) + 0.5, 
                label = f"mean: {format(np.mean(out[layer]), '.2f')}\nstdev: {format(np.std(out[layer]), '.2f')}", rwidth=0.8)
        plt.xlabel("Number of Covering Patches", fontsize = '16')
        plt.ylabel("Number of Points", fontsize = '16')
        plt.title(f"Layer {layer+1}", fontsize = '18')
        plt.legend(fontsize = '16')
        plt.ylim((0, ylim*2/3))
    plt.suptitle(f"PRF {lining}", fontsize = '20')
    if savefig == True:
        plt.savefig(f"Figures/Point_Repetition_Factor_layer_({lining})")
    plt.show()

