from src.coverers.data_structs import * 
from src.coverers.line import * 
from src.readers.reader import *
from src.coverers.wedgecover import *
from src.debug import * 
import numpy as np 
import matplotlib.pyplot as plt 
import time

def wedge_test(lining:str = "makePatches_Projective_center", apexZ0 = 0, z0_spacing = 0.5, ppl = 16, z0_luminousRegion = 15., wedges = [0, 128], lines=1000, v = 'v3', top_layer_cutoff = 50., accept_cutoff = 10., leftRightAlign=True, uniform_N_points = False, acceptance_method = "Analytic", show_acceptance_of_cover=False, savefig=False, figSizeScale=6):
    """Creates acceptance vs z0 plot
    
    Args:
        lining (str, optional): solving method, default is solveS
        apexZ0 (int, optional): the z values the patches are being made projected to beam axis
        z0 (num or list, optional): array of z0 values we are testing over, default is range of -15 to 15 with 0.5 spacing
        n (int, optional): ppl (point per patch per layer), default is 16
        wedges (list, optional): which wedges, enter in list of [starting wedge, ending wedge]
        lines (int, optional): how many line to test acceptance with at each z0 value
        savefig (bool, optional): True to save figure
        uniform_N_points(False or int): number of points in each layer or False to be not uniform
        v (str, optional): version of data, ensure data file is in directory as "wedgeData_{v}_128.txt"
        acceptance_method : choose between 'Analytic' or 'MonteCarlo'
    """

    accept_cutoff = z0_luminousRegion    
    #create list for z values we're testing the acceptance of, number of covers, and PRF
    showZimperfect = False
    if (wedges[1]-wedges[0]) == 1:
        showZimperfect = True
    if (wedges[1]-wedges[0]) > 50:
        show_acceptance_of_cover = False
        z0_spacing = 0.2
    num_covers = []
    PRF = []
    data_string = f'{v} events'
    if uniform_N_points != False:
        data_string = f'Uniform {uniform_N_points} points'
        wedges = [0, 1]

    zInnerLayer = np.arange(-22, 22+z0_spacing, z0_spacing)
    z0Array = np.arange(-z0_luminousRegion, z0_luminousRegion+z0_spacing, z0_spacing)
    mean_list = np.zeros(( wedges[1]-wedges[0], len(z0Array)))
    z0Imperfect = []
    z0OverEfficiency = []
    #read wedgeData file and create environment
    all_data = readFile(f'python/data/wedgeData_{v}_128.txt', wedges[1])
    #loop through all events
    for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
        print('wedge: ', k)
        #convert to existing data format
        env, points = all_data[k] 
        env = Environment(top_layer_lim = top_layer_cutoff, beam_axis_lim=z0_luminousRegion)
        data = DataSet(env)
        if show_acceptance_of_cover:
            plt.figure(figsize = (1.7*z0_luminousRegion/figSizeScale, top_layer_cutoff/figSizeScale))
        if uniform_N_points == False:
            data.importData(points)
        else:
            data.generateUniform([uniform_N_points, uniform_N_points, uniform_N_points, uniform_N_points, uniform_N_points])
        #add the 1 micron boundary points
        data.addBoundaryPoint()
        #solve for cover
        cover = wedgeCover(env, data)
        cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, leftRight=leftRightAlign, show = False)
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
        sleep_time_between_patches = 0.0
        
        #this loop shows quilt of patches  
        for iz, zIn in enumerate(np.array(zInnerLayer)):
            
            list_of_intersections = []
            for patch in cover.patches: 
                list_of_segs = [pgram.crossSection(zIn) for pgram in patch.parallelograms]
                overlap_of_superpoints = intersection(patch.env, list_of_segs, True) 
                list_of_intersections.append(overlap_of_superpoints)
                
            if show_acceptance_of_cover:
                    
                plt.xlabel(r"$z_1$ (cm)", fontsize = 18)
                plt.ylabel(r"$z_{top}$ (cm)", fontsize = 18)
                plt.title("acceptance of cover", fontsize = 18)
                z1Lim = cover.patches[-1].straightLineProjectorFromLayerIJtoK(-top_layer_cutoff,z0_luminousRegion,env.num_layers,0,1)
                plt.axline((z1Lim, -top_layer_cutoff), (env.trapezoid_edges[0], top_layer_cutoff),linewidth=1, color='black')
                plt.axline((-z1Lim, top_layer_cutoff), (-env.trapezoid_edges[0], -top_layer_cutoff),linewidth=1, color='black')

                colors = ["b", "r", "g", "c", "m", "y", "k", "chocolate", "indigo", "springgreen", "orange", "rosybrown", "tomato","olive", "deeppink"]
                    
                col = 0
                for line in list_of_intersections: 
                    #plt.xlim(-z0_luminousRegion,z0_luminousRegion)
                    plt.xlim(-env.trapezoid_edges[0],env.trapezoid_edges[0])
                    plt.ylim(-top_layer_cutoff, top_layer_cutoff)
                    plt.plot([zIn, zIn], [line.min_z5_accepted, line.max_z5_accepted], c=colors[col % len(colors)], alpha=0.3, linewidth=3)
                    col += 1
                    time.sleep(sleep_time_between_patches)

        #these loops calculate acceptance vs z0
        for iz, z0 in enumerate(np.array(z0Array)):
            
            if acceptance_method == "Analytic": 
                
                list_of_z0intersections = []
                list_of_z0intersectionsCopy = []
                for patch in cover.patches: 
                    # convert to a z0 scan when parameter space is (z1,z5), cast shadow from z0 to layer5 of each superpoint 
                    list_of_segs_z0Scan = [
                        lineSegment(
                            min(patch.env.top_layer_lim,max(-patch.env.top_layer_lim,
                                patch.straightLineProjectorFromLayerIJtoK(z0,SuPoint.min,0,spLayer+1,patch.env.num_layers))),
                            max(-patch.env.top_layer_lim,min(patch.env.top_layer_lim,
                                patch.straightLineProjectorFromLayerIJtoK(z0,SuPoint.max,0,spLayer+1,patch.env.num_layers))),
                        ) for spLayer, SuPoint in enumerate(patch.superpoints)]
                    overlap_of_superpoints_z0Scan = intersection(patch.env, list_of_segs_z0Scan, False) 
                    list_of_z0intersections.append(overlap_of_superpoints_z0Scan)

                list_of_z0intersectionsCopy = list_of_z0intersections
                total_measure = unionOfLineSegments(list_of_z0intersections)
                #print('total_measure:',total_measure)
                percentage_accepted = 100.0*total_measure/(2.0 * patch.env.top_layer_lim)
                if (percentage_accepted < 98.0) and (abs(z0) < z0_luminousRegion):
                    print('wedge: ', k, ' underEfficiency percentage_accepted: ', percentage_accepted, ' z0:', z0)
                    z0Imperfect.append(z0)
                    for seg in list_of_z0intersectionsCopy:
                        print('segment:',seg.min_z5_accepted, seg.max_z5_accepted)
                    #for patch in cover.patches:
                        #if (overlap_of_superpoints_z0Scan.min_z5_accepted < overlap_of_superpoints_z0Scan.max_z5_accepted):
                            #print('overlap_of_superpoints_z0Scan:',overlap_of_superpoints_z0Scan.min_z5_accepted,overlap_of_superpoints_z0Scan.max_z5_accepted)
                if (percentage_accepted > 100.0001) and (abs(z0) < z0_luminousRegion):
                    print('wedge: ', k, ' overEfficiency percentage_accepted: ', percentage_accepted, ' z0:', z0)
                    z0OverEfficiency.append(z0)

                if showZimperfect and show_acceptance_of_cover:
                    zTopVal = [-top_layer_cutoff,top_layer_cutoff]
                    for zImperfect in z0Imperfect:
                        z1Left = cover.patches[-1].straightLineProjectorFromLayerIJtoK(-top_layer_cutoff,zImperfect,env.num_layers,0,1)
                        z1Right = cover.patches[-1].straightLineProjectorFromLayerIJtoK(top_layer_cutoff,zImperfect,env.num_layers,0,1)
                        z1Val = [z1Left,z1Right]
                        plt.plot(z1Val, zTopVal, linewidth=0.01, color='moccasin', alpha=0.4)
                    for zOverEff in z0OverEfficiency:
                        z1Left = cover.patches[-1].straightLineProjectorFromLayerIJtoK(-top_layer_cutoff,zOverEff,env.num_layers,0,1)
                        z1Right = cover.patches[-1].straightLineProjectorFromLayerIJtoK(top_layer_cutoff,zOverEff,env.num_layers,0,1)
                        z1Val = [z1Left,z1Right]
                        plt.plot(z1Val, zTopVal, linewidth=0.01, color='lightskyblue', alpha=0.2)
                    
            elif acceptance_method == "MonteCarlo": 
                percentage_accepted = 0 
                
                lg = LineGenerator(env, z)
                test_lines = lg.generateEvenGrid(lines)
                
                for i in range(len(test_lines)): 
                    for patch in cover.patches:
                        if patch.contains(test_lines[i]):
                            percentage_accepted += 1 
                            break
                
                percentage_accepted = 100.0*percentage_accepted/lines
            
            mean_list[ik, iz] = mean_list[ik, iz] + percentage_accepted
            
        if show_acceptance_of_cover: 
            plt.show()
            plt.close()
        
    mean_num = format(np.mean(num_covers), ".1f")
    std_num = format(np.std(num_covers), ".1f")

    #sets minimum for plot
    if type(apexZ0) == float:
        ymin = 95
    elif type(apexZ0) == int:
        ymin = 99.99
    else:
        ymin = 95

    #creates plots and saves them
    plt.scatter(z0Array, np.mean(mean_list, axis = 0), color = 'r', s = 10)
    plt.plot(z0Array, np.mean(mean_list, axis = 0), color = 'k')
    plt.xlabel(r'$z_0$ [cm]', fontsize = 16)
    plt.ylabel('Acceptance (%)',  fontsize = 16)
    plt.ylim(ymin, 100.0)
    plt.xlim(-z0_luminousRegion,z0_luminousRegion)  
    plt.title(f'{lining}', fontsize = 16)
    mask = np.abs(z0Array) <= accept_cutoff
    PRFm = format(np.mean(out), '.2f')
    PRFs = format(np.std(out), '.2f')
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\n" + r'$apexZ_0$' + f" = {apexZ0}, ppl = {ppl}, " + r"$z_{top}$: "+ f"{top_layer_cutoff}\n" + r'$N_{wedges}$ ' + f"= {wedges[1]}, {data_string}\nAverage non-Acceptance [-{accept_cutoff}, {accept_cutoff}]: {int((100.0-np.mean(mean_list[:, mask]))*10000)} ppm"],
        loc = 8, fontsize = 12)
    if savefig == True:
        try:
            at = len(apexZ0)
        except:
            at = 0
        plt.savefig(f"Figures/wedge_test({lining}_{data_string}_{at}_ppl{ppl}_z0{top_layer_cutoff})")
    plt.show()
    #cover.plot()

def unaccepted_lines(apexZ0:list = [-10, 0, 10], wedge_number = 0, line_origin:list = [-5, 5], accepted = False, unaccepted = True, v = 'v3', top_layer_cutoff = 100., uniform_points = False): 
    filepath = f"data/wedgeData_{v}_128.txt"
    f = open(f'data/{v}_patches.txt')
    filedata = readFile(filepath, stop=128, performance=False)
    env, points = filedata[wedge_number]
    env = Environment(top_layer_lim = top_layer_cutoff)
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
    plt.xlim(-top_layer_cutoff, top_layer_cutoff)
    ds.plot(True, False)

    plt.title(datastring, fontsize = 20)
    plt.show()

def minimal_cover_binary_search(lining:str = "makePatches_Projective_center", accept = 0.999, start = 'odd', ppl = 16, wedges = 128, z_top = 50., z0_spacing = 0.5, z0_luminousRegion = 15., v = 'v3', savefig = False):
    if start == 'odd':
        apexZ0 = [0]
        real_solve = [0]
    else:
        apexZ0 = [0.5, 0.5]
        real_solve = []     
    reached = False
    plt.figure(figsize=(10,7))
    z0 = np.arange(-15, 15+z0_spacing, z0_spacing)
    lines = 1000
    left_floor = 0
    left_ceiling = int(len(z0)/2)
    right_floor = int(len(z0)/2 + 1)
    right_ceiling = int(len(z0))
    current_tries_left = []
    current_tries_right = []
    left_middle = left_floor
    right_middle = right_ceiling
    data_string = f"{v} events"
    left_stop = False
    right_stop = False
    while reached == False:
        mean_list = np.zeros((len(z0), wedges))
        num_covers = []
        PRF = []
        file = f'python/data/wedgeData_{v}_128.txt'
        all_data = readFile(file, wedges)

        for k in range(wedges):
            env, points = all_data[k] 
            env = Environment(top_layer_lim=z_top, beam_axis_lim=z0_luminousRegion)
            data = DataSet(env)
            data.importData(points)
            #add the 1 micron points
            data.addBoundaryPoint()
            #solve for cover
            cover = wedgeCover(env, data)
            cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, show = False)
            num_covers.append(cover.n_patches)
            out = [] 
            x_edges = np.array(env.radii)*(env.top_layer_lim-env.beam_axis_lim)/(env.radii[-1]) + env.beam_axis_lim
            for layer in range(env.num_layers): 
                for point in data.array[layer]: 
                    
                    num_in = 0
                    for patch in cover.patches: 
                        if patch.contains_p(point, layer): 
                            num_in += 1
                    if (np.abs(point.z) <= x_edges[layer]+0.1) & (num_in != 0):                       
                        out.append(num_in)
            PRF.append(out)

            for iz, z in enumerate(np.array(z0)):
                '''
                percentage_accepted = 0 
                    
                lg = LineGenerator(env, z)
                test_lines = lg.generateEvenGrid(lines)
                
                for i in range(len(test_lines)): 
                    for patch in cover.patches:
                        if patch.contains(test_lines[i]): 
                            percentage_accepted += 1 
                            break 

                
                percentage_accepted = percentage_accepted/lines
                '''
                list_of_intersections = []
                for patch in cover.patches: 
                    list_of_segs = [pgram.crossSection(z) for pgram in patch.parallelograms]
                    overlap_of_superpoints = intersection(patch.env, list_of_segs) 
                    list_of_intersections.append(overlap_of_superpoints)
                
                total_measure = unionOfLineSegments(list_of_intersections)
                
                percentage_accepted = total_measure/(2.0 * patch.env.top_layer_lim)
                mean_list[iz, k] = mean_list[iz, k] + percentage_accepted

        z0_means = np.mean(mean_list, axis = 1)
        last_apexZ0 = apexZ0
        if np.all(z0_means[:int(len(z0)/2+1)] >= accept):
            left_stop = True
        if np.all(z0_means[int(len(z0)/2+1):]>= accept):
            right_stop = True
        if (left_stop == True) & (right_stop == True):
            reached = True
        if np.all(z0_means[left_middle:int(len(z0)/2)+1] >= accept):
            if left_stop == False:
                current_tries_left.append(left_middle)
            left_ceiling = left_middle
            debug("Michelle", f"yes: {z0[left_middle]}")
        else:
            left_floor = left_middle
        if np.all(z0_means[int(len(z0)/2)+1:right_middle+1] >= accept):
            if right_stop == False:
                current_tries_right.append(right_middle)
            right_floor = right_middle
            debug("Michelle", f"yes: {z0[right_middle]}")
        else:
            right_ceiling = right_middle
        
        if left_ceiling - left_floor <= 1:
            if left_stop == False:
                real_solve = np.append(real_solve, z0[min(current_tries_left)])
                debug("Michelle", f"New left: {z0[min(current_tries_left)]}")
                current_tries_left = []
                left_floor = -1
                left_ceiling = left_middle
        if right_ceiling - right_floor <= 1:
            if right_stop == False:
                real_solve = np.append(real_solve, z0[max(current_tries_right)])
                debug("Michelle", f"New right: {z0[max(current_tries_right)]}")
                current_tries_right = []
                right_floor = right_middle
                right_ceiling = int(len(z0))

        real_solve = np.unique(real_solve)
        left_middle = int(np.round((left_floor+left_ceiling)/2))
        right_middle = int(np.round((right_floor+right_ceiling)/2))
        apexZ0 = np.copy(real_solve)
        apexZ0 = np.append(apexZ0, z0[left_middle])
        apexZ0 = np.append(apexZ0, z0[right_middle])
        apexZ0 = np.unique(apexZ0)
        
        debug("Michelle", f"real: {real_solve}")
        debug("Michelle", f"next try: {apexZ0}")
        debug("Muchang", f"{real_solve}")
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
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPRF: {PRFm}" + r'$\pm$' + f"{PRFs}, " + r"$z_{top}$ = " +f"{z_top}\n" r'$apexZ_0$' + f" = {last_apexZ0}\n" + r'$N_{wedges}$ ' + f"= {wedges}, {data_string}, ppl = {ppl}\nAverage non-acceptance: {np.round((1.0-np.mean(mean_list))*100, 3)}%"],
        loc = 8, fontsize = 12)
    if savefig == True:
        plt.savefig(f"python/Figures/binary_search_{start}_({lining}_{data_string}_ztop_{z_top}_ppl{16}_wedges{wedges}_z_lum{z0_luminousRegion})")
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

def minimal_cover_linear_search(lining:str = "makePatches_Projective_center", accept = 0.999, start = 'odd', ppl = 16, wedges = [0, 128], v = 'v3', z_5=100., savefig = False):

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
        file = f'python/data/wedgeData_{v}_128.txt'
        all_data = readFile(file, wedges[1])

        for ik, k in enumerate(np.arange(wedges[0], wedges[1])):
            env, points = all_data[k] 
            env = Environment(top_layer_lim=z_5)
            data = DataSet(env)
            data.importData(points)
            #add the 1 micron points
            data.addBoundaryPoint()
            #solve for cover
            cover = wedgeCover(env, data)
            cover.solve(apexZ0 = apexZ0, lining=lining, ppl = ppl, show = False)
            num_covers.append(cover.n_patches)
            out = [] 
            x_edges = np.array(env.radii)*(env.top_layer_lim-env.beam_axis_lim)/(env.radii[-1]) + env.beam_axis_lim
            for layer in range(env.num_layers): 
                for point in data.array[layer]: 
                    
                    num_in = 0
                    for patch in cover.patches: 
                        if patch.contains_p(point, layer): 
                            num_in += 1
                    if (np.abs(point.z) <= x_edges[layer]+0.1) & (num_in != 0):
                        out.append(num_in)
            PRF.append(out)

            for iz, z in enumerate(np.array(z0)):
                '''
                percentage_accepted = 0 
                    
                lg = LineGenerator(env, z)
                test_lines = lg.generateEvenGrid(lines)
                
                for i in range(len(test_lines)): 
                    for patch in cover.patches:
                        if patch.contains(test_lines[i]): 
                            percentage_accepted += 1 
                            break
                percentage_accepted = percentage_accepted/lines
                '''
                list_of_intersections = []
                for patch in cover.patches: 
                    list_of_segs = [pgram.crossSection(z) for pgram in patch.parallelograms]
                    overlap_of_superpoints = intersection(patch.env, list_of_segs) 
                    list_of_intersections.append(overlap_of_superpoints)
                
                total_measure = unionOfLineSegments(list_of_intersections)
                
                percentage_accepted = total_measure/(2.0 * patch.env.top_layer_lim)

                
                
                mean_list[iz, ik] = mean_list[iz, ik] + percentage_accepted
        last_apexZ0 = apexZ0
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
    plt.legend([f"Number of Patches: {mean_num}" + r'$\pm$' + f"{std_num}\nPoint Repetition Factor: {PRFm}" + r'$\pm$' + f"{PRFs}\n" + r'$apexZ_0$' + f" = {last_apexZ0}, ppl = {ppl}\n" + r'$N_{wedges}$ ' + f"= {wedges[1]}, {data_string}\nAverage Acceptance: {np.round(np.mean(mean_list)*100, 2)}%"],
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
        #add the 1 micron points
        data.addBoundaryPoint()
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
        #add the 1 micron points
        data.addBoundaryPoint()
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
        #add the 1 micron points
        data.addBoundaryPoint()
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
        #add the 1 micron points
        data.addBoundaryPoint()
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

