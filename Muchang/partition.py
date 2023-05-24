import numpy as np 
import matplotlib.pyplot as plt 
from generate_data import * 
import math 

data = generate_data(n_pnts=96)    

def partition(array, ll, ul): 
    
    n_points = array.shape[1] 
    layers = array.shape[0]
    
    array = np.sort(array) 
    
    partitions = [[] for _ in range(layers)] 
    
    for layer in range(layers): 
        positives = np.sort([x for x in array[layer] if x > 0])
        negatives = np.sort([x for x in array[layer] if x < 0])[::-1] 
        
        
        middle = []
        for i in range(8): 
            middle.append((positives[0], layers - layer)) 
            positives = np.delete(positives, 0)
            middle.append((negatives[0], layers - layer))
            negatives = np.delete(negatives, 0) 
        partitions[layer].append(middle) 
        
        
        # clear out the positives 
        
        while len(positives) >= 16: 
            positive = [] 
            for i in range(16): 
                positive.append((positives[0], layers - layer)) 
                positives = np.delete(positives, 0)
            partitions[layer].append(positive)

                    
        while len(negatives) >= 16: 
            negative = [] 
            for i in range(16): 
                negative.append((negatives[0], layers - layer)) 
                negatives = np.delete(negatives, 0) 
            partitions[layer].append(negative)
                
    for layer in partitions: 
        print(len(layer))
    
    for p_num in range(len(partitions[0])): 
        scatt = []
        for lay in partitions: 
            scatt += lay[p_num]
        
        plt.scatter(*zip(*scatt), s=5) 
        
    plt.show() 
        
    
    
        
        
                    
                
        
            
        
        
        
    return 
    
    for layer in range(layers): 
        # sort by their absolute values
        array[layer] = sorted(array[layer], key=lambda x : np.abs(x)) 
        
    
    partitions = []
    
    middle_partition = []
    
    for layer in range(layers): 
        
        sorted(array[layer], key=lambda x : np.abs(x)) 
        # take the first 16 out 
        i = 0 
        while i < 16: 
            middle_partition.append((array[layer][i], layers - layer))
            array[layer].remove(array[layer][i])
            i += 1
            
    positive_partitions = [] 
    
    for layer in range(layers): 
        # take the first 16 out 
        i = 0 
        while i < 16 or len(array[layer]) != 0: 
            if array[layer][i] > 0: 
                positive_partitions.append((array[layer][i], layers - layer))
                array[layer].remove(array[layer][i]) 
                i += 1
        
    
            
    print(middle_partition)
            
    
    # partions = []
    
    # small_slices = [] 
    
    
    
    # k = 0
    # while k + 16 <= array.shape[1]: 
    #     relevant = []
    #     for layer in range(layers): 
    #         relevant += [(array[layer][i], layers - layer) for i in range(k, min(k + 16, len(array[layer]) + 1))]
            
    #     k += 16
    #     small_slices.append(relevant)

    
    # for set in small_slices: 
    #     plt.scatter(*zip(*set), s=3) 
        
    # plt.show() 


partition(data, -1, 1) 