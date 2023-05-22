import numpy as np 
import matplotlib.pyplot as plt



def generate_naive_data(ll=-1., ul=1., n_pnts=100, layers=5): 
    
    unf = np.random.rand(layers, n_pnts)
    
    return unf * (ul - ll) + ll


def generate_soph_data(ll=-1., ul=1., n_pnts=100, layers=5): 
    
    first_itcpts = np.random.rand(n_pnts) * (ul - ll) + ll 
    last_itcpts = np.random.rand(n_pnts) * (ul - ll) + ll 

    return np.linspace(first_itcpts, last_itcpts, layers) 


def plot(arr, show_lines=False): 
    num_pnts = 5 * np.ones(arr.shape[1])
    
    for i in range(arr.shape[0]): 
        plt.scatter(arr[i], num_pnts, s=2, c="b")
        num_pnts += 5 * np.ones(arr.shape[1])
    
    plt.yticks(np.arange(0, num_pnts[0] + 1, 5.0))
        
    # plt.show()
    
    if show_lines == True: 
        
        for i in range(arr.shape[1]): 
            plt.plot(arr[:, i], np.linspace(5, 5 * arr.shape[0], 5), c="r")
        
    plt.show() 
    

naive_data = generate_naive_data(n_pnts=100)    
np.savetxt("naive_data.csv", naive_data, delimiter=",")

soph_data = generate_soph_data(n_pnts=100)
np.savetxt("soph_data.csv", naive_data, delimiter=",")