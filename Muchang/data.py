import numpy as np 
import matplotlib.pyplot as plt

class Environment: 
    
    def __init__(self, top_layer_lim = 1.0, bottom_layer_lim = 0.15, layers = 5, radii = 5.0): 
        self.top_layer_lim = top_layer_lim
        self.bottom_layer_lim = bottom_layer_lim
        self.layers = layers 
        self.radii = radii   
    


class DataSet(Environment): 
    
    def __init__(self, env:Environment, n_points = 150): 
        # inherit parent attributes 
        
        self.env = env
        
        self.n_points = n_points
        limits = np.linspace(env.bottom_layer_lim, env.top_layer_lim, env.layers+1)[1:]
        
        array = np.random.uniform(low=-1.0, high=1.0, size=(n_points, env.layers)) 
        
        self.array = np.sort((array * limits).T)
        
        
    def plot(self, show_lines = False): 
        
        nums = np.arange(1, self.env.layers + 1)
        heights = self.env.layers * np.linspace(nums, nums, self.n_points).T
        
        plt.scatter(self.array, heights, c="b", s=2)
        
        max_height = heights[-1][-1]
        
        plt.yticks(np.arange(0, max_height + 1, self.env.layers))
        
        if show_lines == True: 
            plt.plot([0.15, 1.0], [0, max_height], c="r", alpha=0.5)
            plt.plot([-0.15, 0.15], [0.0, 0.0], c="r", alpha=0.5)
            plt.plot([-0.15, -1.0], [0, max_height], c="r", alpha=0.5)
        
        # plt.show()
       


