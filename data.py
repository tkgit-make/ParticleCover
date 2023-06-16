import numpy as np 
import matplotlib.pyplot as plt

class Environment: 
    
    def __init__(self, top_layer_lim = 100, bottom_layer_lim = 15, layers = 5, radii = 5.0): 
        self.top_layer_lim = top_layer_lim
        self.bottom_layer_lim = bottom_layer_lim
        self.layers = layers 
        self.radii = radii   
    

class DataSet(Environment): 
    
    def __init__(self, env:Environment, n_points = 150, equal_spacing = False): 
        # inherit parent attributes 
        self.env = env
        self.n_points = n_points
        limits = np.linspace(env.bottom_layer_lim, env.top_layer_lim, env.layers+1)[1:]
        
        if equal_spacing == False: 
            
            array = np.random.uniform(low=-1.0, high=1.0, size=(n_points, env.layers)) 
            self.array = np.sort((array * limits).T)
            
        else: 
            pnts = np.linspace(-1., 1., n_points)
            col = np.ones(5)
            array, _ = np.meshgrid(pnts, col) 
            self.array = array * limits[:, None] 

    def input_data(self, wedge):
        pass
            
        
    def plot(self, show_lines = False): 
        
        nums = np.arange(1, self.env.layers + 1)
        heights = self.env.layers * np.linspace(nums, nums, self.n_points).T
        
        plt.scatter(self.array, heights, c="g", s=2)
        
        max_height = heights[-1][-1]
        
        plt.yticks(np.arange(0, max_height + 1, self.env.layers))
        
        if show_lines == True: 
            plt.plot([15, 100], [0, max_height], c="r", alpha=0.5)
            plt.plot([-15, 15], [0.0, 0.0], c="r", alpha=0.5)
            plt.plot([-15, -100], [0, max_height], c="r", alpha=0.5)
        
        #plt.show()
       


