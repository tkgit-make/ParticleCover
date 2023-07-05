import numpy as np 
import matplotlib.pyplot as plt

class Point:
    # represents a point 
    
    def __init__(self, layer_num, radius, phi, z):
        self.layer_num = layer_num  # index of which layer it is in (1~5)
        self.radius = radius        # the radius 
        self.phi = phi              # angle phi
        self.z = z                  # the thing we're interested in

class Environment: 
    
    def __init__(self, 
                 top_layer_lim:float = 100.0,       
                 bottom_layer_lim:float = 15.0, 
                 num_layers:int = 5,            # number of layers
                 radii:list = [5., 10., 15., 20., 25.]             # distance between adjacent layers
                ): 
        if top_layer_lim < bottom_layer_lim: 
            raise Exception("The top layer limits cannot be smaller than the bottom layer limits.")
        
        self.top_layer_lim = top_layer_lim
        self.bottom_layer_lim = bottom_layer_lim
        
        if len(radii) != num_layers: 
            raise Exception("The radii do not match the number of layers.")
        
        self.num_layers = num_layers 
        self.radii = sorted(radii)
        

class DataSet(): 
    
    def __init__(self, env:Environment): 
        self.env = env 
        self.array = [[] for _ in range(env.num_layers)]
        self.n_points = [0 for _ in range(env.num_layers)]
        self.total_points = 0
        
    def importData(self, data_array:list): 
        # puts a list of Point objects into DataSet structure
        
        # data array should be a (uniterated) list of Point objects 
        self.total_points = len(data_array) 
        
        for point in data_array: 
            self.array[point.layer_num - 1].append(point)
            
        ln = 0
        for layer in self.array: 
            layer.sort(key=lambda point: point.z) 
            self.n_points[ln] = len(layer)
            ln += 1
            
    def generateUniform(self, n_points:list):  # MAY BE BROKEN 
        
        if len(n_points) != len(self.array): 
            raise Exception("The n_points argument should be of form [*, *, ..., *]. ")
        else: 
            self.n_points = n_points
        
        limits_per_layer = np.linspace(self.env.bottom_layer_lim, 
                                       self.env.top_layer_lim, 
                                       self.env.num_layers+1)[1:]
        
        
        for ln in range(self.env.num_layers): 
            layer_arr = list(np.linspace(-limits_per_layer[ln], limits_per_layer[ln], self.n_points[ln]))
            points = [Point(ln, self.env.radii[ln], 0.0, z) for z in layer_arr]
            self.array[ln] = points
            
    def generateRandom(self, n_points:list):   # MAY BE BROKEN 
        
        if len(n_points) != len(self.array): 
            raise Exception("The n_points argument should be of form [*, *, ..., *]. ")
        else: 
            self.n_points = n_points
         
        
        limits_per_layer = np.linspace(self.env.bottom_layer_lim, 
                                       self.env.top_layer_lim, 
                                       self.env.num_layers+1)[1:]
        
        
        for ln in range(self.env.num_layers): 
            layer_arr = list(np.sort(np.random.uniform(low=-limits_per_layer[ln], high=limits_per_layer[ln], size=self.n_points[ln])))
            points = [Point(ln, self.env.radii[ln], 0.0, z) for z in layer_arr]
            self.array[ln] = points


    def plot(self, show_lines = False, show = False): 
        
        # Plot grey lines 
        for radius in self.env.radii: 
            plt.plot([-self.env.top_layer_lim, self.env.top_layer_lim], 
                     [radius, radius], 
                     color=(0.5, 0.5, 0.5, 0.5), 
                     linewidth=1)
            
        
        coords = [(point.z, point.radius) for layer in self.array for point in layer]
        max_height = self.env.radii[-1]
    
        plt.scatter(*zip(*coords), c="b", s=3)
        
        # X Y Labels
        plt.xlabel('z [cm]', fontsize = 16)
        plt.ylabel('r [cm]',  fontsize = 16)
        plt.yticks(np.arange(0, max_height + 1, self.env.num_layers))
        plt.title(f'Scatter Plot of Space Points', fontsize = 16)
    
        if show_lines == True: 
            
            plt.plot([self.env.bottom_layer_lim, self.env.top_layer_lim], [0.0, max_height], c="r", alpha=0.5)
            plt.plot([-self.env.bottom_layer_lim, self.env.bottom_layer_lim], [0.0, 0.0], c="r", alpha=0.5)
            plt.plot([-self.env.bottom_layer_lim, -self.env.top_layer_lim], [0.0, max_height], c="r", alpha=0.5)
        
        if show == True:
            plt.show()
        
    # Muchang - "I don't know what this is so I kept it here. "
    # def add(self, offset = 0.1):

    #     for i, value in enumerate([[-32, 32],[-49, 49],[-66, 66],[-83, 83],[-100, 100]]):
    #         phi0 = self.array[i][0].phi
    #         self.array[i].insert(0,SpacePoint(int(i+1), int((i+1)*5), phi0, value[0]-offset))
    #         self.array[i].append(SpacePoint(int(i+1), int((i+1)*5), phi0, value[1]+offset))
    #         self.n_points[i] = int(self.n_points[i]+2)
       