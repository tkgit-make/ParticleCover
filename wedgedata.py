from data import *
from spaceptcollection import *

import numpy as np

class WedgeData(DataSet):
    def __init__(self, env:Environment, n_points, data_arr):
        self.env = env
        self.n_points = n_points
        self.array = data_arr # data_arr: nested list of respective SpacePoints per layer for self.env.layers
    
    def plot(self, show_lines = False, show = True):
        z_arr = [pt.z for pt in np.hstack(self.array)]
        r_arr = np.array([])

        for i in range(0, self.env.layers):
            r_arr = np.concatenate([r_arr, (i + 1) * 5 * np.ones(self.n_points[i])])
        
        plt.scatter(z_arr, r_arr, c = "g", s = 2)

        max_radius = r_arr[-1]
        plt.yticks(np.arange(0, max_radius + 1, self.env.layers))
        plt.xlabel('z [cm]', fontsize = 16)
        plt.ylabel('r [cm]',  fontsize = 16)
        plt.title(f'Wedge Data for v3 Data', fontsize = 16)

        if show_lines:
            plt.plot([self.env.bottom_layer_lim, self.env.top_layer_lim], [0, max_radius], c="r", alpha=0.5)
            plt.plot([-self.env.bottom_layer_lim, self.env.bottom_layer_lim], [0.0, 0.0], c="r", alpha=0.5)
            plt.plot([-self.env.bottom_layer_lim, -self.env.top_layer_lim], [0, max_radius], c="r", alpha=0.5)
        if show == True:
            plt.show()


    def add(self, offset = 0.1):

        for i, value in enumerate([[-32, 32],[-49, 49],[-66, 66],[-83, 83],[-100, 100]]):
            phi0 = self.array[i][0].phi
            self.array[i].insert(0,SpacePoint(int(i+1), int((i+1)*5), phi0, value[0]-offset))
            self.array[i].append(SpacePoint(int(i+1), int((i+1)*5), phi0, value[1]+offset))
            self.n_points[i] = int(self.n_points[i]+2)

