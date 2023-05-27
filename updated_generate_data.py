import numpy as np
import matplotlib.pyplot as plt

def generate_naive_data(l_bound = -1., h_bound = 1., r0_xmin = -0.15, r0_xmax = 0.15, layers = 5, radii = 5., size = 100):
    n_slope = (radii * layers) / (l_bound - r0_xmin)
    p_slope = (radii * layers) / (h_bound - r0_xmax)

    fin_arr = np.array([])

    for i in range(1, layers + 1):
        cur_y_coord = radii * i

        cur_xmin = (cur_y_coord - r0_xmin) / n_slope
        cur_xmax = (cur_y_coord - r0_xmax) / p_slope

        arr = np.random.uniform(cur_xmin, cur_xmax, size)

        if i == 1:
            fin_arr = arr
        else:
            fin_arr = np.vstack((fin_arr, arr))

    return fin_arr

def plot(arr, radii = 5.0):
    num_pnts = 5 * np.ones(arr.shape[1])
    
    for i in range(arr.shape[0]):
        plt.scatter(arr[i], num_pnts, s=2, c="b")
        num_pnts += 5 * np.ones(arr.shape[1])
    
    plt.yticks(np.arange(0, num_pnts[0] + 1, radii))
    
    plt.show()

naive_data = generate_naive_data(l_bound = -1., h_bound = 1., r0_xmin = -0.15, r0_xmax = 0.15, layers = 5, radii = 5., size = 100)

plot(naive_data, radii = 5.0)