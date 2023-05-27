import numpy as np

def partition(data_file, layers, pts_per_layer, radii):
    data = np.genfromtxt(data_file, delimiter = ',')
    sorted_data = np.sort(data)
    cur_min_idx = np.zeros(shape = (5, 1), dtype = int)
    limit_slope = None

    while True:
        for i in range(0, layers):
            partition_end_idx = cur_min_idx[i][0] + 15
            cur_slope = ((i + 1) * radii) / sorted_data[i][partition_end_idx]
            
            if limit_slope is None:
                limit_slope = cur_slope
            else:
                if cur_slope > 0 and limit_slope > 0:
                    limit_slope = max(limit_slope, cur_slope)
                else:
                    limit_slope = min(limit_slope, cur_slope)
                    
# partition("updated_data_1.csv", 5, 16, 5.)