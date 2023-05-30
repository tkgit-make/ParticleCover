import numpy as np 
import matplotlib.pyplot as plt 
from updated_generate_data import *

def repeated(covers):

    """Creates covers
    
    Args:
        covers (numpy array): 3D numpy array of first set of covers
    
    Returns:
        numpy array: full set of covers
    """

    acover = np.zeros((1, 5, 16))
    loops = len(covers) - 1
    mins = []
    print(loops)

    for i in np.arange(5):
        amin = (covers[loops, i, 15]/(5*(i+1)))
        mins.append(amin)

    min_index = np.argmin(np.array(mins))
    #find closest points for all other rows
    for i in np.arange(5):
        closest_index = np.argmin(np.abs(sorted_data[i]/(5*(i+1)) - covers[loops, min_index, 15]/(min_index+1)/5))
        
        if closest_index >= 150 - 16:
            acover[0, :, 0:16] = sorted_data[:, -17:-1]
            covers = np.append(covers, acover, axis = 0)
            return covers
        
        else:
            acover[0, i, 0:16] = sorted_data[i, closest_index-1:closest_index+15] 


        
    covers = np.append(covers, acover, axis = 0)
    
    return repeated(covers)
    
    
def generate_cover(data):
    """Generates covers
    
    Args:
        covers (numpy array): 3D numpy array of sorted data
    
    Returns:
        numpy array: full set of covers
    """

    covers = np.zeros((1, 5, 16))
    xshape = np.zeros_like(data[0, 42:58])
                    
    #generate first cover
    for row in np.arange(5):
        #min_index = np.argmin(np.abs(sorted_data[row]))
        covers[0, row, 0:16] = data[row, 0:16]
        #repeated(data, covers)

    return repeated(covers)


if __name__ == "__main__":
	naive_data = generate_naive_data(l_bound = -1., h_bound = 1., r0_xmin = -0.15, r0_xmax = 0.15, layers = 5, radii = 5., size = 150)
	plot(naive_data, radii = 5.0)
	sorted_data = np.sort(naive_data)

	colorset = ['orange', 'k', 'b', 'r', 'y', 'm']
	final = generate_cover(sorted_data)
	for j in range(len(final)):
		for i in range(5):
			plt.scatter(final[j, i, :], np.full_like(final[j, i, :], (i+1)*5), s = 10, color = colorset[j%6])
	plt.xlim(-1, 1)
	plt.ylim(0, 26)