import numpy as np 

arr = np.linspace(1, 11, 6)
print(arr)
print(arr[1])
arr = np.delete(arr, 1)
print(arr)
print(arr[1])