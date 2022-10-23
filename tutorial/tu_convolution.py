"""
In this question, convolution is used to compute the 
laplacian of concentration function when solving the 
differential equation numerically (using Euler's method)
"""

import numpy as np
from scipy.ndimage import convolve


data = np.ones(10)
data = np.stack((data, 2*data, 3*data))
data = np.stack([data, data, data])
# Discrete Laplacian operator (Laplacian matrix)
kernel = np.array([
    [-1, -1, -1],
    [-1,  8, -1],
    [-1, -1, -1]
])

new_data = convolve(data, kernel, mode="wrap")

print(data)
print(new_data)

