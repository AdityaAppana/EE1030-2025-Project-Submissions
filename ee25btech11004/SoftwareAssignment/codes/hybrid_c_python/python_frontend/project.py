import numpy as np
from matplotlib.image import imread
import ctypes
import matplotlib.pyplot as plt

A = imread('/Users/adityaappana/Desktop/IIT Hyderabad/Sem1 EE29/EE1030 - Matrix Theory/ee1030-2025/ee25btech11004/SoftwareProject/SoftwareProject/greyscale.png')
X = np.mean(A,-1);

X = X.astype(np.double)
rows,cols = X.shape

blank = np.zeros((rows,cols), dtype = np.double)


c_lib = ctypes.CDLL('./svdshared.so')

# 2. Define the C type for a pointer to a double
c_double_ptr = ctypes.POINTER(ctypes.c_double)

# 3. Create a C-style array of pointers to doubles
# This will be equivalent to 'double**' in C


c_lib.svd.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags="C_CONTIGUOUS"), # float** og
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags="C_CONTIGUOUS"), # float** sum
    ctypes.c_int,   # int rows
    ctypes.c_int,   # int columns
    ctypes.c_int    # int rank
]

rank = int(input(f"Enter rank for approximation"))


c_lib.svd.restype = None


returnarrr = c_lib.svd(X,blank,rows,cols,rank)

img = plt.imshow(blank)
img.set_cmap('gray')
plt.show()



