This can generate a 2D (surface) or 3D (volume) fit for a data set with missing
values (represented by zero values). So best suited for matrices with
positive values only. Returns a solution spanning the size of the original matrix. 
This function performs Tikhonov Regularization (Ridge Regression) on a matrix. 
A finite difference operator is used. 

X = tikReg3D(data,smoothing) where X is the fit solution, data is the matrix with
missing (zero) values to be fit, and smoothing is the level of smoothing
in the fit, and how weighted the tikhonov matrix is. Higher smoothing values means a smoother solution. 

min(AX-b+(lambda)*TX) note: lambda == smoothing
X = ((A(A^t)+(T^t)T)^-1)(A^t)b 
^t is transposed

X == vector, smoothed solution
b == vector, input data
A == matrix A, operator