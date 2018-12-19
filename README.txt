X == fit solution
data == matrix with missing (zero) values to be fit
smoothing == weight of tikhonov matrix, higher smoothing values means a flatter/smoothly varying solution 

This can generate a 2D (surface) or 3D (volume) fit for a data set with missing
values (represented by zero values). So best suited for matrices with
positive values only. Returns a solution spanning the size of the original matrix. 
Tikhonov regularization (ridge regression) with a finite difference operator is used. 

min(AX-b+(lambda)*TX) note: lambda == smoothing
X = ((A(A^t)+(T^t)T)^-1)(A^t)b 
^t is transposed

X == vector, smoothed solution
b == vector, input data
A == matrix A, operator