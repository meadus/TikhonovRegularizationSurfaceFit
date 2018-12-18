function X = tikReg3D(data,smoothing)
% This can generate a 3D volume which fits to a data set with missing
% values (represented by zero values). So best suited for matrices with
% positive values only. Returns a solution spanning the size of the
% original matrix. 
% This function performs Tikhonov Regularization (Ridge Regression) on a 3D matrix. 
% A finite difference operator is used.

% X = tikReg3D(data,smoothing) where X is the fit solution, data is the matrix with
% missing (zero) values to be fit, and smoothing is the level of smoothing
% in the fit, and how weighted the tikhonov matrix is. Higher smoothing values means a smoother solution. 

% min(AX-b+(lambda)*TX) note: lambda == smoothing
% X = ((A(A^t)+(T^t)T)^-1)(A^t)b 
% ^t is transposed

% X == vector, smoothed solution
% b == vector, input data
% A == matrix A, operator

s = size(data);
ny = s(1);
nx = s(2);
nz = s(3);

b = data(data(:)>0); %rhs data, assuming 0 values are to be excluded from the fit
bind = find(data(:)); %rhs location in full grid

nb = length(b);
ngrid = length(data(:));

%Holds the information for the location of each b value in the full grid
%(bInd) while having a row corresponding to each b value.
A = sparse((1:nb)',bind, ones(nb,1),nb,ngrid);


%Generating the Tikhonov Matrix:
%Uses finite difference approximations to enforce smoothness in the
%solution

%difference approximation in y
[i,j,k] = meshgrid(1:nx,2:(ny-1),1:nz);
ind = j(:) + ny*(i(:)-1)+ny*nx*(k(:)-1);
len = length(ind);

T1 = sparse(repmat(ind,1,3), [ind-1,ind,ind+1], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);

%difference approximation in x
[i,j,k] = meshgrid(2:(nx-1),1:ny,1:nz);
ind = j(:) + ny*(i(:)-1)+ny*nx*(k(:)-1);
len = length(ind);

T2 = sparse(repmat(ind,1,3), [ind-ny,ind,ind+ny], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);

%difference approximation in z
[i,j,k] = meshgrid(1:nx,1:ny,2:(nz-1));
ind = j(:) + ny*(i(:)-1)+ny*nx*(k(:)-1);
len = length(ind);

T3 = sparse(repmat(ind,1,3), [ind-ny*nx,ind,ind+ny*nx], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);


%Combining regularization (tikhonov) matrices
T = [T1;T2;T3];

%appending zeros to the rhs
b = [b;zeros(size(T,1),1)];

%solving the minimization problem (tikhonov regularization solution)
clear T1 T2 T3
X = reshape(([A;smoothing*T]'*[A;smoothing*T])\([A;smoothing*T]'*b),ny,nx,nz); %Here matlab's \ solver is used to solve the minimization problem.

end



