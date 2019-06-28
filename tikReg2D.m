function [X,A,T] = tikReg2D(slice,lambda)
% tikReg2D() generates a surface to fit to the data in "slice". Based on
% John D'Errico's Gridfit. Zeros are ignored from the fit, tikhonov
% regularization allows for fitting a surface over data with large holes or
% missing data.
%
% slice == data to fit a surface over (will fit to non-zero data in the
% array), 2D matrix
% lambda == smoothing paramter, the higher the value, the smoother the fit
%
% W. Quinn Meadus, June 2019

s = size(slice);
ny = s(1);
nx = s(2);

%Tikhonov regularizer
%min(Ax-b+(lambda)*Tx)

b = slice(slice(:)>0); %rhs data, assuming 0 values are to be excluded from the fit
bind = find(slice(:)); %rhs location in full grid

nb = length(b);
ngrid = length(slice(:));

%Holds the information for the location of each b value in the full grid
%(bInd) while having a row corresponding to each b value.
A = sparse((1:nb)',bind, ones(nb,1),nb,ngrid);


%difference approximation in y
[i,j] = meshgrid(1:nx,2:(ny-1));
ind = j(:) + ny*(i(:)-1);
len = length(ind);

T2 = sparse(repmat(ind,1,3), [ind-1,ind,ind+1], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);

%difference approximation in x
[i,j] = meshgrid(2:(nx-1),1:ny);
ind = j(:) + ny*(i(:)-1);
len = length(ind);

T1 = sparse(repmat(ind,1,3), [ind-ny,ind,ind+ny], [-1*ones(len,1),2*ones(len,1),-1*ones(len,1)], ngrid,ngrid);

%Combining regularization (tikhonov) matrices
T = [T1;T2];

%appending zeros to the rhs
b = [b;zeros(size(T,1),1)];

%solving the minimization problem (tikhonov regularization solution)
AT = [A;lambda*T];
X = reshape((AT'*AT)\(AT'*b),ny,nx);




