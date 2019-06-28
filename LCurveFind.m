function [spf,ind,xL,yL] = LCurveFind(slice)
% LCurveFind() finds a smoothing parameter to balance a tikhonov regularization problem
% based on the L-curve method as described in the 1992 paper. This function
% was used to determine a proper smoothing value for the lung water image
% set normalization. Effectively finds a balacance between error from
% fitting to the original data and error from oversmoothing.
%
% [spf,ind,xL,yL] = LCurveFind(slice)
% slice == slice to be normalized, zero values are excluded from the fit
% spf == final smoothing parameter (lambda)
% ind === index of array of possible smoothing parameters
% xL, yL == x and y for the L-curve
% plot(xL,yL,xL(ind),yL(ind),'*') to see what point of the curve was chosen
%
% W. Quinn Meadus, June 2019

slice = abs(slice);
ind = slice>0; %values less than zero will be excluded from the fitting process
b = slice(ind);

%Used a nonlinear array of smoothing values to test
%Each of these has to be tested so having too many makes for an extremely
%slow process. These numbers just happened to fit my problem.
s0 = 0.5:0.25:19;
s1 = 20:10:490;
s2 = 500:100:5000;
s3 = 6000:1000:20000;
sp = [s0,s1,s2,s3];

for i = 1:length(sp)
    %Fitting the surface with a specific smoothing parameter
    [g,A,T] = tikReg2D(slice,sp(i));
    X = g(:);
    
    %calculating the errors with that parameter
    rn(i) = norm(A*X-b,'fro');
    ssn(i) = norm(T*X,'fro');
end

%interpolating the curve
xx = 0.5:0.25:50000;
xL = log(rn);
xL = spline(sp,xL,xx);
yL = log(ssn);
yL = spline(sp,yL,xx);

%calculating maximum curvature, derivatives
xL1 = gradient(xL);
yL1 = gradient(yL);

xL2 = del2(xL);
yL2 = del2(yL);

k = (xL2.*yL1-xL1.*yL2)./(xL1.^2+yL1.^2).^1.5; %final curvature equations
[~,ind] = min(k);


spf = xx(ind); %final smoothing parameter (lambda) at the point of maximum (negative) curvature.
end







