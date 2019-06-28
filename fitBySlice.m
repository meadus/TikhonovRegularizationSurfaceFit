function corr = fitBySlice(img,normalization_mask,smoo)
% corr = fitBySlice(img,normalization_mask,smoo) generates a smoothed
% surface for each slice of a 3D volume. Divide the original image by this
% to 'normalize' the image. 
%
% img == lung image (3D array) to be normalized
% normalization_mask == logical array defining what pixels are to be
% included in the fit (generally muscle tissue and liver, excludes lungs
% and empty space/ dark regions)
% smoo == smothing parameter (for tikhonov regularization)
%
% W. Quinn Meadus, June 2019

data = img.*normalization_mask;

s = size(data);
corr = zeros(s);

for i = 1:s(3)
    corr(:,:,i) = tikReg2D(data(:,:,i),smoo);
end


end