function [imFuse] = imageFusion(Im)
% IMAGEFUSION utilises image fusion to develop an image with the
% information present in all frames within the sequence

% Size of the image passed to the imFuse matrix
sz = size(Im(:,:,:));
imFuse = zeros(sz(1),sz(2));
% Merging images through all frames
for i = 1:size(Im(:,:,:))
    iFrame = Im(:,:,i);
    iFrame = im2double(iFrame);
    imFuse = imFuse + iFrame; 
    % UNCOMMENT TO TEST
    %   imagesc(imFuse);
    %   colormap gray; 
    %   title(num2str(i));
    %   colorbar
    %   pause(); 
end
% Converting double matrix to Uint8 matrix to only permit integer value and
% filtering to remaining noise
imFuse = uint8(255 * mat2gray(imFuse));
imFuse = medfilt2(imFuse);
end