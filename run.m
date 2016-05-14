close all; 
% MAIN 
% Call Image Fusion function and return fused DICOM sequence
imFuse = imageFusion(Im);
% Input x,y coordinates of the region of interest
imshow(imFuse);
[x,y] = ginput(4);
x = round(x); y = round(y); % Round user input coordinates 
roiREC = round([min(x), min(y), max(x)-min(x), max(y)-min(y)]);
inputIm = imcrop(imFuse,roiREC);

% UNCOMMENT TO TEST
%   subplot(1,2,1)
%   imshow(imFuse),title('Original Image')
%   subplot(1,2,2)
%   imshow(inputIm),title('Region of Interest');

% Call region growing function 
seedReg = regionGrow(inputIm,round(0.5*size(inputIm,2)), round(size(inputIm,1)),17);
% Return boundary of the binary cardial walls 
[B] = bwboundaries(seedReg,4);

% Call Peak Boundary and return curve fits of diastole and systole with
% extracted x and y coordinates of fits
[dFit, sFit,bX, bY,xp,yp] = peakBoundary(imFuse,inputIm, B, x, y);

% Call results to return end-diatolic volume and end-systolic volume and
% pertenant testing results
[volED, volES] = results(Im,dFit, xp,yp, sFit, bX, bY, H_D, H_S);
clc; 

disp('Analysing...'); disp('=====================================')
% Compute ejection fraction
sv = volED - volES; % Define stroke volume
EF = (sv/volED) * 100;

% Normalised EF values in cm
scale = 0.0431; % Scale defined after using SureScan_Get_Micron_Per_Pixel
EDVcm = volED * scale^3;
ESVcm = volES * scale^3;
% Build string and display on console
pDVol = ['End-Diastolic Volume: ',num2str(EDVcm),'cm']; disp(pDVol);
pSVol = ['End-Systolic Volume: ',num2str(ESVcm),'cm']; disp(pSVol);
pEF = ['Ejection Fraction: ',num2str(EF),'%']; disp(pEF);

% UNCOMMENT TO DEMONSTRATE 
% n = 2;
% while n ~= 0 || n ~= 1
%     n = input('Enter 0 to view Diastolic Fit or enter 1 for Systolic Fit: ');
%     switch n
%         case 0
%             for i = 1:size(Im(:,:,:,:))
%                 clf;
%                 imshow(Im(:,:,:,i)), title(num2str(i)); hold on;
%                 plot(dFit);
%                 ax = gca; legend(ax,'off');
%                 colormap gray;
%                 for ii = 1:length(B)
%                     boundary = B{ii};
%                     plot(boundary(:,2) + min([x(1) x(2)]), boundary(:,1) + min(y(3), y(4)), 'g', 'LineWidth', 2);
%                     %line([x(1),x(2)], [y(1),y(2)]);
%                 end
%                 pause(0.01);
%             end
%             disp(volume);
%         case 1
%             for i = 1:size(Im(:,:,:,:))
%                 clf;
%                 imshow(Im(:,:,:,i)), title(num2str(i)); hold on;
%                 plot(sFit);
%                 ax = gca; legend(ax,'off');
%                 colormap gray;
%                 pause(0.01);
%             end
%         otherwise
%             disp('Incorrect value input');
%     end
%   end