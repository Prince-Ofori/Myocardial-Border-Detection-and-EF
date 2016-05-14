function [dFit,sFit,X,Y,xp,yp] = peakBoundary(imGlobal,imLocal, B, x, y)
close all;
%PEAKBOUNDARY Computes the a secondary cardial boundary based on peak
%intensity from the ventricular wall
%   Utilising the previously defined boundary, this function computes a
%   secondary boundary by finding the highest intensity region.
X = []; Y = [];
szLoc = size(imLocal);
% Origin defines the centroid of the line extension operation through the center x of the local image 
origin = [0.5*size(imLocal,2), round(size(imLocal,1))];

% Define the (X,Y) values of the boundary
for i = 1:length(B)
    boundary = B{i};
    X = [X; boundary(:,2)];
    Y = [Y; boundary(:,1)];
end
% Plot line extensions from the centroid of the local image through
% boundary points
figure,imshow(imGlobal),hold on;
for p = 1:length(X)
    [xx, yy] = lineExtend([origin(1), origin(2)],[X(p), Y(p)]);
    
    index = find(yy > origin(2));
    yy(index) = [];xx(index) = [];

    
    % Find spurious boundary values within line plots and exclude them
    % from final line plot
    index = find(yy > Y(p));
    xx(index) = []; yy(index) = [];
    
    index = find(xx > szLoc(2));
    xx(index) = [];yy(index) = [];
    
    index = find(xx < 0);
    xx(index) = [];yy(index) = [];
    
    index = find(yy < 0);
    xx(index) = [];yy(index) = [];
    
    index = find(yy > szLoc(1));
    xx(index) = [];yy(index) = [];
    xx = round(xx);yy = round(yy);
    
    % Find elements index elements smaller than 1 to and setting them to a
    % valid index value
    indices = find(xx < 1); xx(indices) = 1;
    indices = find(yy < 1); yy(indices) = 1;   
    pixels = [];
    for px = 1:length(xx)
        pixels = [pixels imLocal(yy(px),xx(px))];
    end
    % Return the indices of the highest intensity point within the line
    % plot 
    [intensity, index] = max(pixels);
    %figure,imshow(imGlobal), hold on; 
    M = [xx(index) + x(1), yy(index) + y(3)];
    mF(p) = {M};
    
    % UNCOMMENT TO TEST
    %   plot(xx + x(1), yy + y(3), 'c.-');
    %   plot(X(p) + x(1), Y(p) + y(3), 'b*');
    %   plot(origin(1) + x(1), origin(2) + y(3), 'r*');
    %   plot(M(1),M(2), 'r*');
    %   plot(X + x(1),Y + y(3), 'c.');
    %   figure,stem(pixels);
    %   pause(); %close gcf;
    
end
[xp, yp] = extractXY(mF);
%   [volX, VolY] = get_volume(xp, yp);
[dFit] = curveFit(xp, yp);
[sFit] = curveFit(X + x(1), Y + y(3));
end