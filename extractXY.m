function [xp, yp] = extractXY(M)
% EXTRACTXY returns the X and Y values returned from the PEAKBOUNDARY function's 
% output cell array

% Preallocate memory for the x and y coordinates 
xp = []; yp = [];
% Iterate through the values of M and return values of X and Y 
for c = 1:length(M)
    xy = M{c};
    xp = [xp xy(1)];
    yp = [yp xy(2)];
end
end