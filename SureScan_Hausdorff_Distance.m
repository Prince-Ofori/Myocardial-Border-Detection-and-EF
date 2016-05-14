function [hd, ind1, ind2] = SureScan_Hausdorff_Distance(pts1, pts2)
%HAUSDORFFDISTANCE  Hausdorff distance between two point sets
%
%   HD = hausdorffDistance(PTS1, PTS2)
%   Computes the Hausdorff distance between the two point sets PTS1 and
%   PTS2. The Hausdorf distance can be used to compare two shapes. 
%
%   The distance between a point x and a set Y is given by:
%     d(x, Y) = inf { d(x,y) | y in Y }
%   The distance between two non empty sets X and Y is given by:
%     d(X, Y) = sup { d(x,Y) | x in X }
%   The Hausdorff distance between sets X and Y distance is defined as the
%   maximum of d(X,Y) and d(Y,X):
%     HD(X,Y) = max { d(X,Y), d(Y,X) }
%
%
%   Example
%   % Compute Hausdorff distance between an ellipse and a rectangle
%     % first define two shapes
%     rect = resamplePolygon(orientedBoxToPolygon([20 30 80 40 30]), 60);
%     poly = ellipseToPolygon([20 30 40 20 30], 500);
%     % display the shapes
%     figure; hold on
%     drawPolygon(poly, 'b');
%     drawPolygon(rect, 'g');
%     axis equal;
%     % compute hausdorff distance
%     [hd ind1 ind2] = hausdorffDistance(poly, rect);
%     p1h = poly(ind1, :);
%     p2h = rect(ind2, :);
%     drawPoint([p1h;p2h], 'mo');
%     drawEdge([p1h p2h], 'm')
%
%   See also
%   minDistancePoints
%
%   References
%   http://en.wikipedia.org/wiki/Hausdorff_distance
%
% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2012-05-04,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2012 INRA - Cepia Software Platform.

% distance from pts1 to pts2
[dists1, ind12] = minDistancePoints(pts1, pts2);
[max1, ind11] = max(dists1);

% distance from pts2 to pts1
[dists2, ind22] = minDistancePoints(pts2, pts1);
[max2, ind21] = max(dists2);

% keep the max of the two distances
hd = max(max1, max2);

% keep the rigt indices
if max1 > max2
    ind1 = ind11;
    ind2 = ind12(ind11);
else
    ind1 = ind22(ind21);
    ind2 = ind21;
end

end





function varargout = minDistancePoints(p1, varargin)
%MINDISTANCEPOINTS Minimal distance between several points
%
%   DIST = minDistancePoints(PTS)
%   Returns the minimum distance between all couple of points in PTS. PTS
%   is a N-by-D array of values, N being the number of points and D the
%   dimension of the points.
%
%   DIST = minDistancePoints(PTS1, PTS2)
%   Computes for each point in PTS1 the minimal distance to every point of
%   PTS2. PTS1 and PTS2 are N-by-D arrays, where N is the number of points,
%   and D is the dimension. Dimension must be the same for both arrays, but
%   number of points can be different.
%   The result is an array the same length as PTS1.
%
%
%   DIST = minDistancePoints(..., NORM)
%   Uses a user-specified norm. NORM=2 means euclidean norm (the default), 
%   NORM=1 is the Manhattan (or "taxi-driver") distance.
%   Increasing NORM growing up reduces the minimal distance, with a limit
%   to the biggest coordinate difference among dimensions. 
%   
%
%   [DIST I J] = minDistancePoints(PTS)
%   Returns indices I and J of the 2 points which are the closest. DIST
%   verifies relation:
%   DIST = distancePoints(PTS(I,:), PTS(J,:));
%
%   [DIST J] = minDistancePoints(PTS1, PTS2, ...)
%   Also returns the indices of points which are the closest. J has the
%   same size as DIST. It verifies relation: 
%   DIST(I) = distancePoints(PTS1(I,:), PTS2(J,:));
%   for I comprised between 1 and the number of rows in PTS1.
%
%
%   Examples:
%   % minimal distance between random planar points
%       points = rand(20,2)*100;
%       minDist = minDistancePoints(points);
%
%   % minimal distance between random space points
%       points = rand(30,3)*100;
%       [minDist ind1 ind2] = minDistancePoints(points);
%       minDist
%       distancePoints(points(ind1, :), points(ind2, :))
%   % results should be the same
%
%   % minimal distance between 2 sets of points
%       points1 = rand(30,2)*100;
%       points2 = rand(30,2)*100;
%       [minDists inds] = minDistancePoints(points1, points2);
%       minDists(10)
%       distancePoints(points1(10, :), points2(inds(10), :))
%   % results should be the same
%   
%   See Also
%   points2d, distancePoints, nndist, hausdorffDistance
%
%   ---------
%   author: David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 15/06/2004.
%

% HISTORY:
% 22/06/2005 compute sqrt only at the end (faster), and change behaviour
%   for 2 inputs: compute min distance for each point in PTS1. 
%   Also add support for different norms.
% 15/08/2005 make difference when 1 array or 2 arrays of points
% 25/10/2006 also returns indices of closest points
% 30/10/2006 generalize to points of any dimension
% 28/08/2007 code cleanup, add comments and help


%% Initialisations

% default norm (euclidean)
n = 2;

% flag for processing of all points
allPoints = false;

% process input variables
if isempty(varargin)
    % specify only one array of points, not the norm
    p2 = p1;
    
elseif length(varargin) == 1
    var = varargin{1};
    if length(var) > 1       
        % specify two arrays of points
        p2  = var;
        allPoints = true;
    else
        % specify array of points and the norm
        n   = var;
        p2  = p1;
    end
    
else
    % specify two array of points and the norm
    p2  = varargin{1};
    n   = varargin{2};
    allPoints = true;
end


% number of points in each array
n1  = size(p1, 1);
n2  = size(p2, 1);

% dimension of points
d   = size(p1, 2);


%% Computation of distances

% allocate memory
dist = zeros(n1, n2);

% different behaviour depending on the norm used
if n == 2
    % Compute euclidian distance (default case).
    % Compute difference of coordinate for each pair of point and for each
    % dimension. Result "dist" is a n1-by-n2 array. 
    % in 2D: dist = dx.*dx + dy.*dy;
    for i = 1:d
        dist = dist + (repmat(p1(:,i), [1 n2])-repmat(p2(:,i)', [n1 1])).^2;
    end
    
    % compute minimal distance:
    if ~allPoints
        % either on all couple of points
        mat = repmat((1:n1)', [1 n1]);
        ind = mat < mat';
        [minSqDist, ind] = min(dist(ind));
    else
        % or for each point of P1
        [minSqDist, ind] = min(dist, [], 2);
    end
    
    % convert squared distance to distance
    minDist = sqrt(minSqDist);
    
elseif n == inf
    % infinite norm corresponds to maximum absolute value of differences
    % in 2D: dist = max(abs(dx) + max(abs(dy));
    for i = 1:d
        dist = max(dist, abs(p1(:,i)-p2(:,i)));
    end
    
else
    % compute distance using the specified norm.
    % in 2D: dist = power(abs(dx), n) + power(abs(dy), n);
    for i = 1:d
        dist = dist + power((abs(repmat(p1(:,i), [1 n2])-repmat(p2(:,i)', [n1 1]))), n);
    end

    % compute minimal distance
    if ~allPoints
        % either on all couple of points
        mat = repmat((1:n1)', [1 n1]);
        ind = mat < mat';
        [minSqDist, ind] = min(dist(ind));
    else
        % or for each point of P1
        [minSqDist, ind] = min(dist, [], 2);
    end

    % convert squared distance to distance
    minDist = power(minSqDist, 1/n);
    
end

if ~allPoints
    % convert index in array to row and column subindices.
    % This uses the fact that index are sorted in a triangular matrix,
    % with the last index of each column being a so-called triangular
    % number
    ind2 = ceil((-1+sqrt(8*ind+1))/2);
    ind1 = ind - ind2*(ind2-1)/2;
    ind2 = ind2 + 1;
end


%% format output parameters

% format output depending on number of asked parameters
if nargout <= 1
    varargout{1} = minDist;
    
elseif nargout == 2
    % If two arrays are asked, 'ind' is an array of indices, one for each
    % point in PTS1, corresponding to the result in minDist
    varargout{1} = minDist;
    varargout{2} = ind;
    
elseif nargout == 3
    % If only one array is asked, minDist is a scalar, ind1 and ind2 are 2
    % indices corresponding to the closest points.
    varargout{1} = minDist;
    varargout{2} = ind1;
    varargout{3} = ind2;
end

end