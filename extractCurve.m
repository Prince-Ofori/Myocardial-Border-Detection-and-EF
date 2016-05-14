function [truX, truY] = extractCurve(fit, x, y)
% EXTRACTCURVE returns the X and Y values within a cfit object 
%figure; h = plot( fit, x,y); hold on

for p = 1:length(x)
    curveX = x(p);
    curveY = fit.p1*curveX^5 + fit.p2*curveX^4 + fit.p3*curveX^3 + fit.p4*curveX^2 + fit.p5*curveX + fit.p6;
    curveXY = [curveX, curveY];
    coord(p) = {curveXY};
    %plot(curveX,curveY,'g.')
end
% Extract coordinantes within the cell array 
truX = []; truY = [];
for c = 1:length(coord)
    xy = coord{c};
    truX = [truX xy(1)];
    truY = [truY xy(2)];
end
% Round output true x and y coordinates 
truX = round(truX); truY = round(truY);
end