function [xx, yy] = lineExtend(startP, endP)
% Extend line through start and end points 

lnGradient = (endP(2)-startP(2))/(endP(1)-startP(1));
limit = 10000;
if(abs(lnGradient) == Inf) % Vertical line
    yy = -limit:limit; % Line limit between -10,000 to +10,000
    xx = startP(1)*ones(1,length(yy));
elseif(lnGradient == 0) % Horizontal Line
    xx = -limit:limit;
    yy = startP(2)*ones(1,length(xx));
else %Computes diagonal line 
    xx = -limit:limit;
    b = startP(2) - (lnGradient * startP(1));
    yy = lnGradient*xx+b;
end
end