function [VolumeED,VolumeES] = results(Im,dFit, xED,yED, sFit, xES, yES, H_D, H_S)
figure, imshow(Im(:,:,:,1)); hold on; plot(xED,yED, 'c*');
[xpick,ypick] = ginput(2); close gcf;

clear xED yED;
xED = round(xpick(1)):round(xpick(2));
for p = 1:length(xED)
    x = xED(p);
    yED(p) = dFit.p1*x^5 + dFit.p2*x^4 + dFit.p3*x^3 + dFit.p4*x^2 + dFit.p5*x + dFit.p6;
end
if ( yED(1) < ypick(1) )
    yED = [ypick(1):-1:yED(1)  yED];
    xED = [xED(1)*ones(1,(length(yED)-length(xED)))  xED];
end
if ( yED(end) < ypick(2) )
    yED = [yED yED(end):ypick(2)];
    xED = [xED xED(end)*ones(1,(length(yED)-length(xED))) ];
end
index = find( yED > ypick(1) ); yED(index) = []; xED(index) = [];
index = find( yED > ypick(2) ); yED(index) = []; xED(index) = [];

[xED,yED, VolumeED] = get_volume_LV(xED',yED');

clear xES yES;
xES = round(xpick(1)):round(xpick(2));
for p = 1:length(xES)
    x = xES(p);
    yES(p) = sFit.p1*x^5 + sFit.p2*x^4 + sFit.p3*x^3 + sFit.p4*x^2 + sFit.p5*x + sFit.p6;
end

if ( yES(1) < ypick(1) )
    yES = [ypick(1):-1:yES(1)  yES];
    xES = [xES(1)*ones(1,(length(yES)-length(xES)))  xES];
end
if ( yES(end) < ypick(2) )
    yES = [yES yES(end):ypick(2)];
    xES = [xES xES(end)*ones(1,(length(yES)-length(xES))) ];
end
index = find( yES > ypick(1) ); yES(index) = []; xES(index) = [];
index = find( yES > ypick(2) ); yES(index) = []; xES(index) = [];

[xES,yES, VolumeES] = get_volume_LV(xES',yES');

gold_ED1 = H_D{1};
gold_ED2 = H_D{2};
gold_ES1 = H_S{1};
gold_ES2 = H_S{2};


close all
for fr = 1:50 %size(Im,4)
    hold off;
    imshow(Im(:,:,:,fr));
    hold on;
    plot(xED,yED,'b.');
    plot(xES,yES,'r.');
    plot(gold_ED1(:,1),gold_ED1(:,2),'y.');
    plot(gold_ED2(:,1),gold_ED2(:,2),'y.');
    plot(gold_ES1(:,1),gold_ES1(:,2),'c.');
    plot(gold_ES2(:,1),gold_ES2(:,2),'c.');
    pause(0.01);
end


clc
[hd, ~, ~] = SureScan_Hausdorff_Distance(gold_ED1, gold_ED2)
[hd, ~, ~] = SureScan_Hausdorff_Distance(gold_ED1, [xED yED])
[hd, ~, ~] = SureScan_Hausdorff_Distance(gold_ED2, [xED yED])


[hd, ~, ~] = SureScan_Hausdorff_Distance(gold_ES1, gold_ES2)
[hd, ~, ~] = SureScan_Hausdorff_Distance(gold_ES1, [xES yES])
[hd, ~, ~] = SureScan_Hausdorff_Distance(gold_ES2, [xES yES])

VolumeED
VolumeES
end