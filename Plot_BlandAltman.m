function [meanDiff, stdDiff] = Plot_BlandAltman(A,B, PlotType, Xlabel, Ylabel)

meanAB = (A+B)./2;
difff = (A-B); % ./ (meanAB/100);
meanDiff = mean(difff);
stdDiff = std(difff);

meanp2D = meanDiff+2*stdDiff;
meanm2D = meanDiff-2*stdDiff;
n = length(difff);
minD = min(0,min(meanAB));
maxD = max(meanAB)*1.1;

if strcmp( PlotType , 'full' )
    fig = figure; set(fig,'units','normalized','outerposition',[0 0 1 1]);
end
plot(meanAB,difff,'O','MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5)



% scatter(meanAB,difff,'filled','SizeData',20, 'MarkerFaceColor', 'r');
hold on;
line([minD; maxD],[meanp2D meanp2D], 'Color', 'b','LineStyle','-','LineWidth',1);
text(minD+0.05,meanp2D+0.5,'Mean + 2SD','fontsize',12);
hold on;
line([minD; maxD],[meanm2D meanm2D], 'Color', 'b','LineStyle','-','LineWidth',1);
text(minD+0.05,meanm2D-0.5,'Mean - 2SD','fontsize',12);
hold on;
line([minD; maxD],[meanDiff meanDiff], 'Color', 'k','LineStyle','-','LineWidth',1);
line([minD; maxD],[0        0       ], 'Color', 'k','LineStyle','--','LineWidth',1);

if nargin > 4
    xlabel(['( ' Xlabel ' + ' Ylabel ' )/2'],'fontsize',12, 'FontWeight', 'Bold');
    ylabel([Xlabel ' - ' Ylabel],'fontsize',12, 'FontWeight', 'Bold');
end

set(gca,'FontWeight','bold')


% [~,~,ci,~] = ttest(difff,0,'Alpha',0.05)
end