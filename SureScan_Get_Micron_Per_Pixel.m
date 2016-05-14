function d = SureScan_Get_Micron_Per_Pixel(Im)

% clc;
fig1 = figure;

imshow(Im);
% imshow(fliplr(Im));
set(fig1,'units','normalized','outerposition',[0 0 1 1]);
title('Zoom in the figure and press any key...', 'Fontsize',14);
display('Zoom in the figure and press any key...');
pause;

[x,y] = ginput(2);
close gcf;

d  = sqrt ( (x(1)-x(2))^2 + (y(1)-y(2))^2 );

prompt = {'Enter distance in cm:'};
dlg_title = 'Input';
numlines = 1;
defaultanswer = {'5'};
distance = inputdlg(prompt, dlg_title, numlines, defaultanswer);



d = str2num(distance{1}) * 10000 / d;
end