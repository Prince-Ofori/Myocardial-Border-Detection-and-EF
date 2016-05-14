function seed = regionGrow(i,x,y,tolerance )
%REGIONGROW function applies a single seed region growth using a
%thresholding value within tolerance parameter
%mean based on the tolerance value passed to the function

%If a tolerance parameter is not passed set tolerance to 1
if exist('tolerance', 'var')== 0
    tolerance = 1;
end
%If x and y coordinates of the seed pixel are not passed get input from the
%user and round the nearest interger 
if (x == 0 || y == 0 )
    imshow(i);
    [x,y] = ginput(1);
    %Round x and y to the nearest integer
    x = round(x);
    y = round(y);
end
%Create seed point within a logical matrix(Logical Indexing)
seed = false(size(i,1),size(i,2));
seedIn = seed;
seed(uint8(y),uint8(x)) = 1;
while(sum(seed(:)) ~= sum(seedIn(:)))
    seedIn = seed;
    %Evaluate intensity at seed point and calculate region mean
    seg = i(seed);
    % UNCOMMENT TO TEST
    %   subplot(1,2,1); imagesc(seed);
    %   subplot(1,2,2); imagesc(seg); pause(0.1);
    meanSeg = mean(seg);
    %Grow region using a disk morphological growth
    growSeed = imdilate(seed,strel('disk',1,0)) - seed;
    %Evaluate new intensity within the growSeed region
    intenseSearch = find(growSeed);
    newIntense = i(intenseSearch);
    seed(intenseSearch(newIntense > meanSeg - tolerance & newIntense < meanSeg + tolerance)) = 1;  
end
%Outline the boundary of the segmented region
se = strel('disk',1,0);
imErode = imerode(seed,se);
seed = seed - imErode;
end

