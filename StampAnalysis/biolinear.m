%function to get the subpixel interpolation to accuratly extract the
%intensity information.

%variable [image]: the input image intensity matrix
%variable [surround_pix]; the radius of surrounding pixels used to estimate
%the background
%variable [x,y]: the coordinates of the centroid



function [intensity] = biolinear(image,surround_pix,x,y)
    x0 = x-0.5;
    x1 = fix(x+0.5);
    x2 = x+0.5;
    
    y0 = y-0.5;
    y1 = fix(y+0.5);
    y2 = y+0.5;
    
    area1 = (x1-x0)*(y2-y1)*image(fix(y2),fix(x0));
    area2 = (x2-x1)*(y2-y1)*image(fix(y2),fix(x2));
    area3 = (x2-x1)*(y1-y0)*image(fix(y0),fix(x2));
    area4 = (x1-x0)*(y1-y0)*image(fix(y0),fix(x0));
    
    %estimate the background by 32*32
    rangeMinX = x1 -surround_pix;
    rangeMaxX = x1 +surround_pix;
    rangeMinY = y1 -surround_pix;
    rangeMaxY = y1 +surround_pix;
    if(x1<(surround_pix+1))
        rangeMinX = 1;
    end
    
    if(x1>(size(image,2)-surround_pix))
        rangeMaxX = size(image,2);
    end
    
    if(y1<(surround_pix+1))
        rangeMinY = 1;
    end
    
    if(y1>(size(image,1)-surround_pix))
        rangeMaxY = size(image,1);
    end
    estimateArea = image(rangeMinY:rangeMaxY,rangeMinX:rangeMaxX);
    minlist = sort(estimateArea(:));
    min4 = mean(minlist(1:4));
    
    intensity = area1+area2+area3+area4-min4;
