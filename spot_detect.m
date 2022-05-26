function [spots] = spot_detect(image,cycle,chanel,filter)

       I_ori = image;
        %regionSize were defined as squres such as 64*64
       I = I_ori;
       Iblur = imgaussfilt(I);
       I_new = I - Iblur;

       thresh = median(I_new(:)) + std2(I_new);
       I_new = I_new > thresh;
       I_new = bwareaopen(I_new,filter);

       ttspots = regionprops(I_new,'Centroid','Area');

       % result is in the following format
       % area, position x, position y (position x and y are for the centroid"
       result = reshape(struct2array(ttspots), [], numel(ttspots)).';
       result(:,2) = result(:,2);
       result(:,3) = result(:,3);
       
       xlimUP = size(image,1)-0.5;
       ylimUP = size(image,2)-0.5;
       
       result = result(result(:,2)>=1.5 & result(:,3)>=1.5 & result(:,2)<=ylimUP & result(:,3)<=xlimUP,:);
      
       Intensity = zeros(size(result,1),1);
       for j = 1:size(result,1)
            xd = result(j,2);
            yd = result(j,3);
            Intensity(j,1) = biolinear(image,16,xd,yd);
       end
       
       spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result Intensity];
