
%method the find the maximum true spots 

function [spots] = spot_finding(image,cycle,chanel,algorithm)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % local maxima search of connectivity 4
  % Arthur 7/17/2019
  
  % image registration
  I_ori = image;
  
if(strcmp(algorithm, 'localMaxima'))
      % local max search won't give any area
      % this hack defines all area as 4
      % so all clusters would survive later size filtering
      hacking_area = 6;

      % get the local max map
      lm_map = imregionalmax(I_ori, 4);
      lm_map = im2double(lm_map);

      % find threshold = 1 std above mean
%       median_itst = median(I_ori);
%       sigma = std2(I_ori);
%       itst_thrsh = median_itst + sigma;

      % apply local max and threshold on original image
      local_maxima = lm_map .* I_ori;
      I_new = local_maxima > 0;

      % collect all points
      spot_centers = regionprops(I_new, 'Centroid');
      num_clusters = numel(spot_centers);
      reformated_spot_centers = reshape(struct2array(spot_centers), [], numel(spot_centers))';
      spots = horzcat(repelem(cycle,  num_clusters)', repelem(chanel, num_clusters)', repelem(hacking_area, num_clusters)', reformated_spot_centers);
end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(algorithm, 'globalThesh'))
   I = I_ori;
   Iblur = imgaussfilt(I);
   I_new = I - Iblur;
   
   
   thresh = median(I_new(:)) + 2*std2(I_new);
   I_new = I_new > thresh;
   I_new = bwareaopen(I_new,2);
 
   ttspots = regionprops(I_new,'Centroid','Area');
   
   % result is in the following format
   % area, position x, position y (position x and y are for the centroid"
   result = reshape(struct2array(ttspots), [], numel(ttspots)).';
   result(:,2) = result(:,2);
   result(:,3) = result(:,3);
   
   spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result];
end


if(strcmp(algorithm, 'localThesh'))
    %regionSize were defined as squres such as 64*64
   regionSize = 32;
   I = I_ori;
   Iblur = imgaussfilt(I);
   I_new = I - Iblur;
   for i = 1:(1024/regionSize)
       for j = 1:(1024/regionSize)
           stRow = (i-1)*regionSize+1;
           edRow = i*regionSize;
           stCol = (j-1)*regionSize+1;
           edCol = j*regionSize;
           Inew_region = I_new(stRow:edRow,stCol:edCol);
           thresh = median(Inew_region(:)) + 2*std2(Inew_region);
           I_new(stRow:edRow,stCol:edCol) = I_new(stRow:edRow,stCol:edCol) > thresh;
           
       end
   end
   I_new = bwareaopen(I_new,2);
   ttspots = regionprops(I_new,'Centroid','Area');
   
   % result is in the following format
   % area, position x, position y (position x and y are for the centroid"
   result = reshape(struct2array(ttspots), [], numel(ttspots)).';
   result(:,2) = result(:,2);
   result(:,3) = result(:,3);
   
   spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result];
end


if(strcmp(algorithm, 'localBkg'))
    %regionSize were defined as squres such as 64*64
   regionSize = 64;
   breaks = 100;
   bkgSNR = 2.8284;
   I = I_ori;
   I_new = I;
   for i = 1:(1024/regionSize)
       for j = 1:(1024/regionSize)
           stRow = (i-1)*regionSize+1;
           edRow = i*regionSize;
           stCol = (j-1)*regionSize+1;
           edCol = j*regionSize;
           Inew_region = I_new(stRow:edRow,stCol:edCol);
           intenAll = Inew_region(:);
  
           [weightedbin,edges] = histcounts(intenAll,breaks);
           [pks,locs] = findpeaks(weightedbin,'MinPeakProminence',5);
           select = [pks;locs]';
           thresh = bkgSNR*edges(select(1,2));
           I_new(stRow:edRow,stCol:edCol) = I_new(stRow:edRow,stCol:edCol) > thresh;
           
       end
   end
   I_new = bwareaopen(I_new,2);
   ttspots = regionprops(I_new,'Centroid','Area');
   
   % result is in the following format
   % area, position x, position y (position x and y are for the centroid"
   result = reshape(struct2array(ttspots), [], numel(ttspots)).';
   result(:,2) = result(:,2);
   result(:,3) = result(:,3);
   
   spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result];
end

if(strcmp(algorithm, 'ltlaplacian'))
    %regionSize were defined as squres such as 64*64
    
   regionSize = 64;
   lap = [-1 -1 -1; -1 8 -1; -1 -1 -1]; %laplacian filter
   image1 = I_ori;
   imageLap = imfilter(image1, lap, 'conv'); 
   I = image1 + imageLap;

   I_new = I;
   for i = 1:(1024/regionSize)
       for j = 1:(1024/regionSize)
           stRow = (i-1)*regionSize+1;
           edRow = i*regionSize;
           stCol = (j-1)*regionSize+1;
           edCol = j*regionSize;
           Inew_region = I_new(stRow:edRow,stCol:edCol);
           thresh = median(Inew_region(:)) +2*std2(Inew_region);
           I_new(stRow:edRow,stCol:edCol) = I_new(stRow:edRow,stCol:edCol) > thresh;
           
       end
   end
   I_new = bwareaopen(I_new,2);
   ttspots = regionprops(I_new,'Centroid','Area');
   
   % result is in the following format
   % area, position x, position y (position x and y are for the centroid"
   result = reshape(struct2array(ttspots), [], numel(ttspots)).';
   result(:,2) = result(:,2);
   result(:,3) = result(:,3);
   
   spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result];
end


if(strcmp(algorithm, 'ladaptG'))
   %using the matlab function adaptThresh method
   regionSize = 64;
   sensitivity = 0.18;
   I = I_ori;
   %Iblur = imgaussfilt(I);
   I_new = I;
   for i = 1:(1024/regionSize)
       for j = 1:(1024/regionSize)
           stRow = (i-1)*regionSize+1;
           edRow = i*regionSize;
           stCol = (j-1)*regionSize+1;
           edCol = j*regionSize;
           Inew_region = I_new(stRow:edRow,stCol:edCol);
           T = adaptthresh(Inew_region,sensitivity,'Statistic','gaussian');
           I_new(stRow:edRow,stCol:edCol) = imbinarize(I_new(stRow:edRow,stCol:edCol),T);
       end
   end
   I_new = bwareaopen(I_new,2);
   ttspots = regionprops(I_new,'Centroid','Area');
   
   % result is in the following format
   % area, position x, position y (position x and y are for the centroid"
   result = reshape(struct2array(ttspots), [], numel(ttspots)).';
   result(:,2) = result(:,2);
   result(:,3) = result(:,3);
   
   spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result];
end

if(strcmp(algorithm, 'localTheshS'))
    %regionSize were defined as squres such as 64*64
   regionSize = 359;
   I = I_ori;
   Iblur = imgaussfilt(I);
   I_new = I - Iblur;
   for i = 1:(718/regionSize)
       for j = 1:(718/regionSize)
           stRow = (i-1)*regionSize+1;
           edRow = i*regionSize;
           stCol = (j-1)*regionSize+1;
           edCol = j*regionSize;
           Inew_region = I_new(stRow:edRow,stCol:edCol);
           thresh = median(Inew_region(:)) + 2*std2(Inew_region);
           I_new(stRow:edRow,stCol:edCol) = I_new(stRow:edRow,stCol:edCol) > thresh;
           
       end
   end
   I_new = bwareaopen(I_new,2);
   ttspots = regionprops(I_new,'Centroid','Area');
   
   % result is in the following format
   % area, position x, position y (position x and y are for the centroid"
   result = reshape(struct2array(ttspots), [], numel(ttspots)).';
   result(:,2) = result(:,2);
   result(:,3) = result(:,3);
   
   spots = [repelem(cycle,size(result,1))' repelem(chanel,size(result,1))' result];
end
