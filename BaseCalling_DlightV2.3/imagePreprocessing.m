


%image pre-processing to pre-sharp all the images 

function [images] = imagePreprocessing(inputDir,Fov, totalCyc,size,channel,filter)
    images = zeros(totalCyc,size,size);
    for i = 1:totalCyc
        j = i;
        %j = i+1;
        img = [inputDir,'Incorp_',num2str(j),'/',Fov,'/img_000000000_',channel,'_000.tif'];
        %images(i,:,:) = imadjust(im2double(imread(img)));
        if(strcmp(filter, 'Raw'))
            images(i,:,:) = im2double(imread(img));
        end
        
        
        %as illumina indicated in RTA, perform Laplacian presharpening
        if(strcmp(filter, 'Laplacian'))
            lap = [-1 -1 -1; -1 32 -1; -1 -1 -1]; %laplacian filter
            image1 = im2double(imread(img));
            imageLap = imfilter(image1, lap, 'conv'); 

            %// Change - Adding to original image now
            imageSharpened = image1 + imageLap;
            images(i,:,:) = imageSharpened;
        end
        
        %DOI: 10.21917/ijivp.2016.0180
        %combine Laplacian filter with sobel gradient
        if(strcmp(filter, 'Sobel'))
            lap = [-1 -1 -1; -1 8 -1; -1 -1 -1]; %laplacian filter
            image1 = im2double(imread(img));
            imageLap = imfilter(image1, lap, 'conv'); 
            imageSharpened = image1 + imageLap;
            
            imageSobel = imgradient(image1);
            imageAvgSobel = imfilter(imageSobel, ones(5)/25, 'symmetric');
            imageMask = imageAvgSobel.*imageSharpened;
            imageFinalS = image1 + imageMask;
            
            images(i,:,:) = imageFinalS;
        end
        
    end

