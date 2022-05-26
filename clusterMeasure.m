clc; clear; close all;



%loading the image for analysis
stamp2 = imadjust(imread('/media/xiaonan/T7/slideStamp_final/processed/Copy 2_Z09_crop.tif'));
stamp10 = imadjust(imread('/media/xiaonan/T7/slideStamp_final/processed/Copy 10_Z08_crop.tif'));
stamp50 = imadjust(imread('/media/xiaonan/T7/slideStamp_final/processed/Copy 50_Z09_crop.tif'));
stamp2 = im2double(stamp2);
stamp10 = im2double(stamp10);
stamp50 = im2double(stamp50);



%register the image to the stamp50, 
%stamp50's quality is better than the other two 
transM2 = imregcorr(stamp2,stamp50,'similarity');
trans_stamp2 = imwarp(stamp2,transM,'OutputView',imref2d(size(stamp50)));
imshowpair(stamp50, trans_stamp2);

transM10 = imregcorr(stamp10,stamp50,'similarity');
trans_stamp10 = imwarp(stamp10,transM10,'OutputView',imref2d(size(stamp50)));
imshowpair(stamp50, trans_stamp10);


%detect the spots and keep the center and the intensity
spotStamp2 = spot_detect(trans_stamp2, 1, 3, 3);
spotStamp10 = spot_detect(trans_stamp10, 1, 3, 3);
spotStamp50 = spot_detect(stamp50, 1,3,3);
imshow(stamp50);
hold;
scatter(spotStamp50(:,4),spotStamp50(:,5));

csvwrite('/media/xiaonan/T7/slideStamp_final/processed/stamp2_cluster.csv',spotStamp2);
csvwrite('/media/xiaonan/T7/slideStamp_final/processed/stamp10_cluster.csv',spotStamp10);
csvwrite('/media/xiaonan/T7/slideStamp_final/processed/stamp50_cluster.csv',spotStamp50);
