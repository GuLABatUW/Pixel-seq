#whole image analysis
probeCon = im2double(imread('./beforeCleavage.tif'));

probeTaqI =  im2double(imread('./afterCleavage.tif'));

csvwrite('probeTaqI.csv',probeTaqI(:));
csvwrite('probeCon.csv',probeCon(:));
