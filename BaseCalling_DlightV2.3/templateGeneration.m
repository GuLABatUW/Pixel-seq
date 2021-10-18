            


function [tempImage] = templateGeneration(SC3,SC5,SNI,STX,template,imageSize,Cover)

NC3 = reshape(SC3(template,:,:),[imageSize imageSize]);
NC5 = reshape(SC5(template,:,:),[imageSize imageSize]);
NNI = reshape(SNI(template,:,:),[imageSize imageSize]);
NTX = reshape(STX(template,:,:),[imageSize imageSize]);
logFilter = -fspecial('log',5);
tNC3 = imfilter(NC3,logFilter);
tNC5 = imfilter(NC5,logFilter);
tNNI = imfilter(NNI,logFilter);
tNTX = imfilter(NTX,logFilter);

SC5N = imfuse(tNNI,tNC5,'blend','Scaling','joint');
SC3T = imfuse(tNC3,tNTX,'blend','Scaling','joint');

if(Cover>0.7)
    tempImage = SC3T;
else
    tempImage = imadjust(imfuse(SC5N,SC3T,'blend','Scaling','joint'));
end