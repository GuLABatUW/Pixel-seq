clc; clear; close all;
delete(gcp('nocreate'));
warning off;
%###############################################################################
% Updated based on the DlightP2_1, Speed optimlization,
% Change in the template image generation, for base calling 
%###############################################################################
%parameter the sequencing and the folder store all the images for one FOV
totalCycle = 24;
imageSize = 1024;
Fov = '4-Pos_003_003';
Image_dir = '/run/user/1000/gvfs/smb-share:server=gu-lab-nas.local,share=image%20file/SPT_20210818_0816_S6_36cyc/';
nthread = 8;

%Cycles used to generate template, the quality could be manually checked.
%channel-crosstalk correction for template chasity value calculation was
%performed using default matrix./1:8;4;0.55;3;25;1.5;3.0
templateCycles = 1:8;
maxNoPurityBase = 4;
chasityCutoff = 0.52;

%image preprocessing Filter could be 'Raw', 'Laplacian', 'Sobel'
Filter = 'Raw';
transMo = [1 0 0;0 1 0;0 0 1];

%The spot finding algorithm can be setted as 'globalThesh' or 'localMaxima'
%or 'localThesh' or 'localBkg' or 'ltlaplacian' or 'ladaptG'
subregionRegNum = 1;
spotAlgorithm = 'localThesh';

%image registration control, 1=eable or 0=disable
imagereg = 1;

%radius threshold to merge the spots
%these thresholds are closed correlated with cluster size and cluster
%density. Generally,these parameters are setted according to the following
%rules:
%rs_minimum should be the median value for the distribution from the background noise to the maximum single spot.
%rs_maximum should be the median value for single spots extending to fusion spots.
%r1 should be the mean value for the same-base distance distribution
%r2 should be the mean value for the diff-base distance distribution
rs_minimum = 2;
rs_maximum = 25;
r1 = 1.5;
r2 = 3.0;
pr = 4.0; %neighbor radius
Rc = 0.65; %neighbor correlation cutoff
Neighbor = 'F'; %calculate neighbor or not with 'T' or 'F'
% image pre-processing to pre-sharp all the images and save them 
SC3 = imagePreprocessing(Image_dir,Fov,totalCycle,imageSize,'Cy3',Filter);
SC5 = imagePreprocessing(Image_dir,Fov,totalCycle,imageSize,'Cy5',Filter);
SNI = imagePreprocessing(Image_dir,Fov,totalCycle,imageSize,'NIR',Filter);
STX = imagePreprocessing(Image_dir,Fov,totalCycle,imageSize,'TxR',Filter);

%generate template image for the selected cycle
Sfirst_template = templateGeneration(SC3,SC5,SNI,STX,templateCycles(2),imageSize,0.65);
%tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%image registration
tmatrix = zeros(totalCycle,3,3);
%initiate the registered image
rSC5 = zeros(totalCycle,imageSize,imageSize);
rSNI = zeros(totalCycle,imageSize,imageSize);
rSC3 = zeros(totalCycle,imageSize,imageSize);
rSTX = zeros(totalCycle,imageSize,imageSize);
for i=1:totalCycle
    TSC5 = reshape(SC5(i,:,:),imageSize,imageSize);
    TSC3 = reshape(SC3(i,:,:),imageSize,imageSize);
    TSNI = reshape(SNI(i,:,:),imageSize,imageSize);
    TSTX = reshape(STX(i,:,:),imageSize,imageSize);
    regChannel = 0;
	if(imagereg == 1)
        sharedtemplate = Sfirst_template;
        regSC3 = TSC3;
        regSC5 = TSC5;
        regSNI = TSNI;
        regSTX = TSTX;
    
        %using channel to registration and same transformation matrix for all
        lastwarn('');
        [optimizer, metric] = imregconfig('multimodal');
        %transM = imregtform(regSC3,sharedtemplate,'translation', optimizer, metric);
        transM = imregcorr(regSC3,sharedtemplate,'translation');
        if(lastwarn())
        else %if(abs(transM.T(1,1)-1)<0.005 && abs(transM.T(2,2)-1)<0.005 && abs(transM.T(3,3)-1)<0.005)
            regChannel = 1;
        end

        if(regChannel == 0)
            lastwarn('');
            %transM = imregtform(regSTX,sharedtemplate,'translation', optimizer, metric);
            transM = imregcorr(regSC5,sharedtemplate,'translation');
            if(lastwarn())
            else
                regChannel = 1;  
            end
        end
        
        if(regChannel == 0)
            lastwarn('');
            %transM = imregtform(regSC5,sharedtemplate,'translation', optimizer, metric);
            transM = imregcorr(regSNI,sharedtemplate,'translation');
            if(lastwarn())
            else
                regChannel = 1;  
            end
        end
    
        if(regChannel == 0)
            lastwarn('');
            %transM = imregtform(regSNI,sharedtemplate,'translation', optimizer, metric);
            transM = imregcorr(regSTX,sharedtemplate,'translation');
            if(lastwarn())
                if(i==1)
                    transM.T = transMo;
                else
                    %if failed assigned to the previous transform matrix
                    transM.T = reshape(tmatrix(i-1,:,:),3,3);
                end
            end
        end
    
        tmatrix(i,:,:) = transM.T;
        SPC3 = imwarp(TSC3,transM,'OutputView',imref2d(size(sharedtemplate)));
        SPC5 = imwarp(TSC5,transM,'OutputView',imref2d(size(sharedtemplate)));
        SPNI = imwarp(TSNI,transM,'OutputView',imref2d(size(sharedtemplate)));
        SPTX = imwarp(TSTX,transM,'OutputView',imref2d(size(sharedtemplate)));
  else
      SPC3 = TSC3;
      SPC5 = TSC5;
      SPNI = TSNI;
      SPTX = TSTX;
  end
  
  rSC5(i,:,:) = SPC5;
  rSNI(i,:,:) = SPNI;
  rSC3(i,:,:) = SPC3;
  rSTX(i,:,:) = SPTX;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spoting finding for each image
%the default maximum spots for each cycle 1500k/mm2
spotAll = zeros(size(templateCycles,2)*1500000,5);
totSpots = 0;
CurrentStoreSpots = 0;
for i = templateCycles
  
  %generate transform matrix
  TSC5 = reshape(rSC5(i,:,:),imageSize,imageSize);
  TSC3 = reshape(rSC3(i,:,:),imageSize,imageSize);
  TSNI = reshape(rSNI(i,:,:),imageSize,imageSize);
  TSTX = reshape(rSTX(i,:,:),imageSize,imageSize);
  
  %store the cycle and channel information
  spotCurrentC3 = spot_finding(TSC3,i,3,spotAlgorithm);
  spotCurrentC5 = spot_finding(TSC5,i,5,spotAlgorithm);
  spotCurrentNI = spot_finding(TSNI,i,7,spotAlgorithm);
  spotCurrentTX = spot_finding(TSTX,i,9,spotAlgorithm);
  CurrentStoreSpots = size(spotCurrentC3,1)+size(spotCurrentC5,1)+size(spotCurrentNI,1)+size(spotCurrentTX,1);
  
  totSpots=totSpots+CurrentStoreSpots;
  spotAll((totSpots+1):(totSpots+CurrentStoreSpots),:) = [spotCurrentC3;spotCurrentC5;spotCurrentNI;spotCurrentTX];
end

%filter the spots rules based on NON-REPEAT and ISOLATED CLUSTERs
filterspot = spotAll(spotAll(:,3)>=rs_minimum & spotAll(:,3)<=rs_maximum,:);


%filter the spots with more than two template cycles (chastity < 0.6)
alltempSpot = chasityCheckn(filterspot,rSC3,rSC5,rSNI,rSTX,imageSize,templateCycles,chasityCutoff,maxNoPurityBase,nthread,r1,r2,'profileA1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bilinear introplate intensity identification
%intensity extraction based spots template
%this step will also perform background subtraction

intensityAll = zeros(size(alltempSpot,1),totalCycle*4+3);
%define the region of each spots
intensityPixel = zeros(size(alltempSpot,1)*25,totalCycle*4+3);
for j = 1:size(alltempSpot,1)
      xcord = alltempSpot(j,4);
      ycord = alltempSpot(j,5);
      if(xcord>1023.5 || ycord>1023.5 || xcord<1.5 || ycord<1.5)
          continue
      end
      region = (fix(xcord/256.1)+1)*10+(fix(ycord/256.1)+1);
      intensityAll(j,1) = intensityAll(j,1)+region;
      intensityAll(j,2:3) = alltempSpot(j,4:5);
      aIN = (j-1)*25+1;
      bIN = j*25;
      xIN = fix(xcord)-2:fix(xcord)+2;
      yIN = fix(ycord)-2:fix(ycord)+2;
      intensityPixel(aIN:bIN,2:3) = combvec(xIN,yIN)';
      
end

%build the neighbor pixel matrix
intensityUPixel = unique(intensityPixel,'rows');
intensityUnPixel = intensityUPixel(intensityUPixel(:,2)>0 & intensityUPixel(:,3)>0 & intensityUPixel(:,2)<=1024 & intensityUPixel(:,3)<=1024,:);
%define single value index
sIN = (intensityUnPixel(:,2)-1)*1024+intensityUnPixel(:,3);


%extract the template cycle intensity
intensityCy3 = zeros(size(alltempSpot,1),totalCycle);
intensityCy5 = zeros(size(alltempSpot,1),totalCycle);
intensityNIR = zeros(size(alltempSpot,1),totalCycle);
intensityTXR = zeros(size(alltempSpot,1),totalCycle);

pixelCy3 = zeros(size(intensityUnPixel,1),totalCycle);
pixelCy5 = zeros(size(intensityUnPixel,1),totalCycle);
pixelNIR = zeros(size(intensityUnPixel,1),totalCycle);
pixelTXR = zeros(size(intensityUnPixel,1),totalCycle);


for i = templateCycles
    TSC5 = reshape(rSC5(i,:,:),imageSize,imageSize);
    TSC3 = reshape(rSC3(i,:,:),imageSize,imageSize);
    TSNI = reshape(rSNI(i,:,:),imageSize,imageSize);
    TSTX = reshape(rSTX(i,:,:),imageSize,imageSize);
    
    intensityCy5(:,i) = alltempSpot(:,4*i+3);
    intensityNIR(:,i) = alltempSpot(:,4*i+4);
    intensityCy3(:,i) = alltempSpot(:,4*i+5);
    intensityTXR(:,i) = alltempSpot(:,4*i+6);
    
    pixelCy3(:,i) = TSC3(sIN);
    pixelCy5(:,i) = TSC5(sIN);
    pixelNIR(:,i) = TSNI(sIN);
    pixelTXR(:,i) = TSTX(sIN);
end



parpool('profileA2',nthread);
parfor i=(max(templateCycles)+1):totalCycle
  
    TSC5 = reshape(rSC5(i,:,:),imageSize,imageSize);
    TSC3 = reshape(rSC3(i,:,:),imageSize,imageSize);
    TSNI = reshape(rSNI(i,:,:),imageSize,imageSize);
    TSTX = reshape(rSTX(i,:,:),imageSize,imageSize);
    
    temIntensityC3 = zeros(size(alltempSpot,1),1);
    temIntensityC5 = zeros(size(alltempSpot,1),1);
    temIntensityNI = zeros(size(alltempSpot,1),1);
    temIntensityTX = zeros(size(alltempSpot,1),1);
    for j = 1:size(alltempSpot,1)
        xd = alltempSpot(j,4)+0.5;
        yd = alltempSpot(j,5)+0.5;
        temIntensityC3(j,1) = biolinear(TSC3,16,xd,yd);
        temIntensityC5(j,1) = biolinear(TSC5,16,xd,yd);
        temIntensityNI(j,1) = biolinear(TSNI,16,xd,yd);
        temIntensityTX(j,1) = biolinear(TSTX,16,xd,yd);
    end
    intensityCy3(:,i) = temIntensityC3(:,1);
    intensityCy5(:,i) = temIntensityC5(:,1);
    intensityNIR(:,i) = temIntensityNI(:,1);
    intensityTXR(:,i) = temIntensityTX(:,1);
    
    pixelCy3(:,i) = TSC3(sIN);
    pixelCy5(:,i) = TSC5(sIN);
    pixelNIR(:,i) = TSNI(sIN);
    pixelTXR(:,i) = TSTX(sIN);
    
end

intensityAll(:,(1:totalCycle)*4) = intensityCy5;
intensityAll(:,(1:totalCycle)*4+1) = intensityNIR;
intensityAll(:,(1:totalCycle)*4+2) = intensityCy3;
intensityAll(:,(1:totalCycle)*4+3) = intensityTXR;

intensityUnPixel(:,(1:totalCycle)*4) = pixelCy5;
intensityUnPixel(:,(1:totalCycle)*4+1) = pixelNIR;
intensityUnPixel(:,(1:totalCycle)*4+2) = pixelCy3;
intensityUnPixel(:,(1:totalCycle)*4+3) = pixelTXR;

delete(gcp('nocreate'));


%%%%%%%%%%%%%%%%%%%%%%%combine four channel intensity into one



%filter the region 00 induced by the program
intensityFormal= intensityAll(intensityAll(:,1)>0,:);

intensityFormalR = intensityFormal(intensityFormal(:,8)<quantile(intensityFormal(:,8),0.999) ...
    & intensityFormal(:,9)<quantile(intensityFormal(:,9),0.999) ...
    & intensityFormal(:,10)<quantile(intensityFormal(:,10),0.999) ...
    & intensityFormal(:,11)<quantile(intensityFormal(:,11),0.999),:);
%intensityFormalR(:,46) = intensityFormalR(:,46)*0.8;
intensityFormalRP = intensityFormalR(randperm(size(intensityFormalR,1)),:);
csvwrite('/home/xiaonan/lab508/mockIllumina/basecallingTest/intensityFormal.csv',intensityFormalR);

%create neighbor projection by correlation
if(strcmp(Neighbor, 'T'))
    idx = rangesearch(intensityFormalR(:,2:3),intensityUnPixel(:,2:3),pr);
    pixelProjection = zeros(size(intensityUnPixel,1),5);
    pixelProjection(:,1:2) = intensityUnPixel(:,2:3);
    for i = 1:size(idx,1)
        b = idx{i};
        if(~isempty(b))
            tempCorr = 0;
            for j = 1:size(b,2)
                corrP_S = corr(intensityUnPixel(i,4:(totalCycle*4+3))', intensityFormalR(b(j),4:(totalCycle*4+3))');
                if(corrP_S >= Rc && corrP_S>tempCorr)
                    pixelProjection(i,3:4) = intensityFormalR(b(j),2:3);
                    pixelProjection(i,5) = b(j);
                    tempCorr = corrP_S;
                end
            end
        end
    end
    pixelFProjection = pixelProjection(pixelProjection(:,3)>0 & pixelProjection(:,4)>0,:);
    csvwrite('/home/xiaonan/lab508/mockIllumina/basecallingTest/intensityNeighbor.csv',pixelFProjection);
end

if(strcmp(Neighbor, 'F'))
    csvwrite('/home/xiaonan/lab508/mockIllumina/basecallingTest/intensityNeighbor.csv',intensityFormalR(:,2:3));
end
%timeElapsed = toc
%tmp = intensityFormalR(intensityFormalR(:,2)>200 & intensityFormalR(:,3)>200,:);
%intensityTmp = tmp(randperm(size(tmp,1)),:);

%csvwrite('intensityFormal.csv',intensityTmp);

