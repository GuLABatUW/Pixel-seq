function [alltempSpot] = chasityCheckn(definedspots,checkImgSC3,checkImgSC5,checkImgSNI,checkImgSTX,imgSize,selectedCycles,Pcutoff,maxNoPurity,nTh,r1,r2,pfname)


%considering the raw spots are not quality-reliable, color crosstalk
%correction will be performed by default matrix

%crosstalkM = [1 1.8479 0.5415 0.4633;0.1695 1 0.1098 0.1590;0.0295 0.0739 1 0.5711;0.0301 0.1416 0.3666 1];
templateIntensity = zeros(size(definedspots,1),size(selectedCycles,2)*4+3);
storeSpots = zeros(size(definedspots,1),size(definedspots,2)+1);
storeSpots(:,1:5) = definedspots;
%define the region of each spots
for j = 1:size(definedspots,1)
      xcord = definedspots(j,4);
      ycord = definedspots(j,5);
      if(xcord>1023.5 || ycord>1023.5 || xcord<1.5 || ycord<1.5)
          continue
      end
      region = (fix(xcord/256.1)+1)*10+(fix(ycord/256.1)+1);
      templateIntensity(j,1) = templateIntensity(j,1)+region;
      templateIntensity(j,2:3) = definedspots(j,4:5);
      storeSpots(j,6) = storeSpots(j,6) + region;
end

%filter the 00 region induced by the program

templateClean= templateIntensity(templateIntensity(:,1)>0,:);
spotsClean = storeSpots(storeSpots(:,6)>0,1:5);


mycluster = parpool(pfname,nTh);
intensityCy3 = zeros(size(templateClean,1),size(selectedCycles,2));
intensityCy5 = zeros(size(templateClean,1),size(selectedCycles,2));
intensityNIR = zeros(size(templateClean,1),size(selectedCycles,2));
intensityTXR = zeros(size(templateClean,1),size(selectedCycles,2));

parfor i=selectedCycles
    TSC3 = reshape(checkImgSC3(i,:,:),imgSize,imgSize);
    TSC5 = reshape(checkImgSC5(i,:,:),imgSize,imgSize);
    TSNI = reshape(checkImgSNI(i,:,:),imgSize,imgSize);
    TSTX = reshape(checkImgSTX(i,:,:),imgSize,imgSize);
  
    temIntensityC3 = zeros(size(templateClean,1),1);
    temIntensityC5 = zeros(size(templateClean,1),1);
    temIntensityNI = zeros(size(templateClean,1),1);
    temIntensityTX = zeros(size(templateClean,1),1);
    for j = 1:size(templateClean,1)
        xd = templateClean(j,2)+0.5;
        yd = templateClean(j,3)+0.5;
        temIntensityC3(j,1) = biolinear(TSC3,16,xd,yd);
        temIntensityC5(j,1) = biolinear(TSC5,16,xd,yd);
        temIntensityNI(j,1) = biolinear(TSNI,16,xd,yd);
        temIntensityTX(j,1) = biolinear(TSTX,16,xd,yd);
    end
    intensityCy3(:,i) = temIntensityC3(:,1);
    intensityCy5(:,i) = temIntensityC5(:,1);
    intensityNIR(:,i) = temIntensityNI(:,1);
    intensityTXR(:,i) = temIntensityTX(:,1);
end
delete(gcp('nocreate'));
  
templateClean(:,(1:size(selectedCycles,2))*4) = intensityCy5(:,selectedCycles);
templateClean(:,(1:size(selectedCycles,2))*4+1) = intensityNIR(:,selectedCycles);
templateClean(:,(1:size(selectedCycles,2))*4+2) = intensityCy3(:,selectedCycles);
templateClean(:,(1:size(selectedCycles,2))*4+3) = intensityTXR(:,selectedCycles);

rawTemplateIntensity =  templateClean; 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%spatial normalization
  for i = 4:size(templateClean,2)
    %step 1 find the brightest region 
    meanIntenByRegion = zeros(4,4);
    for j = 1:4
        for k = 1:4
            meanIntenByRegion(j,k) = mean(templateClean(templateClean(:,1)==((j*10)+k),i));
        end
    end
    [row,col] = find(ismember(meanIntenByRegion,max(meanIntenByRegion(:))));
    brightInd = row*10+col;
    brightRegion = templateClean(templateClean(:,1)==brightInd,i);
    
    %prepare to have spatial normalization
    for j = 1:4
        for k = 1:4
            currentInd = j*10+k;
            if(currentInd==brightInd)
                continue;
            end
            
            currentRegion = templateClean(templateClean(:,1)==currentInd,i);

            normRegion = quantile(brightRegion(:),tiedrank(currentRegion)/length(currentRegion));
      
            %put it back
            templateClean(templateClean(:,1)==currentInd,i) = normRegion;
        end
    end
  end
      
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%crosstalk correction
  for i = 1:size(selectedCycles,2)
    st = (i-1)*4+4;
    ed = (i-1)*4+7;
    
    ttC3 = checkImgSC3(i,412:612,412:612);
    ttC5 = checkImgSC5(i,412:612,412:612);
    ttTX = checkImgSTX(i,412:612,412:612);
    ttNI = checkImgSNI(i,412:612,412:612);
    
    templateIntensity = [ttC5(:) ttNI(:) ttC3(:) ttTX(:)];
    crosstalkM = zeros(4,4);
    for j = 1:4
        for k = 1:4
            if j==k
                crosstalkM(j,k) = 1;
            else
                crosstalkM(j,k) = weightedHist(templateIntensity(:,j),templateIntensity(:,k),100,1);
            end
        end
    end
    templateClean(:,st:ed) = templateClean(:,st:ed)/crosstalkM;
  end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%channel correction
    for i = 4:size(templateClean,2)
        current = templateClean(:,i);
        percentile05 = quantile(current,0.05);
        percentile95 = quantile(current,0.95);
        current(current<percentile05) = 0;
        current(current>percentile95) = 100;
        current(current>=percentile05 & current<=percentile95) = (current(current>=percentile05 & current<=percentile95)-percentile05)/(percentile95-percentile05)*98+1;
        templateClean(:,i) = current;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Chasity value calculation
    chasityNumOfFirstFour = zeros(size(templateClean,1),1);
    %chasityAll = zeros(size(templateClean,1),size(selectedCycles,2));
    for j = 1:size(templateClean,1)
        chty = selectedCycles;
        for i = 1:size(selectedCycles,2)
            st = (i-1)*4+4;
            ed = (i-1)*4+7;
            sortI = sort(templateClean(j,st:ed));
            chty(i) = sortI(4)/(sortI(3)+sortI(4));
            %chasityAll(j,i) = chty(i);
        end
        chasityNumOfFirstFour(j) = size(chty(chty<Pcutoff),2);
    end
    
    cleansp = [spotsClean chasityNumOfFirstFour rawTemplateIntensity(:,4:35)];
    clean = cleansp(cleansp(:,6)<maxNoPurity,:);
    cleanSpot = clean;
  
    %count the spots for each cycles and each channel
    spotStat = zeros(size(selectedCycles,2),6);
    for i = 1:size(selectedCycles,2)
        spotStat(i,1) = selectedCycles(1,i);
        spotStat(i,2) = size(cleanSpot(cleanSpot(:,1)==selectedCycles(1,i) & cleanSpot(:,2)==3,:),1);
        spotStat(i,3) = size(cleanSpot(cleanSpot(:,1)==selectedCycles(1,i) & cleanSpot(:,2)==5,:),1);
        spotStat(i,4) = size(cleanSpot(cleanSpot(:,1)==selectedCycles(1,i) & cleanSpot(:,2)==7,:),1);
        spotStat(i,5) = size(cleanSpot(cleanSpot(:,1)==selectedCycles(1,i) & cleanSpot(:,2)==9,:),1);
        spotStat(i,6) = spotStat(i,3) + spotStat(i,4);
    end

    %sort and find the gold and siliver cycle
    sortSpot = sortrows(spotStat,6,'descend');

    %merge the spots by the setting rules
    registerSpot = zeros(1000000,38);
    %step by step merge
    %step 1: merge C(NIR - 7) with A(cy5 - 5) from golden cycle

    goldenCycleCy3 = cleanSpot(cleanSpot(:,1)==sortSpot(1,1) & cleanSpot(:,2)==3,:);
    goldenCycleCy5 = cleanSpot(cleanSpot(:,1)==sortSpot(1,1) & cleanSpot(:,2)==5,:);
    goldenCycleNIR = cleanSpot(cleanSpot(:,1)==sortSpot(1,1) & cleanSpot(:,2)==7,:);
    goldenCycleTXR = cleanSpot(cleanSpot(:,1)==sortSpot(1,1) & cleanSpot(:,2)==9,:);

    %siliverCycleCy3 = cleanSpot(cleanSpot(:,1)==sortSpot(2,1) & cleanSpot(:,2)==3,:);
    siliverCycleCy5 = cleanSpot(cleanSpot(:,1)==sortSpot(2,1) & cleanSpot(:,2)==5,:);
    siliverCycleNIR = cleanSpot(cleanSpot(:,1)==sortSpot(2,1) & cleanSpot(:,2)==7,:);
    %siliverCycleTXR = cleanSpot(cleanSpot(:,1)==sortSpot(2,1) & cleanSpot(:,2)==9,:);

    %registered all the spots
    totRegisterSpot = size(goldenCycleCy5,1);
    registerSpot(1:totRegisterSpot,:) = goldenCycleCy5;

    %register NIR against Cy5
    idx = rangesearch(goldenCycleNIR(:,[4,5]),registerSpot(:,[4,5]),r1,'Distance','euclidean');
    points_within_range = unique([idx{:}]);
    nonOverlapSpots_NIR = goldenCycleNIR(setdiff(1:size(goldenCycleNIR,1),points_within_range), :);
    cSpot = size(nonOverlapSpots_NIR,1);
    totRegisterSpot = totRegisterSpot+cSpot;
    registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots_NIR;

    %prepare siliver cycles
    idx = rangesearch(siliverCycleNIR(:,[4,5]),siliverCycleCy5(:,[4,5]),r1);
    points_within_range = unique([idx{:}]);
    nonOverlapSpots_NIR = siliverCycleNIR(setdiff(1:size(siliverCycleNIR,1),points_within_range), :);
    siliverSpots = [siliverCycleCy5;nonOverlapSpots_NIR];

    %register siliver cycles against golden cycles
    regist5 = registerSpot(registerSpot(:,2)==5,:);
    regist7 = registerSpot(registerSpot(:,2)==7,:);
    siliver5 = siliverSpots(siliverSpots(:,2)==5,:);


    %ready for A C
    %A vs A + C
    idx = rangesearch(siliver5(:,[4,5]),regist5(:,[4,5]),r2);
    points_within_range = unique([idx{:}]);
    tmpSpot = siliver5(setdiff(1:size(siliver5,1),points_within_range), :);

    idx = rangesearch(tmpSpot(:,[4,5]),regist7(:,[4,5]),r1);
    points_within_range = unique([idx{:}]);
    nonOverlapSpots = tmpSpot(setdiff(1:size(tmpSpot,1),points_within_range), :);
    cSpot = size(nonOverlapSpots,1);
    totRegisterSpot = totRegisterSpot+cSpot;
    registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;


    %C vs A + C
    regist5 = registerSpot(registerSpot(:,2)==5,:);
    regist7 = registerSpot(registerSpot(:,2)==7,:);
    siliver7 = siliverSpots(siliverSpots(:,2)==7,:);

    idx = rangesearch(siliver7(:,[4,5]),regist7(:,[4,5]),r2);
    points_within_range = unique([idx{:}]);
    tmpSpot = siliver7(setdiff(1:size(siliver7,1),points_within_range), :);

    idx = rangesearch(tmpSpot(:,[4,5]),regist5(:,[4,5]),r1);
    points_within_range = unique([idx{:}]);
    nonOverlapSpots = tmpSpot(setdiff(1:size(tmpSpot,1),points_within_range), :);
    cSpot = size(nonOverlapSpots,1);
    totRegisterSpot = totRegisterSpot+cSpot;
    registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;

    %ready for golden T(TxR - 9) against A + C
    idx = rangesearch(goldenCycleTXR(:,[4,5]),registerSpot(:,[4,5]),r1);
    points_within_range = unique([idx{:}]);
    nonOverlapSpots = goldenCycleTXR(setdiff(1:size(goldenCycleTXR,1),points_within_range), :);
    cSpot = size(nonOverlapSpots,1);
    totRegisterSpot = totRegisterSpot+cSpot;
    registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;

    %ready for golden G(Cy3 - 3) against A+C
    idx = rangesearch(goldenCycleCy3(:,[4,5]),registerSpot(:,[4,5]),r1);
    points_within_range = unique([idx{:}]);
    nonOverlapSpots = goldenCycleCy3(setdiff(1:size(goldenCycleCy3,1),points_within_range), :);
    cSpot = size(nonOverlapSpots,1);
    totRegisterSpot = totRegisterSpot+cSpot;
    registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;

    %read to consider all the other template cycels with setting rules

    for i = 3:size(sortSpot,1)
        currentCycle = sortSpot(i,1);

        cy3_current_spot = cleanSpot(cleanSpot(:,1)==currentCycle & cleanSpot(:,2)==3,:);
        regist3 = registerSpot(registerSpot(:,2)==3,:);
        regist579 = registerSpot(registerSpot(:,2)~=3 & registerSpot(:,2)>0,:);
        idx = rangesearch(cy3_current_spot(:,[4,5]),regist3(:,[4,5]),r2);
        points_within_range = unique([idx{:}]);
        tmpSpot = cy3_current_spot(setdiff(1:size(cy3_current_spot,1),points_within_range), :);

        idx = rangesearch(tmpSpot(:,[4,5]),regist579(:,[4,5]),r1);
        points_within_range = unique([idx{:}]);
        nonOverlapSpots = tmpSpot(setdiff(1:size(tmpSpot,1),points_within_range), :);
        cSpot = size(nonOverlapSpots,1);
        totRegisterSpot = totRegisterSpot+cSpot;
        registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;



        cy5_current_spot = cleanSpot(cleanSpot(:,1)==currentCycle & cleanSpot(:,2)==5,:);
        regist5 = registerSpot(registerSpot(:,2)==5,:);
        regist379 = registerSpot(registerSpot(:,2)~=5 & registerSpot(:,2)>0,:);
        idx = rangesearch(cy5_current_spot(:,[4,5]),regist5(:,[4,5]),r2);
        points_within_range = unique([idx{:}]);
        tmpSpot = cy5_current_spot(setdiff(1:size(cy5_current_spot,1),points_within_range), :);

        idx = rangesearch(tmpSpot(:,[4,5]),regist379(:,[4,5]),r1);
        points_within_range = unique([idx{:}]);
        nonOverlapSpots = tmpSpot(setdiff(1:size(tmpSpot,1),points_within_range), :);
        cSpot = size(nonOverlapSpots,1);
        totRegisterSpot = totRegisterSpot+cSpot;
        registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;


        nir_current_spot = cleanSpot(cleanSpot(:,1)==currentCycle & cleanSpot(:,2)==7,:);
        regist7 = registerSpot(registerSpot(:,2)==7,:);
        regist359 = registerSpot(registerSpot(:,2)~=7 & registerSpot(:,2)>0,:);
        idx = rangesearch(nir_current_spot(:,[4,5]),regist7(:,[4,5]),r2);
        points_within_range = unique([idx{:}]);
        tmpSpot = nir_current_spot(setdiff(1:size(nir_current_spot,1),points_within_range), :);

        idx = rangesearch(tmpSpot(:,[4,5]),regist359(:,[4,5]),r1);
        points_within_range = unique([idx{:}]);
        nonOverlapSpots = tmpSpot(setdiff(1:size(tmpSpot,1),points_within_range), :);
        cSpot = size(nonOverlapSpots,1);
        totRegisterSpot = totRegisterSpot+cSpot;
        registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;

        txr_current_spot = cleanSpot(cleanSpot(:,1)==currentCycle & cleanSpot(:,2)==9,:);
        regist9 = registerSpot(registerSpot(:,2)==9,:);
        regist357 = registerSpot(registerSpot(:,2)~=9 & registerSpot(:,2)>0,:);
        idx = rangesearch(txr_current_spot(:,[4,5]),regist9(:,[4,5]),r2);
        points_within_range = unique([idx{:}]);
        tmpSpot = txr_current_spot(setdiff(1:size(txr_current_spot,1),points_within_range), :);

        idx = rangesearch(tmpSpot(:,[4,5]),regist357(:,[4,5]),r1);
        points_within_range = unique([idx{:}]);
        nonOverlapSpots = tmpSpot(setdiff(1:size(tmpSpot,1),points_within_range), :);
        cSpot = size(nonOverlapSpots,1);
        totRegisterSpot = totRegisterSpot+cSpot;
        registerSpot((totRegisterSpot-cSpot+1):totRegisterSpot,:) = nonOverlapSpots;
    end
    %remove the non-existed registerSpots
    alltempSpot = registerSpot(registerSpot(:,1)>0,:);


