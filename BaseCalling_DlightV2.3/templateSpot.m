

%calculate all the spots with defined algorithm
function [filterspot] = templateSpot(templateCycles,SC3,SC5,SNI,STX,imageSize,spotAlgorithm)
spotAll = zeros(size(templateCycles,2)*1000000,5);
totSpots = 0;
for i = templateCycles
    
  %generate transform matrix
  TPC5 = reshape(SC5(i,:,:),imageSize,imageSize);
  TPC3 = reshape(SC3(i,:,:),imageSize,imageSize);
  TPNI = reshape(SNI(i,:,:),imageSize,imageSize);
  TPTX = reshape(STX(i,:,:),imageSize,imageSize);
  
  %store the cycle and channel information
  spotCurrentC3 = spot_finding(TPC3,i,3,spotAlgorithm);
  spotCurrentC5 = spot_finding(TPC5,i,5,spotAlgorithm);
  spotCurrentNI = spot_finding(TPNI,i,7,spotAlgorithm);
  spotCurrentTX = spot_finding(TPTX,i,9,spotAlgorithm);
  CurrentStoreSpots = size(spotCurrentC3,1)+size(spotCurrentC5,1)+size(spotCurrentNI,1)+size(spotCurrentTX,1);
  
  totSpots=totSpots+CurrentStoreSpots;
  spotAll((totSpots+1):(totSpots+CurrentStoreSpots),:) = [spotCurrentC3;spotCurrentC5;spotCurrentNI;spotCurrentTX];
end
filterspot = spotAll(spotAll(:,3)>0,:);