



function [topCof] = weightedHist(pairedColorLeft,pairedColorRight,breaks,minPP)

  [polTheta, polR] = cart2pol(pairedColorLeft,pairedColorRight);
  
  minPhi = min(polTheta);
  maxPhi = max(polTheta);
  
  interval = (maxPhi-minPhi)/breaks;

  bins = minPhi:interval:maxPhi;
  
  weightedbin = zeros(size(bins,2),1);
  for thetaInd = 1:size(polTheta,1)
    for i = 1:(size(bins,2)-1)
        if(polTheta(thetaInd) >= bins(i) && polTheta(thetaInd) <bins(i+1))
            weightedbin(i) = weightedbin(i)+polR(thetaInd);
        end
    end
  end
  
  %sort the weighted bin to find the two maximum value
  [pks,locs] = findpeaks(weightedbin,'MinPeakProminence',minPP);
  select = [pks,locs];
  %sortedPeak = sortrows(select,1,'descend');
  %cutoff = sortedPeak(1,1)*0.3;
  upperPeak = select(select(:,1)>5,:);
  
  topCof = tan(bins(upperPeak(1,2)));

  
  %topCof = bins(upperPeak(:,2));
  %topCof =sortedPeak;
  
  