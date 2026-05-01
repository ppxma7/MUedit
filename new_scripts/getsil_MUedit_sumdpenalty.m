function sil = getsil_MUedit_sumdpenalty(PulseT, fsamp)

% find number of windows
used = PulseT ~= 0; % logical vector: 1 where signal is nonzero
d = diff([0, used, 0]);  % pad with zeros at both ends
startIdx = find(d == 1);
nrwin = length(startIdx);

[~,spikes] = findpeaks(PulseT,'MinPeakDistance',round(fsamp*0.005)); % 4b: Peak detection

maxSpikes = maxk(PulseT(spikes),8*nrwin);
PulseT = PulseT/mean(maxSpikes(nrwin*2+1:8*nrwin)); % normalization of the MU pulse train
spikes = spikes(PulseT(spikes)<1.2);
% [L,C,sumd,D] = kmeans(PulseT(spikes)',2, 'Replicates', 20); % 4c: Kmean classification
[L,C,sumd,D] = kmeans(PulseT(spikes)',2);
[~, idx2] = max(C); % Spikes should be in the class with the highest centroid
within = sumd(idx2);
between = sum(D(L==idx2, setdiff([1 2],idx2)));
sil = (between-within)/max([within,between]); % Silhouette measure
if sil == 1 && any(sumd == 0) && any(histcounts(L)<2)
    sil = 0;
end%if
if sumd(~ismember([1,2],idx2))>1*sumd(idx2) || abs(diff(C))<0.5
    sil = 0;
end%if