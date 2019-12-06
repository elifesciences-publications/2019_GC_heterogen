function [spikeMat,threshold,aafr] = genspiketrain(GEN,MF)

%%% CREATE INPUT SPIKE TRAINS FOR MOSSY FIBERS %%%
% To control the bursty behaviour of the mossy fibers, the threshold-value
% which binarizes the randomreal spikeMat is dynamic. The function is
% adapting the threshold-value at every step, because the threshold is
% responsible for the probability that a spike occur


% create random start numbers individual for every spike train
% create a Matrix with ones at start numbers
randoms = zeros(GEN.nMF,MF.burstsmax);
startMat = zeros(GEN.nMF,GEN.nSteps);
% 1st half of MFs with bursts exp. decaying
for i = 1:GEN.nMF/2
    nBursts = randi([MF.burstsmin,MF.burstsmax]);
    randoms(i,1:nBursts) = randi(GEN.nSteps,[1,nBursts]);
    startMat(i,nonzeros(randoms(i,:))) = randi([MF.burstFreqMin,MF.burstFreqMax]);     %this threshold is the frequency in Hz later
end
clear i nBursts


% add exponential decay function to startMat
threshold = startMat;
for i = 1:GEN.nMF/2
    threshold(i,1) = 0;
    for j = 2:GEN.nSteps
        threshold(i,j) = threshold(i,j)+threshold(i,j-1) * exp(-GEN.dt/MF.tau_decay);
    end
end
clear i j

% 2nd half of MFs - rate coding with underlying gauss
for i = GEN.nMF/2+1 : GEN.nMF
    peakTime = rand * GEN.tSim; %in s
    peakAmp = sh_randWithinInterval(MF.afrMin,MF.afrMax);
    sig = sh_randWithinInterval(MF.gaussWidthMin,MF.gaussWidthMax);
    for j = 1:GEN.nSteps
        threshold(i,j) = peakAmp * exp(-0.5*((j*GEN.dt - peakTime) / sig)^2);
    end
end
clear i nBursts


% increase the probability of spikes where no decays happen
%threshold = threshold + MF.randomspikeprob;


% in order to retain the average firing rate the threshold function needs
% to be normalized
% for i = 1:GEN.nMF
%     threshold(i,:) = threshold(i,:) / (mean(threshold(i,:))/(MF.afr*(1-MF.afrRandFactor+2*MF.afrRandFactor*rand)*GEN.dt));
% end
% clear i


% create a matrix with random real number and binarize them with the values
% of the threshold function.
spikeMat = rand(GEN.nMF, GEN.nSteps) < threshold*GEN.dt;


% remove spikes when the maximum frequency is exceeded
% calculate the count of minimum steps to do before the next spike is
% allowed to occur
% if floor(1/(MF.fr_max*GEN.dt)) < 1
%     minsteps = 1;
% else
%     minsteps = floor(1/(MF.fr_max*GEN.dt));
% end
% 

% do the actual removing and save the positions where spikes were removed
% ind = zeros(size(spikeMat)); %optional for debugging
% deletioncounter = 0;
% for i = 1:GEN.nMF
%     col = find(spikeMat(i,:));
%     for j = 2:length(col)
%         if col(j) - col(j-1) < minsteps
%             spikeMat(i,col(j)) = 0;
%             deletioncounter = deletioncounter + 1;
%             col(j) = col(j-1);
%             %ind(i,col(j)) = 1; %optional for debugging
%         end
%     end
% end
% clear i j col
%[rows,cols] = find(ind); %optional for debugging
%delpos=[rows,cols*GEN.dt]; %optional for debugging


% calculate the actual average firing rate of every mossy fiber
aafr(1:GEN.nMF) = 0;
for i = 1:GEN.nMF
    aafr(i) = length(nonzeros(spikeMat(i,:))) / GEN.tSim;
end
clear i


% display the results
% if deletioncounter > 0
%     display('WARNING: Spikes were removed due to violation of the max-frequency-limit of mossy fibers.');
%     display(['The average firing rate was reduced by ' num2str(100 - mean(aafr) / MF.afr * 100) ' %']);
% else
%     display('0 spikes were removed due to violation of max-frequency-limit of mossy fibers.');
% end
% display([num2str(GEN.nMF) ' mossy fibers with an average firing rate of ' num2str(mean(aafr)) ' spikes per second have been generated.']);

end