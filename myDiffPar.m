function x = myDiffPar(syn, K)
% K: Konstanten

nPar = size(syn,1);

syn = syn - K.LEARN.fminsearchAdd;

if 0 %K.GEN.nGC >= 2000  
    g_syne = gather(gpuArray(syn) * gpuArray(K.SYN_GCPC.g_syn))';      % @GPU 
else
    g_syne = (syn * K.SYN_GCPC.g_syn)';                                % @CPU 
end

VV = zeros(K.GEN.nSteps,nPar);          % vector for membrane potential
spikes = false(K.GEN.nSteps,nPar);      % vector to record spike times
VV(1,:) = K.PC.E;                       % begin simulation of every GC with membrane potential at leak level

% parfor 
nSteps = K.GEN.nSteps;
vSpike = K.PC.Vspike;
vReset = K.PC.Vreset;
k1 = K.k1;
k2 = K.k2;
k3 = K.k3;
k4 = K.k4;
I = K.I;
vTh = K.PC.Vth;

parfor jj = 1:nPar
    vvPar = VV(:,jj);
    spikesPar = spikes(:,jj);
    gSynePar = g_syne(:,jj);
    for ii = 2:nSteps;

        if vvPar(ii-1) == vSpike         % Last datapoint was a spike?
            vvPar(ii) = vReset;          % reset after spike
        else                             % otherwise integrate
            vvPar(ii) = vvPar(ii-1) * (k2-k3*gSynePar(ii-1)) + k1 + k3 * I(ii-1) + k4 * gSynePar(ii-1);

            if vvPar(ii) > vTh           % if voltage is greater than threshold spike
                vvPar(ii) = vSpike;      % change voltage so we see the spike
                spikesPar(ii) = true;    % record the time of the spike
            end
        end
    end
     VV(:,jj) = vvPar;
     spikes(:,jj) = spikesPar;
end

spikes = spikes';
spikeMatWin = double(spikes(:,K.GEN.WinSteps));   % cut out the window

%%%%%%%%%%%%%%%%%%%%%% Temporal Error %%%%%%%%%%%%%%%%%
% nMatWin = size(spikeMatWin,2);
% pcpos = zeros(1,nMatWin);
% temporalError = zeros(nPar,1);
% 
% parfor n = 1:nPar
%     pcpos = find(spikeMatWin(n,:));
% 
%     temp = find(pcpos < K.pausePosMiddle,1,'last');     % max
%     if isempty(temp) == 1
%         ee1 = K.GEN.tWin(1);
%     else
%         ee1 = K.GEN.dt*pcpos(temp);
%     end
% 
%     temp = find(pcpos >= K.pausePosMiddle,1,'first');   % min
%     if isempty(temp) == 1
%         ee2 = K.GEN.tWin(2);
%     else
%         ee2 = K.GEN.dt*pcpos(temp);
%     end
%     temporalError(n) = K.LEARN.temporalErrorScaling * (abs(ee1-K.ee1Trgt)+abs(ee2-K.ee2Trgt));
% end

%%%%%%%%%%%%%%%%%%%%%%%%%% VAN ROSSUM %%%%%%%%%%%%%%%%%
pccon = conv2(spikeMatWin,K.dd,'full'); 
%diffcon = (pccon(:,1:end-size(K.dd,2)+1) - repmat(K.LEARN.trgtcon,[nPar,1]));
diffcon = (pccon - repmat(K.LEARN.trgtcon,[nPar,1]));
meanpcerror2 = sum(diffcon.^2,2);
      
%%%%%%% Sum of temporal error and ROSSUM %%%%%%%%%%%%%%
x = meanpcerror2;% + temporalError;
end
    

