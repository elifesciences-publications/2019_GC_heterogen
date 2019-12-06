function [spikeMatWin,meanpcerror2,temporalError] = myDiffReturnSpikes(syn, fminsearchAdd, ee1Trgt, ee2Trgt, nSteps, g_syn, I, E, Vspike, Vreset, dt, tau_m, Rm, E_e, Vth, WinSteps, pausePosBegin, pausePosEnd, tWin1, tWin2, temporalErrorScaling, d, trgtWin ,trgtcon)

% global GEN      % for general variables
% %global MF       % for mossy fiber variables
% %global SYN_MFGC % for synapse between mf and gc variables
% global GC       % for granule cell variables
% global SYN_GCPC % for synapse between gc and pc variables
% global PC       % for purkinje cell variables
% global LEARN    % for learning variables

% synWithFinalScale = synWithFinalScale - fminsearchAdd;
% 
% syn = synWithFinalScale;
% syn(length(syn))=[];%remove the final scale value again.
% syn = syn * synWithFinalScale(length(synWithFinalScale)); %multiply all values with final scale value.

syn = syn - fminsearchAdd;

g_syne = sum(repmat(syn',1,nSteps) .* g_syn);

% inhibitory input (PC.g_syni) for purkinje cell is zero for now
%PC.g_syni = zeros(1,GEN.nSteps);

% simulate the purkinje cell as LIF - model
% integrate the input and fire at a certain threshold (PC.V)
% binarize the output to get the spike train of the PC (PC.spikeMat)
% spike-rate adaptaion conductance (PC.gsra) and spike refractory
% conductance (PC.gref) should be zero for now
%[PC.V,PC.spikeMat,PC.gsra,PC.gref] = lif_sra_pc(GEN,PC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  lif_sra_pc START
VV = zeros(1,nSteps);        % Vector for membrane potential
spikes = false(1,nSteps);   % vector to record spike times
% gsra = zeros(1,GEN.nSteps);     % spike-rate adaptaion conductance
% gref = zeros(1,GEN.nSteps);     % spike refractory conductance

VV(:,1) = E;      % begin simulation of every GC with membrane potential at leak level
% gsra(:,1) = 0;      % initially no spike-rate adaptation
% gref(:,1) = 0;      % initially no refractory current


for i = 2:nSteps;
    % Adaptation and refractory current:
    % Next line solves GC.tau_sra*(d gsra/GEN.dt) = -gsra
    %    gsra(i) = gsra(i-1) * (1 - GEN.dt/PC.tau_sra);
    % Next line solves GC.tau_ref*(d gref/GEN.dt) = -gref
    %    gref(i) = gref(i-1) * (1 - GEN.dt/PC.tau_ref);
    
    if (VV(i-1) == Vspike)      % Last datapoint was a spike?
        VV(i) = Vreset;         % reset after spike
    else                            % otherwise integrate
        VV(i) = VV(i-1)...                            % Last voltage...
            +dt/tau_m*...
            (E-VV(i-1)...                           % +membrane timeconstant relaxation to resting membrane
            +Rm*I(i-1)...                            % +GC.Rm * currentinjection
            -(VV(i-1)-E_e)*Rm*g_syne(i-1)...% -excitaroty synaptic conductance * GC.Rm * driving force
            );
    end
    
    if VV(i) > Vth       % if voltage is greater than threshold spike
        VV(i) = Vspike;     % change voltage so we see the spike
        spikes(i) = true;        % record the time of the spike
%        gsra(i) = gsra(i-1) + PC.delta_gsra;   % increase the spike-rate-adaptation conductance
%        gref(i) = gref(i-1) + PC.delta_gref;   % increase the refractory conductance
    end
end
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% cut out the window
%PC.spikeMatWin = PC.spikeMat(GEN.WinSteps);
spikeMatWin = spikes(WinSteps);



%%%%%%%%%%%%%%%%%%%%%%%%%% Temoral Error before and after pause %%%%%%%%%%%%%%%%%
%[ee1,ee2] = getTemporalError(GEN,LEARN,PC.spikeMatWin);
%disp(['ee1=',num2str(ee1),'  ee2=',num2str(ee2),'   ee1Trgt=',num2str(ee1Trgt),'  ee2Trgt=',num2str(ee2Trgt)]);
pcpos = find(spikeMatWin);

pausePosMiddle = round((pausePosBegin + pausePosEnd)/2);

if isempty(find(pcpos < pausePosMiddle)) == 1
    ee1 = tWin1;
else
    ee1 = dt*pcpos(max(find(pcpos < pausePosMiddle)));
end

if isempty(find(pcpos >= pausePosMiddle)) == 1
    ee2 = tWin2;
else
    ee2 = dt*pcpos(min(find(pcpos >= pausePosMiddle)));
end



temporalError = temporalErrorScaling * (abs(ee1-ee1Trgt)+abs(ee2-ee2Trgt));

%%%%%%%%%%%%%%%%%%%%%%%%%% BUILT ROSSUM %%%%%%%%%%%%%%%%%

pccon = conv2(double(spikeMatWin),d,'full');
% pccon(end-size(d,2)+2:end) =[];

% trgtcon = conv(double(trgtWin),d);
% trgtcon(end-size(d,2)+2:end) =[];

diffcon = (pccon - trgtcon).^2;

meanpcerror2 = sum(diffcon);

%display(['   Rossum = ',num2str(LEARN.meanpcerror2),'   temporal error = ',num2str(LEARN.temporalError),]);
        
%%%%%%%%%%%%%%%%%%%%%%%%%% Sum of temporal error and ROSSUM %%%%%%%%%%%%%%%%%
x = meanpcerror2;% + temporalError;
end
