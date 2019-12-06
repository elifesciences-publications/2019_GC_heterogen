function [fval,rossErr,tempErr,rossErrDuringLearn,finalPCSpikes] = main(runthrough,CTRL)
%% DECLARATION

% GEN = struct;       % structure array for general variables
% MF = struct;        % structure array for mossy fiber variables
% SYN_MFGC = struct;  % structure array for synapse between mf and gc variables
% GC = struct;        % structure array for granule cell variables
% SYN_GCPC = struct;  % structure array for synapse between gc and pc variables
% PC = struct;        % structure array for purkinje cell variables
% LEARN = struct;     % structure array for learning variables
% PLOT = struct;      % structure array for plotting variables

%% CONTENT

% 1. INPUT - MODEL SETTINGS
% 2. INPUT - PHYSIOLOGICAL SETTING (via config_physio)
% 3. GENERATING AND PREPROCESSING DATA
% 4. LEARNING LOOP
% 5. PLOTTING

%% 1. INPUT - MODEL SETTINGS
%%% INPUT - MODEL SETTINGS - GENERAL %%%

GEN.tSim = 0.5;                 % simulated time (seconds)
GEN.tWin = [0.1 0.5];           % time window which is used by learning algorithm (seconds)
GEN.dt = CTRL.dt;               % time step in seconds
GEN.nGC = CTRL.nGC;             % number of granule cells
GEN.mfpergc = CTRL.mfpergc;     % how many mossy fibers per granule cell
GEN.gcpermf = 10;               % how many granule cells can be connected to one mossy fiber
GEN.pause = [0.3 0.4];          % pause

%%% INPUT - MODEL SETTINGS - LEARN

% GA vs. LERN, timeLimit: LEARN.trails = 10000 =^ 5.5 s 
LEARN.trials = 2000;                        % how many trials should be evaluated
LEARN.tau_rise = 1e-3;                      % exponential rise timeconstant for error function
LEARN.tau_decay = 10e-3;                    % exponential decay timeconstant for error function
%LEARN.tau_decay_rossum = CTRL.tauRossum;    % exponential decay timeconstant for error function
LEARN.temporalErrorScaling = 10^4;          % i.e. no temporal weight
LEARN.fminsearchAdd = 0;
LEARN.maxIter = 500;
LEARN.numberOfOpti = CTRL.nOpti;
LEARN.RandomChangeInPercent = 10;

%init rossErrDuringLearn
%rossErrDuringLearn = zero(LEARN.trials);
rossErrDuringLearn = 1:LEARN.trials;
rossErrDuringLearn = 0*rossErrDuringLearn;
temporalErrDuringLearn = 1:LEARN.trials;
temporalErrDuringLearn = 0*temporalErrDuringLearn;

%% 2. INPUT - PHYSIOLOGICAL SETTING (via config_physio)
%randomseeds = 2000004:2001999;
%rng(randomseeds(runthrough));   % fixes the seed for random numbers
%rngSetting = rng;
config_physio;  % outsourced config file for physiological settings
%rng(rngSetting);   % set back rng 

%% 3. MAKE GC spiking and EPSPs - GENERATING AND PREPROCESSING DATA
% general calculation
GEN.nMF = GEN.nGC / GEN.gcpermf * GEN.mfpergc;                  % how many mossy fibers needed at least
GEN.nSteps = floor(GEN.tSim/GEN.dt);                            % count of steps depends on simulation time and dt
GEN.WinSteps = ceil(GEN.tWin(1)/GEN.dt)+1:ceil(GEN.tWin(2)/GEN.dt);   % give only steps in the Window
GEN.tVecSim = GEN.dt:GEN.dt:GEN.tSim;                           % create a vector with all time steps in tSim excluding 0
GEN.tVecWin = GEN.tWin(1)+GEN.dt:GEN.dt:GEN.tWin(2);            % create a vector with all time steps in tWindow

disp(['The number of MFs is ',num2str(GEN.nMF)]);

% generate simulated mossy fiber spike trains (MF.spikeMat) and the
% underlying threshold function (MF.threshold) and measure the actual
% average firing rate of every mossy fiber (MF.aafr)
[MF.spikeMat,MF.threshold,MF.aafr] = genspiketrain(GEN,MF);

% simulate the synapse at the end of mossy fibers
% convert mossy fiber spike trains to EPSPs (SYN_MFGC.g_syn)
[SYN_MFGC.g_syn] = synapse_mfgc(GEN,SYN_MFGC,MF);

% wire mossy fibers to granule cells in a random fashion and save the
% indices of the selected mossy fibers for each granule cell (GC.mfsel)
GC.mfsel(GEN.nGC,GEN.mfpergc) = 0;
for n = 1:GEN.nGC
    for m = 1:GEN.mfpergc
        GC.mfsel(n,m)=randi(GEN.nMF);
    end
end

% sum up the selected mossy fiber traces in order to get the excitatory
% input for each granule cell (GC.g_syne)
for n = 1 : GEN.nGC
    GC.g_syne(n,:) = sum(...
        [SYN_MFGC.g_syn(GC.mfsel(n,:),:); zeros(1,size(SYN_MFGC.g_syn,2))]);
end

% inhibitory input (GC.g_syni) for granule cells is zero for now
GC.g_syni = zeros(GEN.nGC,GEN.nSteps);

% simulate the granule cell as LIF - model
% integrate the input of GCs and fire at a certain threshold (GC.V)
% binarize the output to get spike trains of granule cells (GC.spikeMat)
% spike-rate adaptaion conductance (GC.gsra) and spike refractory
% conductance (GC.gref) should be zero for now
[GC.V,GC.spikeMat,GC.gsra,GC.gref] = lif_sra_gc(GEN,GC);

% ###############################################################
%   spikeShift
if CTRL.spikeShift == 1
    spikeshift = round(linspace(0,GC.spikeShift,GEN.nGC));

    for n = 1:GEN.nGC
        temp = GC.spikeMat(n,:);
        temp = circshift(temp',spikeshift(n))';
        if spikeshift(n) > 0
           temp(1:spikeshift(n)) = 0; 
        end
        GC.spikeMat(n,:) = temp;
    end

    if CTRL.spikeMatMix == 1
        rngSetting = rng;
        mixInd = randperm(GEN.nGC);
        GC.spikeMat = GC.spikeMat(mixInd,:);
        rng(rngSetting);  
    end
end

% ###############################################################

% ###############################################################
% impose maximal frequency for each GC 
% remove spikes when the maximum frequency is exceeded
% calculate the count of minimum steps to do before the next spike is
% allowed to occur

% do the actual removing and save the positions where spikes were removed
%ind = zeros(size(spikeMat)); %optional for debugging
deletioncounter = 0;      
for n = 1:GEN.nGC
    if floor(1/(GC.fr_max(n)*GEN.dt)) < 1
        minsteps = 1;
    else
        minsteps = floor(1/(GC.fr_max(n)*GEN.dt));
    end
    col = find(GC.spikeMat(n,:));
    for j = 2:length(col)
        if col(j) - col(j-1) < minsteps
            GC.spikeMat(n,col(j)) = 0;
            deletioncounter = deletioncounter + 1;
            col(j) = col(j-1);
        end
    end
end
% ###############################################################

% simulate the synapse at the end of granule cells
% convert granule cell spike trains to EPSPs (SYN_GCPC.g_syn)
[SYN_GCPC.g_syn] = synapse_gcpc(GEN,SYN_GCPC,GC.spikeMat);

% cut the window out of GC spiketrains (GC.spikeMatWin)
GC.spikeMatWin = GC.spikeMat(:,GEN.WinSteps);

% exclude GCs which do not fire during the windows (no EPSP trace in window)
% since they are not firing, no weight should be attached to them
% GCs that have a spike in window (fireGCs)
fireGCs = find(any(GC.spikeMatWin,2));
disp(['Warning: ',num2str(GEN.nGC-length(fireGCs)),...
    ' GCs are not firing during learning time-window.']);

% initial weights (LEARN.gcweight) and initial total weights
% (LEARN.gcweightstotal) are random between 0 and 2 for firing GCs and 0
% for not-firing GCs
LEARN.gcweightstotal = nan(1,GEN.nGC);
for n = 1 : GEN.nGC
    LEARN.gcweightstotal(n) = 1;
end

LEARN.gcweightstotal(setdiff(1:GEN.nGC,fireGCs)) = 0;

% for checking of MF and GC connetivity and EPSPs
% plot_mf(1,1,GEN,PLOT,MF,SYN_MFGC,GC);
% plot_gc(2,1,GEN,PLOT,SYN_MFGC,SYN_GCPC,GC);
% return;            

%% GENERATE TARGET
% make target random ------------------
LEARN.gcweightstarget = nan(1,GEN.nGC);
for n = 1 :GEN.nGC
    LEARN.gcweightstarget(n) = 2 * rand;
end
%clear i
LEARN.gcweightstarget(setdiff(1:GEN.nGC,fireGCs)) = 0;
% do the same procedure for generating the pc output like in the learning
% loop but save the spike train to "LEARN.trgt"
% the rest will be overwritten during the first trial
PC.g_syne = LEARN.gcweightstarget * SYN_GCPC.g_syn;
% PC.g_syne = sum(repmat(LEARN.gcweightstarget',1,GEN.nSteps) .* ...
%                       SYN_GCPC.g_syn);
PC.g_syni = zeros(1,GEN.nSteps);
[~,LEARN.trgt,~,~] = lif_sra_pc(GEN,PC);
LEARN.trgtWin = LEARN.trgt(GEN.WinSteps);

% make artificial target
LEARN.freq1 = 80;
LEARN.freq2 = 120;
LEARN.trgt = 0 * LEARN.trgt;
LEARN.trgt(1:round(1/(LEARN.freq1*GEN.dt)):round(((GEN.pause(1)+GEN.pause(2))/2)/GEN.dt)) = 1;
LEARN.trgt(round(((GEN.pause(1)+GEN.pause(2))/2)/GEN.dt)+1:round(1/(LEARN.freq2*GEN.dt)):length(LEARN.trgt)) = 1;
LEARN.trgt = logical(LEARN.trgt);

LEARN.trgtWin = LEARN.trgt(GEN.WinSteps);
    
% make Pause ------------------
GEN.pausePosBegin = 1+round(GEN.pause(1)/GEN.dt);
GEN.pausePosEnd = 1+round(GEN.pause(2)/GEN.dt);
% cut the window out of target-spiketrain (LEARN.trgtWin)
LEARN.trgt(GEN.pausePosBegin : GEN.pausePosEnd) = 0;
LEARN.trgtWin = LEARN.trgt(GEN.WinSteps);

% calculate pause begin and end again, but not for the window used for learning within the loop
% values will be used in getTemoralError
GEN.pausePosBegin = 1+round((GEN.pause(1)-GEN.tWin(1))/GEN.dt);
GEN.pausePosEnd = 1+round((GEN.pause(2)-GEN.tWin(1))/GEN.dt);
% The values for the target have to be calulated only once. They should no
% change.
[ee1Trgt,ee2Trgt] = getTemporalError(GEN,LEARN,LEARN.trgtWin);

%% 4. Calulate PC spike sequence once to define all global variables
%     and calcualte the rossum and temporal error

% attach weight to granule cell traces and sum them up in order to get the
% excitatory input for the purkinje cell (PC.g_syne)
PC.g_syne = sum(repmat(LEARN.gcweightstotal',1,GEN.nSteps) .* ...
    SYN_GCPC.g_syn);

% inhibitory input (PC.g_syni) for purkinje cell is zero for now
PC.g_syni = zeros(1,GEN.nSteps);

% simulate the purkinje cell as LIF - model
% integrate the input and fire at a certain threshold (PC.V)
% binarize the output to get the spike train of the PC (PC.spikeMat)
% spike-rate adaptaion conductance (PC.gsra) and spike refractory
% conductance (PC.gref) should be zero for now
[PC.V,PC.spikeMat,PC.gsra,PC.gref] = lif_sra_pc(GEN,PC);

% cut out the window
PC.spikeMatWin = PC.spikeMat(GEN.WinSteps);

% Temporal Error before and after pause 
[LEARN] = geterror(GEN,LEARN,PC.spikeMatWin);
[ee1,ee2] = getTemporalError(GEN,LEARN,PC.spikeMatWin);
%disp(['ee1=',num2str(ee1),'  ee2=',num2str(ee2),'   ee1Trgt=',num2str(ee1Trgt),'  ee2Trgt=',num2str(ee2Trgt)]);
LEARN.temporalError = LEARN.temporalErrorScaling * (abs(ee1-ee1Trgt)+abs(ee2-ee2Trgt));

% calculate ROSSUM error
imp = zeros(1,int64((CTRL.tauRossum/GEN.dt)*3));
imp(2) = 1;
dd = zeros(size(imp));
for n=2:size(imp,2)
    dd(n) = imp(n)+dd(n-1)*exp(-GEN.dt/CTRL.tauRossum);
end
dd(1) = [];
CTRL.dd = dd;

% target is always the same therefore this must be calulated only once:
LEARN.trgtcon = conv2(double(LEARN.trgtWin),CTRL.dd,'full');
% 
% LEARN.pccon = conv2(double(PC.spikeMatWin),CTRL.dd,'full');
% LEARN.diffcon = (LEARN.pccon - LEARN.trgtcon).^2;
% LEARN.meanpcerror2(1) = sum(LEARN.diffcon);

% calculate ROSSUM error
[LEARN.meanpcerror2] = vanRossumErr(PC.spikeMatWin,CTRL.dd,LEARN.trgtcon);

%stuff needed only once:
% apply current during a pulse
I = zeros(size(GEN.tVecSim));         % Vector for currents
for n = PC.tpulse/GEN.dt : (PC.tpulse+PC.lengthpulse)/GEN.dt
    I(round(n)) = PC.Iapp;
end

% Constants for optiAlgo
K.LEARN.fminsearchAdd = LEARN.fminsearchAdd;
K.ee1Trgt = ee1Trgt;
K.ee2Trgt = ee2Trgt;
K.GEN.nSteps = GEN.nSteps;
K.SYN_GCPC.g_syn = SYN_GCPC.g_syn;
K.I = I;
K.PC.E = PC.E;
K.PC.Vspike = PC.Vspike;
K.PC.Vreset = PC.Vreset;
K.GEN.nGC = GEN.nGC;
K.PC.Vth = PC.Vth;
K.GEN.WinSteps = GEN.WinSteps;
K.GEN.pausePosBegin = GEN.pausePosBegin;
K.GEN.pausePosEnd = GEN.pausePosEnd;
K.GEN.tWin(1) = GEN.tWin(1);
K.GEN.tWin(2) = GEN.tWin(2);
K.LEARN.temporalErrorScaling = LEARN.temporalErrorScaling;
K.dd = CTRL.dd ;
K.LEARN.trgtWin = LEARN.trgtWin;
K.LEARN.trgtcon = LEARN.trgtcon;
K.pausePosMiddle = round((K.GEN.pausePosBegin + K.GEN.pausePosEnd )/2);

a = GEN.dt;
b = PC.tau_m;
c = PC.E;
d = PC.Rm;
e = PC.E_e;

K.k1 = (a*c)/b;
K.k2 = 1-(a/b);
K.k3 = (a*d)/b;
K.k4 = (a*d*e)/b;

% parallel computing
switch CTRL.par
    
    case 0
        fun = @(x)myDiff(x,K);
    case 1
        fun = @(x)myDiffPar(x,K);
end

%% Learning based on plasticity rule
disp('START pre-Optimization');

if strcmp(CTRL.opti, 'LEARN')
    tl = tic;
    disp('### Learning based on plasticity rule Optimization ###');
    LEARN.gcweightstotalBest=LEARN.gcweightstotal;
    LEARN.bestErr = 1E6;

    for tr = 1:LEARN.trials

        % attach weight to granule cell traces and sum them up in order to get the
        % excitatory input for the purkinje cell (PC.g_syne)
        PC.g_syne = LEARN.gcweightstotal * SYN_GCPC.g_syn;
        % PC.g_syne = sum(repmat(LEARN.gcweightstotal',1,GEN.nSteps) .* ...
        % SYN_GCPC.g_syn);

        % inhibitory input (PC.g_syni) for purkinje cell is zero for now
        % PC.g_syni(tr,:) = zeros(1,GEN.nSteps);

        % simulate the purkinje cell as LIF - model
        % integrate the input and fire at a certain threshold (PC.V)
        % binarize the output to get the spike train of the PC (PC.spikeMat)
        % spike-rate adaptaion conductance (PC.gsra) and spike refractory
        % conductance (PC.gref) should be zero for now
        [PC.V,PC.spikeMat,PC.gsra,PC.gref] = lif_sra_pc(GEN,PC);

        % cut out the window
        PC.spikeMatWin = PC.spikeMat(GEN.WinSteps);

        % determine the position error (distance from every pc spike to the nearest
        % target spike) and convolve it with the error function
        %[LEARN] = errorfunction(GEN,LEARN,PC.spikeMat(tr,:));
        [LEARN] = geterror(GEN,LEARN,PC.spikeMatWin);

        % calculate the position error for each granule cell in a few steps
        % 1. looking up the value in the error function where spikes of a granule cell occur
        % 2. sum up these values and divide the result by the count of values(mean)
        % 3. if there is no value the result of step 2 is NaN, exchange NaN with 0
        %GCerrors = LEARN.errorfct*GC.spikeMatWin;
        GCerrors = repmat(LEARN.errorfct,GEN.nGC,1).*GC.spikeMatWin;
        LEARN.gcposerror = sum(GCerrors,2)./sum(GCerrors~=0,2);
        LEARN.gcposerror(isnan(LEARN.gcposerror)) = 0;

        LEARN.gcposerror = LEARN.gcposerror';

        LEARN.gcposerror = LEARN.gcposerror .* LEARN.gcweightstotal;

        LEARN.gcposerrorfactor = 1 ./ ((LEARN.gcposerror/1) + 1);

        % calculate the ratio between the number of target spikes to output spikes
        LEARN.appc = length(nonzeros(PC.spikeMatWin));
        LEARN.apratio = ((1+length(nonzeros(LEARN.trgtWin))) / ...
            (1+LEARN.appc));

        LEARN.apratiofactor = LEARN.apratio^(0.05);

        % put them together
        LEARN.gcweightstotalTemp = LEARN.gcweightstotal;
        LEARN.gcweightstotal = LEARN.gcweightstotal .* (LEARN.gcposerrorfactor .* LEARN.apratiofactor);

        % calculate the average position error of all GCs per trial
        LEARN.meangcposerror = mean(LEARN.gcposerror);

        % calculate the mean error for every trial
        LEARN.meanpcerror = mean(LEARN.mindistanceVec);
        % Rossum-Error        
        [LEARN.meanpcerror2] = vanRossumErr(PC.spikeMatWin,CTRL.dd,LEARN.trgtcon);
        
        % Temoral Error before and after pause
        [ee1,ee2] = getTemporalError(GEN,LEARN,PC.spikeMatWin);
        %disp(['ee1=',num2str(ee1),'  ee2=',num2str(ee2),'   ee1Trgt=',num2str(ee1Trgt),'  ee2Trgt=',num2str(ee2Trgt)]);
        pcpos = find(PC.spikeMatWin);
    
        pausePosMiddle = round((GEN.pausePosBegin + GEN.pausePosEnd)/2);
    
        if isempty(find(pcpos < pausePosMiddle)) == 1
            ee1 = GEN.tWin(1);
        else
            ee1 = GEN.dt*pcpos(max(find(pcpos < pausePosMiddle)));
        end
    
        if isempty(find(pcpos >= pausePosMiddle)) == 1
            ee2 = GEN.tWin(2);
        else
            ee2 = GEN.dt*pcpos(min(find(pcpos >= pausePosMiddle)));
        end

        temporalError = LEARN.temporalErrorScaling * (abs(ee1-ee1Trgt)+abs(ee2-ee2Trgt));
        
        if ~mod(tr,200) || tr == 2
            display([num2str(runthrough), ' | ', num2str(tr),'   Best fval = ' num2str(LEARN.bestErr,'%4.1f')]);
        end
        rossErrDuringLearn(tr)=LEARN.meanpcerror2;
        
        if CTRL.show == 1
            plot(rossErrDuringLearn);     
            ylim([0 5000]);
            title(['LEARN | Best fval = ' num2str(LEARN.bestErr,'%4.1f')]);
            ylabel('err');
            xlabel('try');
            grid on;
            drawnow;
        end
        
        temporalErrDuringLearn(tr)=temporalError;

        if (LEARN.meanpcerror2 + 0*temporalError) < LEARN.bestErr
            LEARN.bestErr = (LEARN.meanpcerror2 + 0*temporalError);
            %LEARN.gcweightstotalBest=LEARN.gcweightstotal;
            LEARN.gcweightstotalBest=LEARN.gcweightstotalTemp;
        end   
    end

    tl = toc(tl);
    disp(['Time: ' num2str(tl,'%4.1f') ' s']);
end

%% Optimization
disp('START Optimization');

%intial valus
syn0 = LEARN.fminsearchAdd + LEARN.gcweightstotalBest; % starting values

% save startvalues
LEARN.gcweightstotalStart = syn0;

fval = fun(syn0);
   
funReturnSpikes = @(x)myDiffReturnSpikes(x, LEARN.fminsearchAdd,...
    ee1Trgt, ee2Trgt, GEN.nSteps, SYN_GCPC.g_syn, I, PC.E,...
    PC.Vspike, PC.Vreset, GEN.dt, PC.tau_m, PC.Rm, PC.E_e,...
    PC.Vth, GEN.WinSteps, GEN.pausePosBegin, GEN.pausePosEnd,...
    GEN.tWin(1), GEN.tWin(2), LEARN.temporalErrorScaling, dd,...
    LEARN.trgtWin, LEARN.trgtcon);

[PC.spikeMatWin,LEARN.meanpcerror2,LEARN.temporalError] = funReturnSpikes(syn0);

CTRL.ub = 100;
CTRL.lb = 0;

UB = zeros(GEN.nGC,1) + CTRL.ub;
LB = zeros(GEN.nGC,1) + CTRL.lb; % take care that synaptic weights are not getting negative

% now change start values and do it again
for ooo = 1 : LEARN.numberOfOpti
    
    syn0Before = syn0;
    fvalBefore = fval;
    if ooo > 1      % erstmal beim "LEARN-Optimum" suchen
        for iii = 1 : GEN.nGC           %i.e. the last value of syn0 (scaling of weights) will not be changed rrandomly. It is too senistive.
            syn0(iii) = LEARN.fminsearchAdd +...
                (syn0(iii) - LEARN.fminsearchAdd)...
                * sh_randWithinInterval(1-LEARN.RandomChangeInPercent/100,1+LEARN.RandomChangeInPercent/100);
        end
    end
        
    opts.Display = 'off';
%   opts.TimeLimit = CTRL.preOptiTimeLimit;
    opts.TolMesh = CTRL.psTol;
    opts.TolX = CTRL.psTol;
    opts.TolFun = CTRL.psTol;
    %opts.MaxFunEvals = 1000;
    %opts.MaxIter = 100;
     
    if CTRL.par
        opts.UseParallel = 'always';
        opts.CompletePoll = 'on';
        opts.Vectorized = 'on';
        %opts.MaxIter = 5;
    end

    if CTRL.show == 1
        opts.PlotFcns{1} = @psplotbestf;
        opts.PlotFcns{2} = @psplotfuncount;
    end

    tps = tic;
    [syn0,fval,eFlag,out] = patternsearch(fun,syn0,[],[],[],[],LB,UB,[],opts);
    tps = toc(tps);

    [PC.spikeMatWin,LEARN.meanpcerror2,LEARN.temporalError] = funReturnSpikes(syn0);

    disp(['# of Opti = ',num2str(ooo), '  fval = ' num2str(fval,'%4.1f'),...
        '  Ros. = ',num2str(LEARN.meanpcerror2,'%4.0f'),...
        'tempErr = ',num2str(LEARN.temporalError,...
        '%4.0f'), '  t = ' num2str(tps,'%2.1f') 's']);
    
    if fval > fvalBefore
        syn0 = syn0Before;
        fval = fvalBefore;
        disp('CAVE: fval did not decrease. Reset to previous value.');
    end

end

disp(['# Best', ' fval = ' num2str(fval,'%4.1f')]);
  
%% rest
fval = fun(syn0);
[PC.spikeMatWin,rossErr,tempErr] = funReturnSpikes(syn0);
finalPCSpikes = PC.spikeMatWin;

close all;
sound(1);
end