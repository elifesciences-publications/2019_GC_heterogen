%% INPUT - GENERATION OF MOSSY FIBER SPIKE TRAINS
%MF.afr = 30;                % average firing rate (Hz)
%MF.afrRandFactor = 0.9;

%only for 50% of the MFs
MF.afrMin = 10;                % average firing rate (Hz)
MF.afrMax = 100;                % average firing rate (Hz)
MF.gaussWidthMin = 0.2;     % min width if the gauss for increased firing (s)
MF.gaussWidthMax = 0.5;      % min width if the gauss for increased firing (s)

%only for 50% of the MFs
MF.burstsmin = 1;           % how many bursts will occur at least during tSim
MF.burstsmax = 1;           % how many bursts will occur at most during tSim
MF.burstFreqMin = 600;           % maximum frequency during a burst of spikes (Hz)
MF.burstFreqMax = 1200;           % maximum frequency during a burst of spikes (Hz)
%MF.randomspikeprob = 0.02; % increase random spike occurrence and decrease weight of threshold function
MF.tau_decay = 30e-3;      % timeconstant for the decay of spike probability (seconds)


%% INPUT - SYNAPSE BETWEEN MOSSY FIBER AND GRANULE CELL
SYN_MFGC.g_max = 2.0e-8;      % maximal conductance of the synapse (Siemens) %% 80pA/80mV

SYN_MFGC.p_base = 0.5;      % baseline probability of release
SYN_MFGC.n_base = 0.3;
SYN_MFGC.f = 0.2;           % facilitation per spike
SYN_MFGC.tau_f = 0.012;      % timeconstant for facilitation (seconds)
SYN_MFGC.tau_d = 0.013;       % timeconstant for depression (seconds)

SYN_MFGC.tau_rise = 0.02e-3; % timeconstant for rise
SYN_MFGC.tau_decay = 2e-3; % timeconstant for decay


%% INPUT - LIF GRANULE CELL
GC.E = -80e-3;         % resting membrane potential (V)

switch CTRL.mode
    case 'homogen'
        GC.Rm(1:GEN.nGC) = 625e6; %linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        GC.Vth(1:GEN.nGC) = -39e-3;%linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        GC.Vreset(1:GEN.nGC) = -90e-3;    % AP reset potential (V)
        GC.fr_max(1:GEN.nGC) = 500;     % maximal frequency
    case 'homogenR'
        GC.Rm(1:GEN.nGC) = 800e6; %linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        GC.Vth(1:GEN.nGC) = -41e-3;%linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        GC.Vreset(1:GEN.nGC) = -110e-3;    % AP reset potential (V)
        GC.fr_max(1:GEN.nGC) = 1000;     % maximal frequency
    case 'homogenL'        
        GC.Rm(1:GEN.nGC) = 450e6; %linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        GC.Vth(1:GEN.nGC) = -37e-3;%linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        GC.Vreset(1:GEN.nGC) = -80e-3;    % AP reset potential (V)
        GC.fr_max(1:GEN.nGC) = 100;     % maximal frequency
    case 'heterogen'
        GC.Rm(1:GEN.nGC) = linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        GC.Vth(1:GEN.nGC) = linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        GC.Vreset(1:GEN.nGC) = linspace(-80e-3,-110e-3,GEN.nGC); % AP reset potential (V)
        GC.fr_max(1:GEN.nGC) = linspace(1000,100,GEN.nGC); % maximal frequency
    case 'heterogenMix'
        mixInd = randperm(GEN.nGC);
        temp1 = linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        temp2 = linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        temp3 = linspace(-80e-3,-110e-3,GEN.nGC); % AP reset potential (V)
        temp4 = linspace(1000,100,GEN.nGC); % maximal frequency

        GC.Rm = temp1(mixInd);    
        GC.Vth = temp2(mixInd);   
        GC.Vreset = temp3(mixInd);   
        GC.fr_max = temp4(mixInd);   
    case 'gcHetero' % hetero (var 1)
        GC.Rm(1:GEN.nGC) = linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        GC.Vth(1:GEN.nGC) = linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        GC.Vreset(1:GEN.nGC) = linspace(-80e-3,-110e-3,GEN.nGC); % AP reset potential (V)
        GC.fr_max(1:GEN.nGC) = linspace(1000,100,GEN.nGC); % maximal frequency
    case 'synapseHetero' % homo (var 2)
        GC.Rm(1:GEN.nGC) = 625e6; %linspace(450e6,800e6,GEN.nGC);      % membrane resistance (Ohm)
        GC.Vth(1:GEN.nGC) = -39e-3;%linspace(-37e-3,-41e-3,GEN.nGC);   % AP threshold (V)
        GC.Vreset(1:GEN.nGC) = -90e-3;    % AP reset potential (V)
        GC.fr_max(1:GEN.nGC) = 500;     % maximal frequency
    otherwise
        error('no paramsetup');
end

GC.spikeShift = 3;   % 3 ms

GC.tau_m = GC.Rm*7e-12;                     % membrane timeconstant (s) rm*cm, cm=7e-12
GC.Vspike = 20e-3;     % AP overshoot (V)
GC.E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and refractory current
GC.E_e = 0;
GC.E_i = 0;

GC.tau_sra = 0.05;      % time for spike-rate adaptation to decay
GC.delta_gsra = 0;      % increase in spike-rate conductance per spike

GC.tau_ref = 0.05;      % time for refractory conductance to decay
GC.delta_gref = 0;      % increase in refractory conductance per spike

GC.Iapp = 0;
GC.tpulse = GEN.dt;
GC.lengthpulse = 0;


%% INPUT - SYNAPSE BETWEEN GRANULE CELL AND PURKINJE CELL
SYN_GCPC.p_base = 0.4;           % baseline probability of release
SYN_GCPC.n_base = 0.7;           % prev.: 0.3
SYN_GCPC.f = 0.2;                % facilitation per spike
SYN_GCPC.tau_f = 0.05;           % timeconstant for facilitation (seconds)
SYN_GCPC.tau_d = 0.5;            % timeconstant for depression (seconds)
g_maxTmp = 2.3e-9/(GEN.nGC/100); % "/(GEN.nGC/100)" to have simlar PC firing greq.

switch CTRL.mode
    case 'homogen'
        SYN_GCPC.tau_rise(1:GEN.nGC) = 1e-3;%linspace(0.2e-3,2e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        SYN_GCPC.tau_decay(1:GEN.nGC) = 35e-3;%linspace(20e-3,50e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        SYN_GCPC.g_max(1:GEN.nGC) = linspace(1.27*g_maxTmp,1.27*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
        %the factor 1.27 is choosen to results in apporximatly the same number
        %of PC spike for homogen and inhomogen models
    case 'homogenR'
        SYN_GCPC.tau_rise(1:GEN.nGC) = 2e-3;%linspace(0.2e-3,2e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        SYN_GCPC.tau_decay(1:GEN.nGC) = 70e-3;%linspace(20e-3,50e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        SYN_GCPC.g_max(1:GEN.nGC) = linspace(1.27*g_maxTmp,1.27*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
    case 'homogenL'
        SYN_GCPC.tau_rise(1:GEN.nGC) = 0.5e-3;%linspace(0.2e-3,2e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        SYN_GCPC.tau_decay(1:GEN.nGC) = 17.5e-3;%linspace(20e-3,50e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        SYN_GCPC.g_max(1:GEN.nGC) = linspace(0.5*g_maxTmp,0.5*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
    case 'heterogen'
        % neu mischen, aber alle gleich
        SYN_GCPC.tau_rise(1:GEN.nGC) = linspace(0.5*1e-3,2*1e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        SYN_GCPC.tau_decay(1:GEN.nGC) = linspace(0.5*35e-3,2*35e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        SYN_GCPC.g_max(1:GEN.nGC) = linspace(2*g_maxTmp,0.5*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
    case 'heterogenMix'
        mixInd = randperm(GEN.nGC);
        temp1 = linspace(0.5*1e-3,2*1e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        temp2 = linspace(0.5*35e-3,2*35e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        temp3 = linspace(2*g_maxTmp,0.5*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
        
        SYN_GCPC.tau_rise = temp1(mixInd);
        SYN_GCPC.tau_decay = temp2(mixInd);
        SYN_GCPC.g_max = temp3(mixInd);
    case 'gcHetero' % homo (var 1)
        SYN_GCPC.tau_rise(1:GEN.nGC) = 1e-3;%linspace(0.2e-3,2e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        SYN_GCPC.tau_decay(1:GEN.nGC) = 35e-3;%linspace(20e-3,50e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        SYN_GCPC.g_max(1:GEN.nGC) = linspace(1.27*g_maxTmp,1.27*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
    case 'synapseHetero' % hetero (var 2)
        SYN_GCPC.tau_rise(1:GEN.nGC) = linspace(0.5*1e-3,2*1e-3,GEN.nGC);    % timeconstant for rise  1 ... 5 ms
        SYN_GCPC.tau_decay(1:GEN.nGC) = linspace(0.5*35e-3,2*35e-3,GEN.nGC);   % timeconstant for decay 20 ... 100 ms
        SYN_GCPC.g_max(1:GEN.nGC) = linspace(2*g_maxTmp,0.5*g_maxTmp,GEN.nGC);% maximal conductance of the synapse (Siemens) %% 80pA/80mV
    otherwise
        error('no paramsetup');
end

% hier kommt der Zeitversatz rein zwischen erster und letzter GC
% GC 100 um 1 ms nach rechts und linear runter zur GC 1 auf 0 ms


%% INPUT - LIF PURKINJE CELL
PC.E = -50e-3;         % Resting membrane potential (V)
PC.Rm = CTRL.PCRm;          % Membrane resistance (Ohm)
PC.tau_m = PC.Rm*100e-12; % Membrane timeconstant (s)
PC.Vth = -45e-3;       % AP threshold (V)
PC.Vreset = -55e-3;    % AP reset potential (V)
PC.Vspike = 20e-3;     % AP overshoot (V)
PC.E_K = -80e-3;       % reversal potential for spike-rate adaptation (SRA) and refractory current
PC.E_e = 0;
PC.E_i = 0;

PC.tau_sra = 0.05;          % time for spike-rate adaptation to decay
PC.delta_gsra = 0;          % increase in spike-rate conductance per spike

PC.tau_ref = 0.05;          % time for refractory conductance to decay
PC.delta_gref = 0;          % increase in refractory conductance per spike

PC.Iapp = 0;           % Applied current (A)
PC.tpulse = GEN.dt;
PC.lengthpulse = 0;