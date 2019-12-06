function [] = gcc()
% Matlab code for the simulations in the publication:
% Gradients in the mammalian cerebellar cortex enable Fourier-like transformation and improve storing capacity 

clc;
disp('>GCC<');
% version
ver = '2019';
disp(['Version ' ver]);

simTime = tic;  % time of simulation

% parameters 
CTRL.opti = 'LEARN';        % preOptiAlgo
CTRL.nOpti = 5;             % number of optimizations of the optiAlgo
CTRL.psTol = 1E-2;
CTRL.tauRossum = 30E-3;     
CTRL.nGC = 100;             % number of granule cells

CTRL.spikeShift = 0;
CTRL.spikeMatMix = 0;
CTRL.PCRm = 15e6;   
CTRL.mfpergc = 2;
CTRL.n = 100;              % number of runs
CTRL.offsetN = 0;
CTRL.dt = 1E-3;

CTRL.show = 0;  % show some results during runtime | 0,1
CTRL.par = 1;   % parallel computing | 0,1

% mode switch 
switch 2
    case 1,  CTRL.mode = 'homogen';
    case 2,  CTRL.mode = 'heterogen';
    case 3,  CTRL.mode = 'heterogenMix';
    case 4,  CTRL.mode = 'gcHetero';        % var 1
    case 5,  CTRL.mode = 'synapseHetero';   % var 2
    case 6,  CTRL.mode = 'homogenR';
    case 7,  CTRL.mode = 'homogenL';
end

nGC = CTRL.nGC;

rng(2019);
%CTRL.rng = rng;

for forN = 1:numel(nGC)
    CTRL.nGC = nGC(forN); 

    tAlgo = tic;    % time of one run
    n = CTRL.n;

    dateandtime = clock;
    dateandtime(6) = round(dateandtime(6));

    for runthrough = (1 : n) + CTRL.offsetN
        disp('########################');
        disp(['Starting ', num2str(runthrough),'. run-through with ',CTRL.mode, ' | spikeShift: ' num2str(CTRL.spikeShift) ' settings...']);
        disp(['nGC = ' num2str(nGC(forN))]);
        
        % #############
        [mm, rr, tt, rrDuringLearn,PCsipkes] = main(runthrough,CTRL);
        % #############
        
        res{1}(runthrough,:) = mm;
        res{2}(runthrough,:) = rr;
        res{3}(runthrough,:) = tt;
        res{4}(runthrough,:) = rrDuringLearn;
        res{5}(runthrough,:) = PCsipkes;
    end

    close all
    tAlgo = toc(tAlgo)
    disp('Measurement finished.')
    
    CTRL.tAlgo = tAlgo;
    data.mm = res{1};
    data.rr = res{2};
    data.tt = res{3};
    data.rrDuringLearn = res{4};
    data.PCspikes = res{5};
   
    data.name = ver;
    data.date = datestr(clock);
    data.info = CTRL;
    dataSAVE{forN} = data;
    clear data;
    clear res;
end

fileName = 'nGC100hetero';
save(['.\results\dataSAVE_' fileName '.mat'],'dataSAVE');
        
toc(simTime)

end


 