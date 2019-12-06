function [g_syn] = synapse_mfgc(GEN,SYN_MFGC,MF)

n = ones(size(MF.spikeMat))*SYN_MFGC.n_base;
p = ones(size(MF.spikeMat))*SYN_MFGC.p_base*(1-SYN_MFGC.p_base); % previously multiplied with f
r = ones(size(MF.spikeMat));
d = zeros(size(MF.spikeMat));


for i = 1:GEN.nMF
    for j = 2:GEN.nSteps
        % Keep track of short term plasticity first:
        p(i,j) = MF.spikeMat(i,j) * SYN_MFGC.f * ( 1 - p(i,j-1) )...
            + p(i,j-1) * exp( -GEN.dt / SYN_MFGC.tau_f );
        n(i,j) = ( 1 - ( 1 - n(i,j-1) ) * exp( -GEN.dt / SYN_MFGC.tau_d ) )...
            - MF.spikeMat(i,j) * p(i,j) * n(i,j-1);
        
        % exponential raise and decay with respect to p and n
        r(i,j) = 1 - MF.spikeMat(i,j) * SYN_MFGC.g_max * p(i,j) * n(i,j)...
            - (1-r(i,j-1)) * exp(-GEN.dt/SYN_MFGC.tau_rise) ;
        d(i,j) = MF.spikeMat(i,j) * SYN_MFGC.g_max * p(i,j) * n(i,j)...
            + d(i,j-1) * exp(-GEN.dt/SYN_MFGC.tau_decay);
    end
    %clear j
end
%clear i

% Calculate the new conductance for the synapse:
g_syn = d + r - 1;

end


%% plotting option (make function standalone first, comment function... and ...end)
% figure(1);
% 
% MFnumber = 63;
% 
% set(gcf,...
%     'NumberTitle','off',...
%     'Name','LEARNING',...
%     'toolbar','figure',...
%     'units','normalized');
% 
% PLOT.ax(11) = subplot(5,1,1); plot(GEN.tVecSim,n(MFnumber,:));
% title('n');
% xlim([0,GEN.tSim]);
% grid on;
% 
% PLOT.ax(12) = subplot(5,1,2); plot(GEN.tVecSim,p(MFnumber,:));
% title('p');
% xlim([0,GEN.tSim]);
% grid on;
% 
% PLOT.ax(13) = subplot(5,1,3); plot(GEN.tVecSim,r(MFnumber,:));
% title('raise');
% xlim([0,GEN.tSim]);
% grid on;
% 
% PLOT.ax(14) = subplot(5,1,4); plot(GEN.tVecSim,d(MFnumber,:));
% title('decay');
% xlim([0,GEN.tSim]);
% grid on;
% 
% PLOT.ax(15) = subplot(5,1,5); plot(GEN.tVecSim,g_syn(MFnumber,:));
% title('g\_syn');
% xlabel('time [s]');
% xlim([0,GEN.tSim]);
% grid on;
% 
% linkaxes(PLOT.ax,'x');