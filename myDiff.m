function x = myDiff(syn, K)
% K: Konstanten

syn = syn - K.LEARN.fminsearchAdd;

g_syne = syn * K.SYN_GCPC.g_syn;

VV = zeros(1,K.GEN.nSteps);         % vector for membrane potential
spikes = false(1,K.GEN.nSteps);     % vector to record spike times

VV(:,1) = K.PC.E;                   % begin simulation of every GC with membrane potential at leak level
spike = logical(false);

for i = 2:K.GEN.nSteps;
    
    if spike                            % Last datapoint was a spike?
        VV(i) = K.PC.Vreset;            % reset after spike
        spike = false;                  % clear spikeFlag
    else                                % otherwise integrate
        VV(i) = VV(i-1) * (K.k2-K.k3*g_syne(i-1)) + K.k1 + K.k3 * K.I(i-1) + K.k4 * g_syne(i-1);
        
        if VV(i) > K.PC.Vth             % if voltage is greater than threshold spike
            VV(i) = K.PC.Vspike;        % change voltage so we see the spike
            spikes(i) = true;           % record the time of the spike
            spike = true;               % set spikeFlag
        end
    end
end

spikeMatWin = double(spikes(K.GEN.WinSteps));   % cut out the window

%%%%%%%%%%%%%%%%%%%%%% Temporal Error %%%%%%%%%%%%%%%%%
% pcpos = find(spikeMatWin);
% 
% temp = find(pcpos < K.pausePosMiddle,1,'last');     % max
% if isempty(temp) == 1
%     ee1 = K.GEN.tWin(1);
% else
%     ee1 = K.GEN.dt*pcpos(temp);
% end
% 
% temp = find(pcpos >= K.pausePosMiddle,1,'first');   % min
% if isempty(temp) == 1
%     ee2 = K.GEN.tWin(2);
% else
%     ee2 = K.GEN.dt*pcpos(temp);
% end
% 
% temporalError = K.LEARN.temporalErrorScaling * (abs(ee1-K.ee1Trgt)+abs(ee2-K.ee2Trgt));

%%%%%%%%%%%%%%%%%%%%%%%%%% VAN ROSSUM %%%%%%%%%%%%%%%%%
pccon = conv2(spikeMatWin,K.dd,'full'); 
%diffcon = (pccon(1:end-size(K.dd,2)+1) - K.LEARN.trgtcon);
diffcon = (pccon - K.LEARN.trgtcon);
meanpcerror2 = diffcon * diffcon';
      
%%%%%%% Sum of temporal error and ROSSUM %%%%%%%%%%%%%%
x = meanpcerror2;% + temporalError;
end

% VV(i) = VV(i-1)...                                  % Last voltage...
%     +K.GEN.dt/K.PC.tau_m*...
%     (K.PC.E-VV(i-1)...                              % +membrane timeconstant relaxation to resting membrane
%     +K.PC.Rm*K.I(i-1)...                            % +GC.Rm * currentinjection
%     -(VV(i-1)-K.PC.E_e)*K.PC.Rm*g_syne(i-1)...      % -excitaroty synaptic conductance * GC.Rm * driving force
%     );

% a = K.GEN.dt;
% b = K.PC.tau_m;
% c = K.PC.E;
% d = K.PC.Rm;
% e = K.PC.E_e;

% VV(i) = VV(i-1) + a / b * (c - VV(i-1) + d * K.I(i-1) - (VV(i-1) - e) * d * g_syne(i-1)); % 5 Multiplikationen

% k1 = (a*c)/b;
% k2 = 1-(a/b);
% k3 = (a*d)/b;
% k4 = (a*d)/(e*b);
% 
% VV(i) = VV(i-1) * (k2-k3*g_syne(i-1)) + k1 + k3 * K.I(i-1) + k4 * g_syne(i-1);            % 4 Multiplikationen       

