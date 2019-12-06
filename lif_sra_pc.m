function [V,spikes,gsra,gref] = lif_sra_pc(GEN,PC)

V = zeros(1,GEN.nSteps);        % Vector for membrane potential
spikes = false(1,GEN.nSteps);   % vector to record spike times
gsra = zeros(1,GEN.nSteps);     % spike-rate adaptaion conductance
gref = zeros(1,GEN.nSteps);     % spike refractory conductance
I = zeros(size(GEN.tVecSim));         % Vector for currents

V(:,1) = PC.E;      % begin simulation of every GC with membrane potential at leak level
gsra(:,1) = 0;      % initially no spike-rate adaptation
gref(:,1) = 0;      % initially no refractory current
spike = logical(false);

% apply current during a pulse
for i = PC.tpulse/GEN.dt : (PC.tpulse+PC.lengthpulse)/GEN.dt
    I(round(i)) = PC.Iapp;
end
%clear i

for i = 2:GEN.nSteps;
    % Adaptation and refractory current:
    % Next line solves GC.tau_sra*(d gsra/GEN.dt) = -gsra
    gsra(i) = gsra(i-1) * (1 - GEN.dt/PC.tau_sra);
    % Next line solves GC.tau_ref*(d gref/GEN.dt) = -gref
    gref(i) = gref(i-1) * (1 - GEN.dt/PC.tau_ref);
    
    %if (V(i-1) == PC.Vspike)      % Last datapoint was a spike?
    if spike                       % Last datapoint was a spike?
        V(i) = PC.Vreset;         % reset after spike
        spike = false;                  % clear spikeFlag
    else                            % otherwise integrate
        V(i) = V(i-1)...                            % Last voltage...
            +GEN.dt/PC.tau_m*...
            (PC.E-V(i-1)...                           % +membrane timeconstant relaxation to resting membrane
            +PC.Rm*I(i-1)...                            % +GC.Rm * currentinjection
            -(V(i-1)-PC.E_K)*PC.Rm*gsra(i-1)...     % -Rate adaptation conductance * GC.Rm * driving force
            +(V(i-1)-PC.E_K)*PC.Rm*gref(i-1)...     % +Refractory conductance * GC.Rm * driving force
            -(V(i-1)-PC.E_e)*PC.Rm*PC.g_syne(i-1)...% -excitaroty synaptic conductance * GC.Rm * driving force
            );
    end
    
    if V(i) > PC.Vth       % if voltage is greater than threshold spike
        V(i) = PC.Vspike;     % change voltage so we see the spike
        spikes(i) = true;        % record the time of the spike
        gsra(i) = gsra(i-1) + PC.delta_gsra;   % increase the spike-rate-adaptation conductance
        gref(i) = gref(i-1) + PC.delta_gref;   % increase the refractory conductance
        spike = true;               % set spikeFlag
    end
end
%clear i

end