function [V,spikes,gsra,gref] = lif_sra_gc(GEN,GC)

V = zeros(GEN.nGC,GEN.nSteps);        % Vector for membrane potential
spikes = false(GEN.nGC,GEN.nSteps);   % vector to record spike times
gsra = zeros(GEN.nGC,GEN.nSteps);     % spike-rate adaptaion conductance
gref = zeros(GEN.nGC,GEN.nSteps);     % spike refractory conductance
I = zeros(size(GEN.tVecSim));         % Vector for currents

V(:,1) = GC.E;      % begin simulation of every GC with membrane potential at leak level
gsra(:,1) = 0;      % initially no spike-rate adaptation
gref(:,1) = 0;      % initially no refractory current

% applie current during a pulse
for i = GC.tpulse/GEN.dt : (GC.tpulse+GC.lengthpulse)/GEN.dt
    I(round(i)) = GC.Iapp;
end
%clear i

for i = 1:GEN.nGC
    for j = 2:GEN.nSteps;
        % Adaptation and refractory current:
        % Next line solves GC.tau_sra*(d gsra/GEN.dt) = -gsra
        gsra(i,j) = gsra(i,j-1) * (1 - GEN.dt/GC.tau_sra);
        % Next line solves GC.tau_ref*(d gref/GEN.dt) = -gref
        gref(i,j) = gref(i,j-1) * (1 - GEN.dt/GC.tau_ref);
        
        if (V(i,j-1) == GC.Vspike)      % Last datapoint was a spike?
            V(i,j) = GC.Vreset(i);         % reset after spike
        else                            % otherwise integrate
            V(i,j) = V(i,j-1)...                            % Last voltage...
                +GEN.dt/GC.tau_m(i)*...
                (GC.E-V(i,j-1)...                           % +membrane timeconstant relaxation to resting membrane
                +GC.Rm(i)*I(j-1)...                            % +GC.Rm * currentinjection
                -(V(i,j-1)-GC.E_K)*GC.Rm(i)*gsra(i,j-1)...     % -Rate adaptation conductance * GC.Rm * driving force
                +(V(i,j-1)-GC.E_K)*GC.Rm(i)*gref(i,j-1)...     % +Refractory conductance * GC.Rm * driving force
                -(V(i,j-1)-GC.E_e)*GC.Rm(i)*GC.g_syne(i,j-1)...% -excitaroty synaptic conductance * GC.Rm * driving force
                -(V(i,j-1)-GC.E_i)*GC.Rm(i)*GC.g_syni(i,j-1)...   % -inhibitory synaptic conductance * GC.Rm * driving force
                );
        end
        
        if V(i,j) > GC.Vth(i)       % if voltage is greater than threshold spike
            V(i,j) = GC.Vspike;     % change voltage so we see the spike
            spikes(i,j) = true;        % record the time of the spike
            gsra(i,j) = gsra(i,j-1) + GC.delta_gsra;   % increase the spike-rate-adaptation conductance
            gref(i,j) = gref(i,j-1) + GC.delta_gref;   % increase the refractory conductance
        end
    end
    %clear j
end
clear i
end

