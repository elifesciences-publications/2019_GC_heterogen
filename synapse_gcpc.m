function [g_syn] = synapse_gcpc(GEN,SYN_GCPC,spikeMat)

n = ones(size(spikeMat))*SYN_GCPC.n_base;
p = ones(size(spikeMat))*SYN_GCPC.p_base*(1-SYN_GCPC.p_base);  % previously multiplied with f
r = ones(size(spikeMat));
d = zeros(size(spikeMat));

for i = 1:size(spikeMat,1)
    for j = 2:size(spikeMat,2)
        % Keep track of short term plasticity first:
        p(i,j) = spikeMat(i,j) * SYN_GCPC.f * ( 1 - p(i,j-1) )...
            + p(i,j-1) * exp( -GEN.dt / SYN_GCPC.tau_f );
        n(i,j) = ( 1 - ( 1 - n(i,j-1) ) * exp( -GEN.dt / SYN_GCPC.tau_d ) )...
            - spikeMat(i,j) * p(i,j) * n(i,j-1);
        
        % exponential raise and decay with respect to p an n
        r(i,j) = 1 - spikeMat(i,j) * SYN_GCPC.g_max(i) * p(i,j) * n(i,j) - (1-r(i,j-1)) * exp(-GEN.dt/SYN_GCPC.tau_rise(i)) ;
        d(i,j) = spikeMat(i,j) * SYN_GCPC.g_max(i) * p(i,j) * n(i,j) + d(i,j-1) * exp(-GEN.dt/SYN_GCPC.tau_decay(i));
    end
    %clear j
end
%clear i

% Calculate the new conductance for the synapse:
g_syn = d + r - 1;
        
end