function [LEARN] = geterror(GEN,LEARN,inputspks)
inputspks = logical(inputspks);
targetspks = logical(LEARN.trgtWin);
% Find indexes of the spikes in the input and target spiketrains:
LEARN.trgtpos = find(targetspks);
LEARN.pcpos = find(inputspks);

% Find distances: (multiplied by dt for normalization)
LEARN.distanceMat = GEN.dt * (abs(repmat(LEARN.pcpos.',[1,length(LEARN.trgtpos)]) - repmat(LEARN.trgtpos,[length(LEARN.pcpos),1])));

% Determine the peaks of the error function:
LEARN.mindistanceVec = min(LEARN.distanceMat,[],2);

% Distribute the errors to a timevector:
errortrace = zeros(size(inputspks));
errortrace(LEARN.pcpos) = LEARN.mindistanceVec';

% Make a kernel for the learning rule:
imp = zeros(1,int64((LEARN.tau_rise+LEARN.tau_decay)/GEN.dt*10));
imp(2) = 1;
r = ones(size(imp));
d = zeros(size(imp));
for i=2:size(imp,2)
    r(i) = (1-(1-r(i-1))*exp(-GEN.dt/LEARN.tau_rise))-imp(i);
    d(i) = imp(i)+d(i-1)*exp(-GEN.dt/LEARN.tau_decay);
    k(i) = r(i)+d(i)-1;
end
k(1) = [];
k=fliplr(k);
converrortrace = conv(errortrace,k);
converrortrace(1:size(k,2)-1) =[];

LEARN.errorfct = converrortrace;

end
