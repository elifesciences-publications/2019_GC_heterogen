function [ee1,ee2] = getTemporalError(GEN,LEARN,inputspks)
inputspks = logical(inputspks);
%targetspks = logical(LEARN.trgtWin);
% Find indexes of the spikes in the input and target spiketrains:
%LEARN.trgtpos = find(targetspks);
LEARN.pcpos = find(inputspks);

pausePosMiddle = round((GEN.pausePosBegin + GEN.pausePosEnd)/2);

if isempty(find(LEARN.pcpos < pausePosMiddle)) == 1
    ee1 = GEN.tWin(1);
else
    ee1 = GEN.dt*LEARN.pcpos(max(find(LEARN.pcpos < pausePosMiddle)));
end

if isempty(find(LEARN.pcpos >= pausePosMiddle)) == 1
    ee2 = GEN.tWin(2);
else
    ee2 = GEN.dt*LEARN.pcpos(min(find(LEARN.pcpos >= pausePosMiddle)));
end

end
