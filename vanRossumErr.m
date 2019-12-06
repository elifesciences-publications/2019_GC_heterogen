function [meanpcerror] = vanRossumErr(spikeMatWin,dd,trgtcon)
pccon = conv2(double(spikeMatWin),dd,'full');
diffcon = (pccon - trgtcon).^2;
meanpcerror = sum(diffcon);
end


