
function tf_res = FiltHilb(data, freqs, freqbloom, srate)

% funciton for filter hilbert transform of data

pnts = size(data,1);
ntrials = size(data,2);

tf_res = zeros(length(freqs),pnts, ntrials);

for f = 1:length(freqs)
    
    tmpdat= reshape(data,1,pnts*ntrials);
 
        
% create filtered noise
frange = [freqs(f)-freqbloom freqs(f)+freqbloom]; % in frequency
order  = round( 3*srate/frange(1) );

% filter kernel
filtkern = fir1(order,frange/(srate/2));

% signals are filtered noise
tf_res(f,:,:) = reshape( filtfilt(filtkern,1,tmpdat), pnts,ntrials);


end