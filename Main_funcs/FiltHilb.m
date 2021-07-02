
%  This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk

%% Function for narrow band filering of EEG.

% call as tf_res = FiltHilb(data, freqs, freqbloom, srate)

% INPUTS:
% % Where data1 and data2 2d matricies of single channel EEG data e.g., in samples by
% trails fornat. To compute this for all possible channel pairing simply loop throught  channel combinations calling this fucniton each time.
% Freqs is a vector of frequencies that you want to filter the data accordidng to e.g., 2:40 Hz. The function will loop over frequencies
% Freqbloom is the frequency overlap between consecutive bandwidth e.g., at
% 7hz with a freqbloom of 1.5 means that fitler will include data from
% 5.5-8.5Hz. Freqbloom will determine the frequency precision of the data
% Strate is the sampling frequency of the data.

% OUTPUT:
% frequencies by samples by trials matrix of filtered EEG data.


%%

function tf_res = FiltHilb(data, freqs, freqbloom, srate)

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