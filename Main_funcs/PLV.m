

%% This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk



%% function for time varying granger


% Call as [plv_time, plv_trials] = PLV(data1, data2, freqs, srate, winsize)

% INPUTS:
% Where data1 and data2 2d matricies of single channel EEG data e.g., in samples by
% trails format. To compute this for all possible channel pairing simply loop throught
% channel combinations calling this fucniton each time.
% Freqs is a vector of frequencies that you want to filter the data accordidng to e.g., 2:40 Hz. The function will loop over frequencies
% Freqbloom is the frequency overlap between consecutive bandwidth e.g., at
% 7hz with a freqbloom of 1.5 means that fitler will include data from
% 5.5-8.5Hz. Freqbloom will determine the frequency precision of the data
% Winsize is the size of the sliding window used to compute the GC estimate - given in ms.
% Strate is the sampling frequency of the data.


% OUTPUT:
% PLV_time is matrix of frequency x time x trials of plv estimates compute
% over single trials



%% 

function [plv_time, plv_trials] = PLV(data1, data2, freqs, srate, winsize)

pnts = size(data1,1);
ntrials = size(data1,2);

freqoverlap=1.5;

tf_res(1,:,:,:) = FiltHilb(data1, freqs, freqoverlap, srate);
tf_res(2,:,:,:) = FiltHilb(data2, freqs, freqoverlap, srate);

ang_dat = angle(hilbert(permute(tf_res,[3,1,2,4])));
ang_data = permute(ang_dat,[2,3,1,4]);


angles1 = squeeze(ang_data(1,:,:,:));
angles2 = squeeze(ang_data(2,:,:,:));


% plv over trails 
plv_trials =  abs(mean(exp(1i* (angles1-angles2) ),3));


% plv over time within trails 

% initialize
plv_time = zeros(length(freqs),pnts,ntrials);

% loop over time points
for ti=1:pnts
    
    % figure out time points, considering edges
    tidx = max(1,ti-winsize) : min(pnts,ti+winsize);
        
    
    % compute synchronization for this center time points
    plv_time(:,ti,:) = abs(mean(exp(1i*( angles1(:,tidx,:)-angles2(:,tidx,:) )), 2));
end


