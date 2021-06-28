

% This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk

%% 


%% 

function [plv_time, plv_trials] = PLV(data1, data2, freqs, srate, winsize)

pnts = size(data1,1);
ntrials = size(data1,2);

freqoverlap=1.5;

tf_res(1,:,:,:) = FiltHilb(data1, freqs, freqoverlap, srate);
tf_res(2,:,:,:) = FiltHilb(data2, freqs, freqoverlap, srate);

angles1 = squeeze(angle(hilbert( tf_res(1,:,:,:) )));
angles2 = squeeze(angle(hilbert( tf_res(2,:,:,:) )));

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


