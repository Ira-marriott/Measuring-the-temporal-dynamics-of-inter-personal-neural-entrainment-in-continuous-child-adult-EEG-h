
%% This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk



%% function for time varying granger

% Call as corr_ts = tfPow_corr(a,b)

% INPUTS:
% Where data1 and data2 3d matrices of single channel EEG power in the
% format frequency x samples x trials 



% OUTPUT:
% Corr_ts is a frequency x samples x trials matrix of correlations between
% single trial power data of data1 and data2



function corr_ts = tfPow_corr(data1,data2)

% assumes a and b and frequencies by time by trails data

freqs = size(data1,1);
pnts =size(data1,2);
ntrials = size(data1,3);


% correlation at original sample rate
corr_ts = zeros(length(freqs),pnts,ntrials);

for i = 1:freqs
    
    for ti=1:pnts
        
        x = corrcoef(data1(i,ti,:), data2(i,ti,:));
        corr_ts(i,ti,:) = x(2,1);
        
    end
    
end
