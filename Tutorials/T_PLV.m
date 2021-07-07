%% Example for obtaining plv over time/ trials 


% This section of the tutorial relies on the EEGlab structure 
% EEGlab must be downloaded first https://sccn.ucsd.edu/eeglab/download.php


% Load in sample dual EEG data provided or load in your own data
load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;


% intialise variables
srate = EEG.srate;
freqs = 2:40;
winsize= 100;


% define channels manually or user can loop over different channel pairings
chan1 = 1;
chan2 = 16;

% get infant and mum data from specified channels
data1 = squeeze(EEG_infant.data(chan1,:,:));
data2 = squeeze(EEG_mum.data(chan2,:,:));


% main function for this section compute plv over time and/ or trials
% for this function data need to be in samples x trials format - see function for further details
[plv_time, plv_trials] = PLV(data1, data2, freqs, srate, winsize);


%% Plotting

figure
contourf(EEG.times, freqs, plv_trials, 40, 'linecolor','non'); colorbar
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('PLV over trials')
xlim([-400 800])


figure
contourf(EEG.times, freqs, squeeze(mean(plv_time,3)), 40, 'linecolor','non'); colorbar
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('PLV within trials - sliding time window')
xlim([-400 800])