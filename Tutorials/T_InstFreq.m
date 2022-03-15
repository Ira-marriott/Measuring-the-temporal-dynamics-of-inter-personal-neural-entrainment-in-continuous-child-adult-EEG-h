%% Example for obtaining time series of instantaneuous frequency of data 


% This section of the tutorial relies on the EEGlab structure 
% EEGlab must be downloaded first https://sccn.ucsd.edu/eeglab/download.php


% Load in sample dual EEG data provided or load in your own data
load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;

% Initialise variables
center_freq = 7.5;
freqbloom = 1.5;
srate = EEG.srate;
ntrials = EEG.trials;

% define channels manually or user can loop over different channel pairings
chan1 = 1;
chan2 = 10;

% get infant and mum data from specified channels
data1 = double(squeeze(EEG_infant.data(chan1,:,:)));
data2 = double(squeeze(EEG_mum.data(chan2,:,:)));

% filter the data 
filtsig = FiltHilb(data1, center_freq, freqbloom, srate);
filtsig2 = FiltHilb(data2, center_freq, freqbloom, srate);

% obtain hilbert transform of the data
hilsig = hilbert(squeeze(filtsig));
hilsig2 = hilbert(squeeze(filtsig2));

% get instantaneuous frequency of the data for mum and baby
freqsld = InstFreq(hilsig, srate);
freqsld2 = InstFreq(hilsig2, srate);

% median filter result
freqsld = medfilt1(freqsld);
freqsld2 = medfilt1(freqsld2);


%% plotting

figure; 
subplot(211)
plot(EEG.times(1:end-1), freqsld(:,randi(ntrials)), 'linew',2)
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Instantaneous Frequency - infant data')
ylim([0 12])
xlim([-400 800])

subplot(212)
plot(EEG.times(1:end-1), freqsld2(:,randi(ntrials)), 'linew',2)
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Instantaneous Frequency - mum data')
ylim([0 12])
xlim([-400 800])