

%% Example for running time varying granger causality using funcitons from MVGC toolbox

% lets set instead supply an example infant adults dual eeg dataset
load sampleEEGdata.mat;

% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 8;
chan2 = 16;

% set window size and model order for GC prediction.
winsize = 300;
morder = 5;
srate = EEG.srate;

% here we use data contained within the eeglab structure but function will
% work for any data inputted in sampls x trails format

data1 = squeeze(EEG.data(chan1,:,:));
data2 = squeeze(EEG.data(chan2,:,:));


tv_gc = TdGc_3d(data1, data2, winsize, morder, srate);


%% Example for obtaining time series of instantaneuous frequency of data 


load sampleEEGdata.mat;

% for example to look at variance over time in 6-9 activity 

frange = [6 9];
freqs = mean(frange);
freqbloom = 1.5;

chan1 = 8;

data1 = squeeze(EEG.data(chan1,:,:));
srate = EEG.srate;

hilsig = squeeze(FiltHilb(data1,freqs, freqbloom, srate));

freqsld = InstFreq(hilsig, srate);

figure; plot(freqsld(:,1))
