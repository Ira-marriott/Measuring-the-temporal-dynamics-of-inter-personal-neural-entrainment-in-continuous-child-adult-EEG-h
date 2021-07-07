%% Example for running time varying spectral granger causality


% This section of the tutorial relies on the EEGlab structure
% EEGlab must be downloaded first https://sccn.ucsd.edu/eeglab/download.php

% also reliant on functioniliy from the MVGC toolbox
% download here https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html



% Load in sample dual EEG data provided or load in your own data
load infEEG.mat;
% downsample EEG data
EEG = pop_resample(EEG,128);EEG_infant=EEG;


load mumEEG.mat;
% downsample EEG data
EEG = pop_resample(EEG,128); EEG_mum = EEG;


% inintialise variables 
winsize = 300;
morder = 5;
srate = EEG.srate;
freqs = 2:40;


% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 1;
chan2 = 10;


% get infant and mum data from specified channels and ploynomial detrend data
data1 = detrend(squeeze(EEG_infant.data(chan1,:,:)));
data2 = detrend(squeeze(EEG_mum.data(chan2,:,:)));

% next it is a good idea to check if any trials are not stationary
p=0.05;

[kstat_d1,~] = mvgc_kpss(data1',p);
pval = kstat_d1'>p; 
x = find(pval==1);

% and remove any trials that arent stationary
data1(:,x)=[];
data2(:,x)=[];

[kstat_d2,~] = mvgc_kpss(data2',p);
pval2 = kstat_d2'>p; 
y = find(pval2==1);

% and remove any trials that arent stationary
data2(:,y)=[];
data1(:,y)=[];


% main function for this section is GPDC_tf3 whihc computes time varying granger causality 
% for this function data need to be in samples x trials format - see function for further details

tfg = GC_tf3d(data1 ,data2, freqs, srate, morder,winsize);


%% Plotting

figure
subplot(211);contourf(EEG.times, freqs, squeeze(tfg(1,:,:)), 40, 'linecolor','non'); colorbar
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('Spectral GC y->x')
xlim([-400 800])


subplot(212);contourf(EEG.times, freqs, squeeze(tfg(2,:,:)), 40, 'linecolor','non'); cl = colorbar;
cl.Label.String = 'GC estimate';
cl.FontSize = 30;
set(gca, 'fontsize', 30)
xlabel('Time (ms)'),
title('Spectral GC x->y')
cl.Label.String = 'GC estimate';
xlim([-400 800])
