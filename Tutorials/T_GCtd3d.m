%% Example for running time varying granger causality using funcitons from MVGC toolbox


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

% initialise variables
winsize = 300;
morder = 5;
srate = EEG.srate;

% define channels manually or user can loop over different channel pairings
% although this is computationally very slow
chan1 = 8;
chan2 = 16;


% get infant and mum data from specified channels and ploynomial detrend data
data1 = detrend(squeeze(EEG_infant.data(chan1,:,:)));
data2 = detrend(squeeze(EEG_mum.data(chan2,:,:)));


% next it is a good idea to check if any trials are not stationary
p=0.05;

[kstat_d1,~] = mvgc_kpss(data1',p);
pval = kstat_d1'>p;
x = find(pval==1);

% and remove any trials that aren't stationary
data1(:,x)=[];
data2(:,x)=[];

[kstat_d2,~] = mvgc_kpss(data2',p);
pval2 = kstat_d2'>p;
y = find(pval2==1);

% and remove any trials that aren't stationary
data2(:,y)=[];
data1(:,y)=[];


% main function for this section is tv_gc which computes time varying granger causality this can take a while to run for data with a large number of trials
% for this function data need to be in samples x trials format - see function for further details

tv_gc = GC_td3d(data1, data2, winsize, morder, srate);


%% Plotting

figure
plot(EEG.times,squeeze(mean(tv_gc(2,:,:),3)),'b','linew',2)
hold on
plot(EEG.times,squeeze(mean(tv_gc(1,:,:),3)),'g','linew',2)
legend({'x->y','y->x' })
title('Time Domain GC')
set(gca, 'YAxisLocation', 'right', 'fontsize',30)
xlim([-400 800])
hold on
xlabel('Time (ms)')
