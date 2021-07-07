%% Example for running time-frequency power correaltions (spearmans)


% This section of the tutorial relies on the EEGlab structure 
% EEGlab must be downloaded first https://sccn.ucsd.edu/eeglab/download.php


% Load in sample dual EEG data provided or load in your own data 
load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;


% define channels manually or user can loop over different channel pairings
% although this is takes a while
chan1 = 16;
chan2 = 16;


% initilaise variabels that will be needed for later
srate = EEG.srate;
freqbloom = 1.5;
freqs = 2:40;

% get infant and mum data from specified channels
data1 = squeeze(EEG_infant.data(chan1,:,:));
data2 = squeeze(EEG_mum.data(chan2,:,:));

% filter the data 
tf_res=[];
tf_res(1,:,:,:) = FiltHilb(data1, freqs, freqbloom, srate);
tf_res(2,:,:,:) = FiltHilb(data2, freqs, freqbloom, srate);

% define baseline times for power normalisation
bt = dsearchn(EEG.times',[-700 -400]');

% permute data as hilbert function computes hilbert over first dimension
pow = abs(hilbert(tf_res)).^2;
pow2 = abs(hilbert(permute(tf_res,[3,1,2,4]))).^2;

% re permtue data back into original format
pow3 = permute(pow2,[2,3,1,4]);

% decibel conversion/ baseline normalisation of power
baseline_power1 = mean(pow3(1,:,bt(1):bt(2),:),3);  dbconverted1 = 10*log10( bsxfun(@rdivide,pow3(1,:,:,:),baseline_power1));
baseline_power2 = mean(pow3(2,:,bt(1):bt(2),:),3);  dbconverted2 = 10*log10( bsxfun(@rdivide,pow3(2,:,:,:),baseline_power2));


inf_pow_dat = squeeze(dbconverted1(1,:,:,:));
mum_pow_dat = squeeze(dbconverted2(1,:,:,:));

% main funciton for this section calculates single trial tf power correlations
% input data here are frequency x sampels x trails matrix of power
corr_ts = tfPow_corr(inf_pow_dat,mum_pow_dat);


%% Plotting

figure(1);
set(gcf,'color','w');

subplot(1, 2, 1)
contourf(EEG.times, freqs, squeeze(mean(dbconverted1(1,:,:,:),4)) ,40,'linecolor','non');
ylabel('Frequency (hz)')
set(gca, 'fontsize', 30, 'clim', [-5 1])
colorbar;
title('infant brain power')
xlim([-300 800])


subplot(1, 2, 2)
contourf(EEG.times, freqs, squeeze(mean(dbconverted2(1,:,:,:),4)) ,40,'linecolor','non');
set(gca, 'fontsize', 30, 'clim', [-5 1])
dl = colorbar;
title('mum brain power')
xlim([-300 800])
dl.Label.String = 'Power(db)';


figure
set(gcf,'color','w');

contourf(EEG.times , freqs, squeeze(mean(corr_ts,3)),40,'linecolor','non');
title('PC at original sampling rate')
xlabel('Time (ms)'), ylabel('Frequency (hz)')
set(gca, 'clim', [-0.2 0.4])
xlim([-300 800])
colorbar;
set(gca, 'fontsize', 30)


