
%% Example for running time-frequency power correaltions (spearmans)


load sampleEEGdata.mat;

% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 1;
chan2 = 16;

% set window size and model order for GC prediction.
% winsize = 300;
srate = EEG.srate;
freqbloom = 1.5;
freqs = 2:40;

data1 = squeeze(EEG.data(chan1,:,:));
data2 = squeeze(EEG.data(chan2,:,:));

tf_res(1,:,:,:) = FiltHilb(data1, freqs, freqbloom, srate);
tf_res(2,:,:,:) = FiltHilb(data2, freqs, freqbloom, srate);


bt = dsearchn(EEG.times',[-700 -400]');
pow = abs(hilbert(tf_res)).^2;

baseline_power1 = mean(pow(1,:,bt(1):bt(2),:),3);  dbconverted1 = 10*log10( bsxfun(@rdivide,pow(1,:,:,:),baseline_power1));
baseline_power2 = mean(pow(2,:,bt(1):bt(2),:),3);  dbconverted2 = 10*log10( bsxfun(@rdivide,pow(2,:,:,:),baseline_power2));

a = squeeze(dbconverted1(1,:,:,:));
b = squeeze(dbconverted2(1,:,:,:));

% main funciton for this section 
% calculates single trial tf power correlations
corr_ts = tfPow_corr(a,b);



% some plotting
f = figure(1);
set(gcf,'color','w');

subplot(1, 2, 1)
contourf(EEG.times, freqs, squeeze(mean(dbconverted1(1,:,:,:),4)) ,40,'linecolor','non');
ylabel('Frequency (hz)')
set(gca, 'fontsize', 30, 'clim', [-8 8])
colorbar;
title('Signal x power')
xlim([-300 600])


subplot(1, 2, 2)
contourf(EEG.times, freqs, squeeze(mean(dbconverted2(1,:,:,:),4)) ,40,'linecolor','non');
set(gca, 'fontsize', 30, 'clim', [-8 8])
dl = colorbar;
title('Signal y power')
xlim([-300 600])
dl.Label.String = 'Power(db)';


figure
contourf(EEG.times , freqs, squeeze(mean(corr_ts,3)),40,'linecolor','non');

title('PC at original sampling rate')
xlabel('Time (ms)'), ylabel('Frequency (hz)')
set(gca, 'clim', [-0.2 0.4])
xlim([-300 600])
colorbar;
set(gca, 'fontsize', 30)



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


EEG = pop_resample(EEG,128);

% here we use data contained within the eeglab structure but function will
% work for any data inputted in sampls x trails format

data1 = squeeze(EEG.data(chan1,:,1:10));
data2 = squeeze(EEG.data(chan2,:,1:10));

% main function for this section is tv_gc whihc computes time varying
% granger causality 

% this can take a while to run
tv_gc = TdGc_3d(data1, data2, winsize, morder, srate);

figure(1)
% subplot(234);
plot(times2save,tv_gc(2,:),'b','linew',2)
hold on
plot(times2save,tv_gc(1,:),'g','linew',2)
legend({'x->y','y->x' })
title('Time Domain GC')
set(gca, 'YAxisLocation', 'right', 'fontsize',30)
xlim([-400 800])
hold on
xlabel('Time (ms)')
ylim([0 0.2])

%
%% Example for running time varying spectral granger causality

% lets set instead supply an example infant adults dual eeg dataset
load sampleEEGdata.mat;

% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 1;
chan2 = 10;

% set window size and model order for GC prediction.
winsize = 300;
morder = 5;
srate = EEG.srate;
freqs = 2:40;


t2save = -400:200:800;
tmidx = dsearchn(EEG.times',t2save');

EEG = pop_resample(EEG,128);

% here we use data contained within the eeglab structure but function will
% work for any data inputted in sampls x trails format

data(1,:,:) = squeeze(EEG.data(chan1,:,:));
data(2,:,:) = squeeze(EEG.data(chan2,:,:));

% main function for this section is GPDC_tf3 whihc computes time varying
% granger causality 
tfg = GC_tf3d(data(1,:,:) ,data(2,:,:), EEG.times,t2save, freqs, srate, morder,winsize);

figure
subplot(211);contourf(t2save, freqs, squeeze(tfg(1,2,:,:)), 40, 'linecolor','non'); colorbar
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('Spectral GC y->x')


subplot(212);contourf(t2save, freqs, squeeze(tfg(1,1,:,:)), 40, 'linecolor','non'); cl = colorbar;
cl.Label.String = 'GC estimate';
cl.FontSize = 30;
set(gca, 'fontsize', 30)
xlabel('Time (ms)'),
title('Spectral GC x->y')
cl.Label.String = 'GC estimate';


%% Example for obtaining time series of instantaneuous frequency of data 

load sampleEEGdata.mat;

% for example to look at variance over time in 6-9 activity 

center_freq = 7.5;
freqbloom = 1.5;
chan1 = 8;
srate = EEG.srate;

data1 = squeeze(EEG.data(chan1,:,:));

filtsig = FiltHilb(data1, center_freq, freqbloom, srate);

hilsig = squeeze(hilbert(filtsig));

plot(squeeze(real(hilsig(:,1))))

freqsld = InstFreq(hilsig, srate);

figure; plot(freqsld(:,1))







