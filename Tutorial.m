
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
tv_gc = GC_td3d(data1, data2, winsize, morder, srate);

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

freqsld = InstFreq(hilsig, srate);


figure; plot(freqsld(:,1), 'linew',2)
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Instantaneous Frequency')




%% Example for obtaining plv over time/ trials 

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
ntrials = EEG.trials;
pnts = EEG.pnts;
winsize= 100;

data1 = squeeze(EEG.data(chan1,:,:));
data2 = squeeze(EEG.data(chan2,:,:));

% main function for this section compute plv over time and/ or trials
[plv_time, plv_trials] = PLV(data1, data2, freqs, srate, winsize);


figure
contourf(EEG.times, freqs, plv_trials, 40, 'linecolor','non'); colorbar
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('PLV over trials')



figure
contourf(EEG.times, freqs, squeeze(mean(plv_time,3)), 40, 'linecolor','non'); colorbar
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('PLV within trials - sliding time window')


%% Example for obtaining pte

load sampleEEGdata.mat

order = 15;
twin = 200;
freqz=2:40;
ntrials = EEG.trials;

times2save = -700:50:1000;
times2saveidx = dsearchn(EEG.times',times2save');
twin_pnts = round(twin/(1000/EEG.srate));


% narrowband filter data
EEG2 = pop_eegfiltnew(EEG,8,12);


pte = zeros(2, length(times2save));

chan1 = 1;
chan2 = 2;

for ti = 1:ntrials

for i = 1:length(times2save)
    
X = squeeze(EEG2.data([chan1 chan2],times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),ti));

% main fucntion for this section compute phase transfer entropy on narrow
% band filtered data

[dPTE, PTE] = PhaseTE_MF(X', [], []);

pte(1,i,ti) = PTE(1,2);
pte(2,i,ti) = PTE(2,1);

end
end

figure(3)
subplot(211)
plot(EEG.times,squeeze(mean(EEG.data(1,:,:),3)),'b','linew',2 )
hold on
plot(EEG.times,squeeze(mean(EEG.data(2,:,:),3))-3,'g','linew',2 )
xlim([-400 1000])
set(gca, 'fontsize',30,'ytick',[])
legend({'x', 'y'})
title('Raw amplitude over time')



subplot(212)
plot( times2save, squeeze(mean(pte(1,:,:),3)),'b', 'linew',2);
hold on 
plot( times2save, squeeze(mean(pte(2,:,:),3)),'g', 'linew',2);
set(gca, 'fontsize',30)
title('Phase Transfer Entropy')
xlabel('Time (ms)')
xlim([-400 1000])
legend({'x->y', 'y->x'})
ylabel('PTE estimate')



