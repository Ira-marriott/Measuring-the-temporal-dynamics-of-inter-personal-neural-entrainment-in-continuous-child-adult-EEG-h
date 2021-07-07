
%% Example for running time-frequency power correaltions (spearmans)

% load sampleEEGdata.mat;
load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;


% define channels manually or user can loop over different channel pairings
% although this is takes a while
chan1 = 16;
chan2 = 16;

% set window size and model order for GC prediction.
% winsize = 300;
srate = EEG.srate;
freqbloom = 1.5;
freqs = 2:40;

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
set(gca, 'fontsize', 30, 'clim', [-5 1])
colorbar;
title('Signal x power')
xlim([-300 800])


subplot(1, 2, 2)
contourf(EEG.times, freqs, squeeze(mean(dbconverted2(1,:,:,:),4)) ,40,'linecolor','non');
set(gca, 'fontsize', 30, 'clim', [-5 1])
dl = colorbar;
title('Signal y power')
xlim([-300 800])
dl.Label.String = 'Power(db)';


figure
contourf(EEG.times , freqs, squeeze(mean(corr_ts,3)),40,'linecolor','non');
title('PC at original sampling rate')
xlabel('Time (ms)'), ylabel('Frequency (hz)')
set(gca, 'clim', [-0.2 0.4])
xlim([-300 800])
colorbar;
set(gca, 'fontsize', 30)


%% Example for running time varying granger causality using funcitons from MVGC toolbox

% lets set instead supply an example infant adults dual eeg dataset
load infEEG.mat;
% downsample EEG data
EEG = pop_resample(EEG,128);EEG_infant=EEG;


load mumEEG.mat;
% downsample EEG data
EEG = pop_resample(EEG,128); EEG_mum = EEG;

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

% ploynomial detrend data
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


% main function for this section is tv_gc which computes time varying
% granger causality 
% this can take a while to run for data with a large number of trials

tv_gc = GC_td3d(data1, data2, winsize, morder, srate);


figure
% subplot(234);
plot(EEG.times,squeeze(mean(tv_gc(2,:,:),3)),'b','linew',2)
hold on
plot(EEG.times,squeeze(mean(tv_gc(1,:,:),3)),'g','linew',2)
legend({'x->y','y->x' })
title('Time Domain GC')
set(gca, 'YAxisLocation', 'right', 'fontsize',30)
xlim([-400 800])
hold on
xlabel('Time (ms)')
% ylim([0 0.2])


%
%% Example for running time varying spectral granger causality

% lets set instead supply an example infant adults dual eeg dataset
load infEEG.mat;
% downsample EEG data
EEG = pop_resample(EEG,128);EEG_infant=EEG;


load mumEEG.mat;
% downsample EEG data
EEG = pop_resample(EEG,128); EEG_mum = EEG;

% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 1;
chan2 = 10;

% set window size and model order for GC prediction.
winsize = 300;
morder = 5;
srate = EEG.srate;
freqs = 2:40;


% here we use data contained within the eeglab structure but function will
% work for any data inputted in sampls x trails format

% ploynomial detrend data
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

% reshape data to fit into the structure
data(1,:,:) = data1;
data(2,:,:)=data2;

% main function for this section is GPDC_tf3 whihc computes time varying
% granger causality 
tfg = GC_tf3d(data(1,:,:) ,data(2,:,:), freqs, srate, morder,winsize);

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



%% Example for obtaining time series of instantaneuous frequency of data 

load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;

% for example to look at variance over time in 6-9 activity 

center_freq = 7.5;
freqbloom = 1.5;
chan1 = 8;
srate = EEG.srate;

data1 = squeeze(EEG_infant.data(chan1,:,:));
data2 = squeeze(EEG_mum.data(chan2,:,:));

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


% some plotting
figure; subplot(211)

plot(EEG.times(1:end-1), freqsld(:,2), 'linew',2)
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Instantaneous Frequency - infant data')
ylim([0 12])
xlim([-400 800])

subplot(212)
plot(EEG.times(1:end-1), freqsld2(:,2), 'linew',2)
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title('Instantaneous Frequency - mum data')
ylim([0 12])
xlim([-400 800])



%% Example for obtaining plv over time/ trials 

load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;

% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 1;
chan2 = 16;

% set window size and model order for GC prediction.
srate = EEG.srate;
freqs = 2:40;
winsize= 100;

data1 = squeeze(EEG_infant.data(chan1,:,:));
data2 = squeeze(EEG_mum.data(chan2,:,:));


% main function for this section compute plv over time and/ or trials
[plv_time, plv_trials] = PLV(data1, data2, freqs, srate, winsize);


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


%% Example for obtaining pte

load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;

% initialise variables
order = 15;
twin = 200;
freqz=2:40;
ntrials = EEG.trials;
twin_pnts = round(twin/(1000/EEG.srate));
srate = EEG.srate;
freqbloom = 1.5;
pnts = EEG.pnts;
times2save = -700:50:1000;
times2saveidx = dsearchn(EEG.times',times2save');


% define channels manually or user can loop over different channel pairings
% although this is computationally very slow 
chan1 = 1;
chan2 = 16;

data1 = squeeze(EEG_infant.data(chan1,:,:));
data2 = squeeze(EEG_mum.data(chan2,:,:));

% filter the data 
filtsig = FiltHilb(data1, freqz, freqbloom, srate);
filtsig2 = FiltHilb(data2, freqz, freqbloom, srate);

% initialise PTE output
pte = zeros(2, length(freqz), length(times2save), ntrials);

for fi = 1:length(freqz)

for ti = 1:ntrials

for i = 1:length(times2save)
    
    tidx = times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2); 
    X = squeeze(cat(1,filtsig(fi,tidx,ti), filtsig2(fi,tidx,ti))) ;
    
    
% X = squeeze(EEG.data([chan1 chan2],times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),ti));

% main fucntion for this section compute phase transfer entropy on narrow
% band filtered data

[dPTE, PTE] = PhaseTE_MF(X', [], []);

pte(1,fi,i,ti) = PTE(1,2);
pte(2,fi,i,ti) = PTE(2,1);

end
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
plot( times2save, squeeze(mean(pte(1,8,:,:),4)),'b', 'linew',2);
hold on 
plot( times2save, squeeze(mean(pte(2,8,:,:),4)),'g', 'linew',2);
set(gca, 'fontsize',30)
title('Phase Transfer Entropy')
xlabel('Time (ms)')
xlim([-400 1000])
legend({'x->y', 'y->x'})
ylabel('PTE estimate')



