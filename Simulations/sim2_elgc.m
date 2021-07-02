
%  This code accompnies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk  


% also see - Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.

% and see - Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference. Journal of neuroscience methods, 223, 50-68.



%% section 2 - event locked changes in granger causality


load sampleEEGdata.mat
times = EEG.times;

% simulation parameters
sfreq = 7;
ntrials = 10;
srate = 256;
time  = -1:1/srate:1.5;
gwidth = 0.1;
noiseamp = 0.5; % noise standard deviation

signal = zeros(ntrials, length(time));
signal2 = zeros(ntrials, length(time));


for ti = 1 :ntrials
    
    swave  = (3*randn)*sin(2*pi*sfreq*time);
    gausw  = exp( -4*log(2)*(time).^2 / gwidth^2 );
    signal(ti,:) = swave .* gausw;
    
    
    for i=26:length(signal)
        
        signal2(ti,i) = signal(ti,i-25);
        
    end
    
    
end

% create data and multiple channels plus noise
data1 = signal + noiseamp * randn(ntrials,length(time));
data2 = signal2 + noiseamp * randn(ntrials,length(time));

%% fit to eeglab structure for ease of event locked analysis

clear EEG
clear data

data(1,:,:) = data1(:,1:end-1)';
data(2,:,:) = data2(:,1:end-1)';

EEG          = eeg_emptyset;
EEG.data     = data;
EEG.pnts     = size(data,2);
EEG.srate    = srate;
EEG.nbchan   = size(data,1);
EEG.event    = []; % overweighted events not in EEG.event structure
EEG          = eeg_checkset(EEG);
EEG.trials   = size(EEG.data,3);
EEG.times = times;


times2save = -700:200:1000;
dsearchn(EEG.times',times2save')
%

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

morder = 5;
fs = 256;
twin  =  200;
twin_pnts = round(twin/(1000/fs));

times2save = -800:50:1000;
times2saveidx = dsearchn(EEG.times',times2save');


tv_gc = zeros(2, length(times2save));


for i = 1:length(times2save)
    
    X = EEG.data([1 2],times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),:);
    
        % taken from MVGC toolbox 
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    var_info(info,true);
    
    F = autocov_to_pwcgc(G);
    
    tv_gc(1,i) = F(1,2);
    tv_gc(2,i) = F(2,1);

end




%%

pnts= length(EEG.times);
bt = dsearchn(EEG.times',[-700 -400]');

freqs = linspace(2,40,38);
freqbloom = 1.5;


tf = zeros(2,length(freqs),pnts, ntrials);

for f = 1:length(freqs)
    
    tmp= reshape(data,2,pnts*ntrials);
    
    for i=1:2
        
        tmpdat = tmp(i,:);
        
        % create filtered noise
        frange = [freqs(f)-freqbloom freqs(f)+freqbloom]; % in frequency
        order  = round( 3*srate/frange(1) );
        
        % filter kernel
        filtkern = fir1(order,frange/(srate/2));
        
        % signals are filtered noise
        tf(i,f,:,:) = reshape( filtfilt(filtkern,1,tmpdat), pnts,ntrials);
        
        
    end
    
end


pow1 = abs(hilbert(tf(1,:,:,:))).^2;
pow2 = abs(hilbert(tf(2,:,:,:))).^2;


baseline_power1 = mean(pow1(:,:,bt(1):bt(2),:),3);  dbconverted1 = 10*log10( bsxfun(@rdivide,pow1,baseline_power1));
baseline_power2 = mean(pow2(:,:,bt(1):bt(2),:),3);  dbconverted2 = 10*log10( bsxfun(@rdivide,pow2,baseline_power2));

t2save = -400:200:800;
tmidx = dsearchn(EEG.times',t2save');


d2p = squeeze(dbconverted1);
d2p2 = squeeze(dbconverted2);


twin = 400;
tfg = GPDC_tf3(data(1,:,:) ,data(2,:,:), EEG.times,t2save, freqs, srate, morder,twin);

%% plotting for figure 3


figure(1)
set(gcf,'color','w');

subplot(231)
plot(times,squeeze(mean(data1(:,1:end-1),1)),'color','b','linew',2)
hold on
plot(times,squeeze(mean(data2(:,1:end-1),1))-2,'g','linew',2)
% xlim([-0.5 .8])
legend({'x', 'y'})
xlim([-400 800])
set(gca, 'ytick',[], 'fontsize',30)
title('Raw Amplitude')
% axis square


figure(1)
subplot(234);
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

figure(1);
subplot(232);contourf(EEG.times,freqs, squeeze(mean(d2p,3)) ,40, 'linecolor','non');
xlim([-400 800])
colorbar
set(gca,'clim',[-5 15])
set(gca, 'fontsize', 30)
title('Signal x power')
ylabel('Frequency (hz)')


subplot(233);
contourf(EEG.times,freqs, squeeze(mean(d2p2,3)),40, 'linecolor','non');
xlim([-400 800])
colorbar
set(gca,'clim',[-5 15])
set(gca, 'fontsize', 30)
title('Signal y power')


subplot(235);contourf(t2save, freqs, squeeze(tfg(1,2,:,:)), 40, 'linecolor','non'); colorbar
set(gca,'clim',[0 0.2])
set(gca, 'fontsize', 30)
xlabel('Time (ms)'), ylabel('Frequency (hz)')
title('Spectral GC y->x')


subplot(236);contourf(t2save, freqs, squeeze(tfg(1,1,:,:)), 40, 'linecolor','non'); cl = colorbar;
set(gca,'clim',[0 0.2])
cl.Label.String = 'GC estimate';
cl.FontSize = 30;
set(gca, 'fontsize', 30)
xlabel('Time (ms)'),
title('Spectral GC x->y')
cl.Label.String = 'GC estimate';


%% for correlations between granger and power


% dat2c = squeeze(mean(dbconverted1(:,:,times2saveidx,:),4));
% dat2c2 = squeeze(mean(dbconverted2(:,:,times2saveidx,:),4));
% 
% gc2c = squeeze(tfg(:,1,:,:));
% gc2c2 = squeeze(tfg(:,2,:,:));
% 
% 
% [rval, lagz] = power_pdc_ccorr(gc2c2, dat2c, freqs);
% figure;stem(lagz, rval,'b')
