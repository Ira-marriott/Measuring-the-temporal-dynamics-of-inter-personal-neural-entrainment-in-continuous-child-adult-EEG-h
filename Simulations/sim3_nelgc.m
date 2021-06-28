
%  This code accompnies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk  


% also see - Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.

% and see - Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference. Journal of neuroscience methods, 223, 50-68.



%% Simulation 3 non event locked changes in granger causality.


load sampleEEGdata.mat
times = EEG.times;

% simulation parameters
freq1 = 10;
freq2 = 2;
srate = 256;
noiseamp = 0.5;
time  = -1:1/srate:1.5;
ntrials = 10;

clear x
for ti = 1:10
    
    % define x
    x(ti,:) = sin(2*pi*freq1.*time) + sin(2*pi*freq2.*time) + noiseamp * randn(size(time))/2;
    
    % define y by x
    % y = [0 .2];
    
    for i=2:length(x)
        
        y(ti,i) =  x(ti,i-1) + randn/i;
    end
end


%% fit to eeglab structure just for ease of analysis

clear EEG
clear data

data(1,:,:) = y(:,1:end-1)';
data(2,:,:) = x(:,1:end-1)';

EEG          = eeg_emptyset;
EEG.data     = data;
EEG.pnts     = size(data,2);
EEG.srate    = srate;
% EEG.chanlocs = EEG.chanlocs;
EEG.nbchan   = size(data,1);
EEG.event    = []; % overweighted events not in EEG.event structure
EEG          = eeg_checkset(EEG);
EEG.trials   = size(EEG.data,3);
EEG.times = times;


regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

morder = 5;
fs = 256;
twin  =  200;
twin_pnts = round(twin/(1000/fs));

times2save = -800:50:1100;
times2saveidx = dsearchn(EEG.times',times2save');


tv_gc = zeros(2, length(times2save));

for i = 1:length(times2save)
    
    X = EEG.data([1 2],times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),:);
    
    % taken from MVGC toolbox 
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    var_info(info,true);
    %
    %
    F = autocov_to_pwcgc(G);
    
    tv_gc(1,i) = F(1,2);
    tv_gc(2,i) = F(2,1);
    

end

%% plotting for figure 5


figure(1)
set(gcf,'color','w');

subplot(121)
plot(times,x(1,1:end-1),'b','linew',2)
hold on
plot(times,y(1,1:end-1)-3,'g','linew',2)
legend({'x', 'y'})
xlim([-400 1100])
xticklabels({'0.5', '1', '1.5', })
set(gca, 'ytick',[], 'fontsize',30)
title('Raw Amplitude')
xlabel('Time (s)')


figure(1)
subplot(122);plot(times2save,tv_gc(1,:)-7,'b','linew',2)
hold on
plot(times2save,tv_gc(2,:),'g','linew',2)
legend({'x->y','y->x' })
title('Time Domain GC')
set(gca, 'fontsize',30)
ylabel('GC estimate')
xlim([-400 1100])
xticklabels({'0.5', '1', '1.5', })
hold on
xlabel('Time (s)')


%% this section is for the spectral analysis of the above
%
% pnts= length(EEG.times);
% % freqs = 1:40;
% bt = dsearchn(EEG.times',[-700 -400]');
%
% % data1 = erp(1,1:end,:);
% % data2 = erp(2,1:end,:);
% %
% % [~, ~,pow1]=conv_mf2(data(1,:,:), srate, freqs);
% % [~, ~,pow2]=conv_mf2(data(2,:,:), srate, freqs);
%
% freqs = linspace(2,40,38);
% freqbloom = 1.5;
%
%
% tf = zeros(2,length(freqs),pnts, ntrials);
%
% for f = 1:length(freqs)
%
%     tmp= reshape(data,2,pnts*ntrials);
%
%     for i=1:2
%
%         tmpdat = tmp(i,:);
%
% % create filtered noise
% frange = [freqs(f)-freqbloom freqs(f)+freqbloom]; % in frequency
% order  = round( 3*srate/frange(1) );
%
% % filter kernel
% filtkern = fir1(order,frange/(srate/2));
%
% % signals are filtered noise
% tf(i,f,:,:) = reshape( filtfilt(filtkern,1,tmpdat), pnts,ntrials);
%
%
%     end
%
% end
%
%
% pow1 = abs(hilbert(tf(1,:,:,:))).^2;
% pow2 = abs(hilbert(tf(2,:,:,:))).^2;
%
%
%
% baseline_power1 = mean(pow1(:,:,bt(1):bt(2),:),3);  dbconverted1 = 10*log10( bsxfun(@rdivide,pow1,baseline_power1));
% baseline_power2 = mean(pow2(:,:,bt(1):bt(2),:),3);  dbconverted2 = 10*log10( bsxfun(@rdivide,pow2,baseline_power2));
%
% % times = -2500:1000/srate:2500;
% times2save = -400:200:800;
% times2saveidx = dsearchn(EEG.times',times2save');
%
%
%
% d2p = squeeze(dbconverted1);
% d2p2 = squeeze(dbconverted2);
%
%
% figure(1);
% subplot(232);contourf(EEG.times,freqs, squeeze(mean(d2p,3)) ,40, 'linecolor','non');
% xlim([-400 800])
% colorbar
% set(gca,'clim',[-5 15])
% set(gca, 'fontsize', 20)
% title('Signal x power')
%  ylabel('Frequency (hz)')
%
%
% % d2p2 = squeeze(dbconverted2);
% subplot(233);contourf(EEG.times,freqs, squeeze(mean(d2p2,3)),40, 'linecolor','non');
% xlim([-400 800])
% colorbar
% set(gca,'clim',[-5 15])
% set(gca, 'fontsize', 20)
% title('Signal y power')
%
%
%
% order = 5;
% twin = 500;
% freqz=freqs;
%
%
%
% tfg = GPDC_tf3(data(1,:,:) ,data(2,:,:), EEG.times,times2save, freqz, srate, order,twin);
%
%
% subplot(235);contourf(times2save, freqz, squeeze(tfg(1,2,:,:)), 40, 'linecolor','non'); colorbar
% % set(gca,'clim',[0 0.2])
% set(gca, 'fontsize', 20)
% xlabel('Time (ms)'), ylabel('Frequency (hz)')
% title('Spectral GC y->x')
%
%
%
%
% subplot(236);contourf(times2save, freqz, squeeze(tfg(1,1,:,:)), 40, 'linecolor','non'); colorbar
% % set(gca,'clim',[0 0.2])
% set(gca, 'fontsize', 20)
% xlabel('Time (ms)'),
% title('Spectral GC x->y')
%
%
