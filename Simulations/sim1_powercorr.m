
%  This code accompnies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk  


% also see - Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.


%% section 1  - Single trial power correlations


% simulation parameters
sfreq = 7;
ntrials = 100;
srate = 256;
time = -1:1/srate:1.5;
gwidth = 0.1;
noiseamp = 0.5; % noise standard deviation
pnts  = length(time);


signal = zeros(2,length(time),ntrials);

% non phase locked activity
for ti = 1:ntrials
    
    % create signal
    swave  = sin(2*pi*sfreq*time + 2*pi*rand );
    gausw  =  exp( -4*log(2)*(time).^2 / gwidth^2 );
    signal(1,:,ti) = swave .* gausw;
    
    gausw2  =  exp( -4*log(2)*(time-.3).^2 / gwidth^2 );
    signal(2,:,ti) = swave .* gausw2;
    
end

clear data
% create data and multiple channels plus noise
data(1,:,:) = squeeze(signal(1,:,:)) + noiseamp * randn(length(time),ntrials);
data(2,:,:) = squeeze(signal(2,:,:)) + noiseamp * randn(length(time),ntrials);


%%  filter data

freqs = linspace(2,25,23);
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

%% fit to EEG lab template - this will help for event locked analysis

load sampleEEGdata.mat
times = EEG.times;

clear EEG

datax(1,:,:) = data(1,1:end-1,:);
datax(2,:,:) = data(2,1:end-1,:);

EEG          = eeg_emptyset;
EEG.data     = datax;
EEG.pnts     = size(datax,2);
EEG.srate    = srate;
EEG.nbchan   = size(datax,1);
EEG.event    = [];
EEG          = eeg_checkset(EEG);
EEG.trials   = size(EEG.data,3);
EEG.times = times;


bt = dsearchn(EEG.times',[-700 -400]');
pow = abs(hilbert(tf)).^2;

baseline_power1 = mean(pow(1,:,bt(1):bt(2),:),3);  dbconverted1 = 10*log10( bsxfun(@rdivide,pow(1,:,:,:),baseline_power1));
baseline_power2 = mean(pow(2,:,bt(1):bt(2),:),3);  dbconverted2 = 10*log10( bsxfun(@rdivide,pow(2,:,:,:),baseline_power2));

a = squeeze(dbconverted1(1,:,:,:));
b = squeeze(dbconverted2(1,:,:,:));

% correlation at original sample rate
corr_ts = zeros(length(freqs),pnts,ntrials);

for i = 1:length(freqs)
    
    for ti=1:pnts
        
        x = corrcoef(a(i,ti,:), b(i,ti,:));
        corr_ts(i,ti,:) = x(2,1);
        
    end
    
    
end


t2save = -400:100:800;
tidx = dsearchn(EEG.times', t2save');
twin = 300;

corr_ts2 = zeros(length(freqs), length(tidx),ntrials);
ps= zeros(length(freqs), length(tidx),ntrials);


% downsample data
A = a(:,tidx,:);
B = b(:,tidx,:);

% downsampled correlation
for i = 1:length(freqs)
    
    for ti=3:length(tidx)-3
        
        %     x = corrcoef(a(i,tidx(ti)-twin/2:tidx(ti)+twin/2,:), b(i,tidx(ti)-twin/2:tidx(ti)+twin/2,:));
        
        [x,p] = corrcoef(A(i,ti-2,:), B(i,ti+2,:));
        
        corr_ts2(i,ti,:) = x(2,1);
        ps(i,ti,:) = p(2,1);
        
    end
    
    
end


% signficant correaltions
z = squeeze(ps(:,:,1));
q = z<0.05;



% downsample using moving average for plotting only
A2 = movmean(A,8,2);
B2 = movmean(B,8,2);


%% plotting for power correlations - figure 2


f = figure(1);
set(gcf,'color','w');

subplot(3, 2, 1)
contourf(EEG.times, freqs, squeeze(mean(dbconverted1(1,:,1:end-1,:),4)) ,40,'linecolor','non');
ylabel('Frequency (hz)')
set(gca, 'fontsize', 30, 'clim', [-8 8])
colorbar;
title('Signal x power')
xlim([-300 600])


subplot(3, 2, 2)
contourf(EEG.times, freqs, squeeze(mean(dbconverted2(1,:,1:end-1,:),4)) ,40,'linecolor','non');
set(gca, 'fontsize', 30, 'clim', [-8 8])
dl = colorbar;
title('Signal y power')
xlim([-300 600])
dl.Label.String = 'Power(db)';

subplot(3, 2, 3)
contourf(t2save , freqs, squeeze(mean(A2,3)),40,'linecolor','non');
title('Signal x downsampled')
ylabel('Frequency (hz)')
xlim([-300 600])
colorbar;
set(gca, 'fontsize', 30, 'clim', [-5 2])


subplot(3, 2, 4)
contourf(t2save , freqs, squeeze(mean(B2,3)) ,40,'linecolor','non');
title('Signal y downsampled')
xlim([-300 600])
cl = colorbar;
set(gca, 'fontsize', 30, 'clim', [-5 2])
cl.Label.String = 'Power(db)';


subplot(3, 2, 5)
contourf(EEG.times , freqs, squeeze(mean(corr_ts(:,1:end-1,:),3)),40,'linecolor','non');
title('PC at original sampling rate')
xlabel('Time (ms)'), ylabel('Frequency (hz)')
set(gca, 'clim', [-0.2 0.4])
xlim([-300 600])
colorbar;
set(gca, 'fontsize', 30)


subplot(3, 2, 6)
contourf(t2save , freqs, squeeze(mean(corr_ts2,3)),40,'linecolor','non');
title('PC between downsampled signals')
xlabel('Time (ms)')
xlim([-300 600])
cl = colorbar;
cl.Label.String = 'Correlation';
set(gca, 'fontsize', 30,'clim', [-0.2 0.4])

