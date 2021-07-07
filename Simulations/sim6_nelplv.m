

%  This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk


% also see - Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.

% and see - Barnett, L., & Seth, A. K. (2014). The MVGC multivariate Granger causality toolbox: a new approach to Granger-causal inference. Journal of neuroscience methods, 223, 50-68.



%% simulation of frequency change to increased plv

% signal properties
load sampleEEGdata.mat
EEG = pop_resample(EEG,1000);

srate =  1000;
% time  = (0:pnts-1)/srate;
time  = -1:1/srate:1.5;
pnts = length(time);
noiseamp = 0.5;
ntrials = 2;

% time series of instantaneous frequencies, which is created by
% interpolating a few random datapoints
x = linspace(6,9,10);
freqTS = abs( interp1(linspace(time(1),time(end),10),x,time,'spline') );

% this is the construction of a frequency-varying sine wave
centfreq = mean(x);
r = 2*pi*rand(1,length(time));

k        = (centfreq/srate)*2*pi/centfreq;
Osig   = sin(2*pi.*centfreq.*time + k*cumsum(freqTS-centfreq)) + noiseamp * randn(ntrials,length(time));


y = linspace(12,9,10);
freqTS = abs( interp1(linspace(time(1),time(end),10),y,time,'spline') );

% this is the construction of a frequency-varying sine wave
centfreq = mean(y);
k        = (centfreq/srate)*2*pi/centfreq;
Osig2   = sin(2*pi.*centfreq.*time + k*cumsum(freqTS-centfreq)) + noiseamp * randn(ntrials,length(time));


frange = [6 9]; % in frequency
frange2 = [9 12]; % in frequency

order  = round(1*srate/frange(1) );
order2  = round(1*srate/frange2(1) );

% % filter kernel
filtkern = fir1(order,frange/(srate/2));
filtkern2 = fir1(order,frange2/(srate/2));

% % signals are filtered noise
sig = filtfilt(filtkern,1,Osig');
sig2 = filtfilt(filtkern2,1,Osig2');

% hilbert transform of signal
hilsig = hilbert(sig);

% hilbert transform of signal
hilsig2 = hilbert(sig2);


%%

angles1 = angle(hilsig) ;
angles2 = angle(hilsig2) ;


% analysis parameters
winsize = 100; % in points

% initialize
synchTS = zeros(1,pnts);

plv =  abs(mean(exp(1i* (angles1-angles2) ),2));

% loop over time points
for ti = 1:pnts
    
    %     figure out time points, considering edges
    
    tidx = max(1,ti-winsize) : min(pnts,ti+winsize);
    
    
    %     compute synchronization for this center time points
    synchTS(ti) = abs(mean(exp(1i*( angles1(tidx)-angles2(tidx) ))));
end

%% plotting - see figure 6


figure(3), clf
set(gcf,'color','w');

subplot(411)
plot(EEG.times,Osig(:,1:end-1), 'b')
hold on
plot(EEG.times,Osig2(:,1:end-1),'g')
ylabel('Amplitude')
% set(gca,'ylim',[-1.1 1.1])
set(gca,'ytick',[],'fontsize',30)
title('Raw amplitude over time')
xlim([-400 1100])


subplot(412)
plot(EEG.times,angle(hilsig(1:end-1,1)), 'b', 'linew',2)
hold on
plot(EEG.times,angle(hilsig2(1:end-1,1)),'g', 'linew',2)
ylabel('Phase')
set(gca, 'fontsize',30)
set(gca,'ytick',[],'xtick',[])
title('Phase angle over time')
xlim([-400 1100])

freqsld = InstFreq(hilsig, srate);
freqsld2 = InstFreq(hilsig2, srate);


subplot(413)
plot(EEG.times,freqsld(:,1),'b', 'linew',2)
hold on
plot(EEG.times,freqsld2(:,1),'g', 'linew',2)
% hold on
% plot(time,freqTS,'r')
set(gca,'ylim',[0 15] , 'fontsize',30,'xtick',[])
legend({'Signal x';'Signal y'})
title('Frequency over time')
ylabel('Frequency hz')
xlim([-400 1100])


figure(3)
subplot(414)
plot(EEG.times,synchTS(1:end-1),'k', 'linew',2)
ylabel('PLV')
set(gca, 'fontsize',30)
title('Synchronisation over time')
xlim([-400 1100])
xlabel('Time(s)')
xticklabels({'0.5', '1', '1.5', })



%%
