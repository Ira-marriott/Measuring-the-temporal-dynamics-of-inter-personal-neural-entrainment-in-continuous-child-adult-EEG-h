
%  This code accompnies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk  

% also see - Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.



%% 

load sampleEEGdata.mat;

% simulation parameters 
f = 7;
srate = 256;
t  = -1:1/srate:1.5;
ntrials = 100;
d = zeros( length(t), ntrials);
s = zeros(length(t), ntrials);
noiseamp=0.5;

pnts = length(t);


for i = 1:ntrials

r = 2*pi*rand(1,320);

s(1:320,i) = sin(2*pi*f*t(1:320)+r); % constant phase

s(321:640,i) =  sin(2*pi*f*t(321:640)  );
% 

g = 2*pi*rand(1,320);

d(1:320,i) = sin(2*pi*f*t(1:320)+g); % constant phase

d(321:640,i) = sin(2*pi*f*t(321:640 ) ); % 

end


signal1 = s + noiseamp * randn(length(t),ntrials);

signal2 = d + noiseamp * randn(length(t),ntrials);

figure(3)
set(gcf,'color','w');


%% calculate phase synchronization

% extract angles from Hilbert transform
angles1 = angle(hilbert( signal1 ));
angles2 = angle(hilbert( signal2 ));

% show phase angle time series
% subplot(312)
% plot(EEG.times,angles1(1:end-1,1:5),'b')
% hold on
% plot(EEG.times,angles2(1:end-1,1:5),'g')
% ylabel('Phase angle')
% title('Phase angle over time')
% xlim([-400 1000])
% set(gca, 'ytick', [], 'fontsize', 30)


%% synchronization over time (in windows)

% analysis parameters
winsize = 100; % in points

% initialize
synchTS = zeros(1,pnts);

plv =  abs(mean(exp(1i* (angles1-angles2) ),2));

% loop over time points
% for ti=1:ntrials
%     
%     % figure out time points, considering edges
%     tidx = max(1,ti-winsize) : min(pnts,ti+winsize);
%     
%     
%     abs(mean(exp(1i* (phasediffs1) ),4));
%     
%     
%     
%     % compute synchronization for this center time points
%     synchTS(ti) = abs(mean(exp(1i*( angles1(tidx)-angles2(tidx) ))));
% end


%% plotting for figure 4(a and c)

subplot(221)
% plot(EEG.times,squeeze(mean(signal1(1:end-1,:),2)),'b', EEG.times,squeeze(mean(signal2(1:end-1,:),2)),'g')
plot(EEG.times,squeeze(mean(signal1(1:end-1,:),2)),'b','linew',2)
hold on 
plot(EEG.times,squeeze(mean(signal2(1:end-1,:),2))-3,'g','linew',2)

ylabel('Amplitude')
title('Raw amplitude over time')

frange = [6 9]; % in frequency
order  = round( 3*srate/frange(1) );
noiseamp = 0.5;

% filter kernel
filtkern = fir1(order,frange/(srate/2));

% signals are filtered noise
signal1 = filtfilt(filtkern,1,signal1);
signal2 = filtfilt(filtkern,1,signal2);
xlim([-400 1000])
set(gca, 'ytick', [], 'fontsize', 30)



% and plot
subplot(223)
plot(EEG.times,plv(1:end-1),'k','linew',2)
set(gca,'ylim',[0 1])
xlabel('Time (ms)'), ylabel('PLV')
title('Phase locking')
xlim([-400 1000])
set(gca, 'fontsize', 30)
hold on


%% for accompanying spectral analysis


% load sampleEEGdata.mat
% times = EEG.times;
% 
% clear EEG
% clear data
% 
% data(1,:,:) = signal1(1:end-1,:);
% data(2,:,:) = signal2(1:end-1,:);
% 
% EEG          = eeg_emptyset;
% EEG.data     = data;
% EEG.pnts     = size(data,2);
% EEG.srate    = srate;
% % EEG.chanlocs = EEG.chanlocs;
% EEG.nbchan   = size(data,1);
% EEG.event    = []; % overweighted events not in EEG.event structure
% EEG          = eeg_checkset(EEG);
% EEG.trials   = size(EEG.data,3);
% EEG.times = times;
% 
% 
% 
% freqs = linspace(2,40,38);
% freqbloom = 1.5;
% data = EEG.data;
% tf = zeros(2,length(freqs),EEG.pnts, ntrials);
% 
% 
% for f = 1:length(freqs)
%     
%     tmp= reshape(data,2,EEG.pnts*ntrials);
%     
%     for i=1:2
%         
%         tmpdat = tmp(i,:);
%         
% % create filtered noise
% frange = [freqs(f)-freqbloom freqs(f)+freqbloom]; % in frequency
% order  = round( 5*srate/frange(1) );
% 
% % filter kernel
% filtkern = fir1(order,frange/(srate/2));
% 
% % signals are filtered noise
% tf(i,f,:,:) = reshape( filtfilt(filtkern,1,tmpdat), EEG.pnts,ntrials);
% 
% 
%     end
%     
% end
% 
% angles1 = squeeze(angle(tf(1,:,:,:)));
% angles2 = squeeze(angle(tf(2,:,:,:)));
% 
% 
% plv = zeros(length(freqs), EEG.pnts);
% 
% for fi = 1:length(freqs)
%     
%     plv(fi,:,:) =  abs(mean(exp(1i* (angles1(fi,:,:)-angles2(fi,:,:)) ),3));
%     
% end
% 
% 
% for fi = 1:length(freqs)
%   
% % loop over time points
% for ti = 1:pnts
%     
% %     figure out time points, considering edges
%    
%     tidx = max(1,ti-winsize) : min(pnts,ti+winsize);
%      
% %     phasediffs1  = 
%     
% %     abs(mean(exp(1i* (phasediffs1) ),4));
%     
%    
% %     compute synchronization for this center time points
%     synchTS(fi,ti) = abs(mean(exp(1i*( angles1(tidx)-angles2(tidx) ))));
% end
% end
% 
% plv = movmean(plv, [51 51]);
% times2save = -700:100:1000;
% times2saveidx = dsearchn(EEG.times',times2save');
% 
% % figure(3)
% % subplot(325)
% % contourf(times2save, freqs, plv(:,times2saveidx) ,40 , 'linecolor', 'non'); 
% % colorbar
% % set(gca, 'fontsize',30)
% % title('Phase Transfer Entropy')
% % ylabel('Frequency (hz)')
% % xlabel('Time (ms)')
% 
% 
% figure(3)
% subplot(325)
% contourf(EEG.times, freqs, plv ,40 , 'linecolor', 'non'); 
% colorbar
% set(gca, 'fontsize',30)
% title('Phase Transfer Entropy')
% ylabel('Frequency (hz)')
% xlabel('Time (ms)')
% 
