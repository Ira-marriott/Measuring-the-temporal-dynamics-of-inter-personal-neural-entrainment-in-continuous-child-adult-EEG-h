
%  This code accompnies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk  


% also see - Cohen, M. X. (2014). Analyzing neural time series data: theory and practice. MIT press.


%% simulation for event locked changes in phase transfer entropy

% simulation parameters
ntrials = 100;
f = 7;
srate = 256;
t  = -1:1/srate:1.5;
d = zeros( length(t), ntrials);
s = zeros(length(t), ntrials);
noiseamp=0.5;
pnts = length(t);


for i = 1:ntrials

r = 2*pi*rand(1);

s(1:320,i) = sin(2*pi*f*t(1:320)+r); % constant phase

s(321:640,i) =  sin(2*pi*f*t(321:640)  );
% 


d(1:370,i) = sin(2*pi*f*t(1:370)+r); % constant phase

d(371:640,i) = sin(2*pi*f*t(371:640 ) ); % 

end


data1 = s + noiseamp * randn(length(t),ntrials);
data2 = d + noiseamp * randn(length(t),ntrials);


%% load into eeglab structure

load sampleEEGdata.mat
times = EEG.times;

clear EEG
clear data

data(1,:,:) = data1(1:end-1,:);
data(2,:,:) = data2(1:end-1,:);

EEG          = eeg_emptyset;
EEG.data     = data;
EEG.pnts     = size(data,2);
EEG.srate    = srate;
EEG.nbchan   = size(data,1);
EEG.event    = []; 
EEG          = eeg_checkset(EEG);
EEG.trials   = size(EEG.data,3);
EEG.times = times;


order = 15;
twin = 200;
freqz=2:40;

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

[dPTE, PTE] = PhaseTE_MF(X', [], []);

pte(1,i,ti) = PTE(1,2);
pte(2,i,ti) = PTE(2,1);

end
end

%% plotting for figure 4(b and d)


figure(3)
subplot(222)
plot(EEG.times,squeeze(mean(EEG.data(1,:,:),3)),'b','linew',2 )
hold on
plot(EEG.times,squeeze(mean(EEG.data(2,:,:),3))-3,'g','linew',2 )
xlim([-400 1000])
set(gca, 'fontsize',30,'ytick',[])
legend({'x', 'y'})
title('Raw amplitude over time')



subplot(224)
plot( times2save, squeeze(mean(pte(1,:,:),3)),'b', 'linew',2);
hold on 
plot( times2save, squeeze(mean(pte(2,:,:),3)),'g', 'linew',2);
set(gca, 'fontsize',30)
title('Phase Transfer Entropy')
xlabel('Time (ms)')
xlim([-400 1000])
legend({'x->y', 'y->x'})
ylabel('PTE estimate')


%% for time frequency pte - not used in main text 


freqs = linspace(2,40,38);
freqbloom = 1.5;
data = EEG.data;
tf = zeros(2,length(freqs),EEG.pnts, ntrials);


for f = 1:length(freqs)
    
    tmp= reshape(data,2,EEG.pnts*ntrials);
    
    for i=1:2
        
        tmpdat = tmp(i,:);
        
% create filtered noise
frange = [freqs(f)-freqbloom freqs(f)+freqbloom]; % in frequency
order  = round( 3*srate/frange(1) );

% filter kernel
filtkern = fir1(order,frange/(srate/2));

% signals are filtered noise
tf(i,f,:,:) = reshape( filtfilt(filtkern,1,tmpdat), EEG.pnts,ntrials);


    end
    
end

%% for time frequency pte - not used in main text 

pte = zeros(2,length(freqs), length(times2save), ntrials);


for fi =1:length(freqs)
    
    for ti = 1:ntrials

for i = 1:length(times2save)
    
X = squeeze(tf([chan1 chan2],fi,times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2),ti));

[dPTE, PTE] = PhaseTE_MF(X', [], []);

pte(1,fi,i,ti) = PTE(1,2);
pte(2,fi,i,ti) = PTE(2,1);

end
    end
end

plv =  abs(mean(exp(1i* (angle(hilbert(EEG2.data(1,:,:)))-angle(hilbert(EEG2.data(2,:,:))) )),3));



baseline_pte =  mean(pte(:,:,3:7,:),3);  
pteZ =  bsxfun(@minus,pte,baseline_pte);





