%% Example for obtaining pte


% This section of the tutorial relies on the EEGlab structure
% EEGlab must be downloaded first https://sccn.ucsd.edu/eeglab/download.php


% Load in sample dual EEG data provided or load in your own data
load infEEG.mat; EEG_infant=EEG;
load mumEEG.mat; EEG_mum = EEG;


% initialise variables
order = 15;
twin = 200;
freqz=2:20;
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

% get infant and mum data from specified channels
data1 = double(squeeze(EEG_infant.data(chan1,:,:)));
data2 = double(squeeze(EEG_mum.data(chan2,:,:)));

% filter the data
filtsig = FiltHilb(data1, freqz, freqbloom, srate);
filtsig2 = FiltHilb(data2, freqz, freqbloom, srate);



% initialise PTE output
pte = zeros(2, length(freqz), length(times2save), ntrials);

for fi = 1:length(freqz)
    %     loop over frequencies
    
    for ti = 1:ntrials
        %     loop over trials
        
        for i = 1:length(times2save)
            %     loop over specified times
            
            %     define time indx
            tidx = times2saveidx(i)-floor(twin_pnts/2):times2saveidx(i)+floor(twin_pnts/2);
            
%             get mum and infant data from that time segment 
            X = squeeze(cat(1,filtsig(fi,tidx,ti), filtsig2(fi,tidx,ti))) ;
                        
            % main fucntion for this section compute phase transfer entropy on narrow band filtered data
%             data is channels (2) x samples format
            [dPTE, PTE] = PhaseTE_MF(X', [], []);
            
            pte(1,fi,i,ti) = PTE(1,2);
            pte(2,fi,i,ti) = PTE(2,1);
            
        end
    end
end


%% Plotting

figure(3)
subplot(211)
plot(EEG.times,squeeze(mean(EEG.data(1,:,:),3)),'b','linew',2 )
hold on
plot(EEG.times,squeeze(mean(EEG.data(2,:,:),3))-3,'g','linew',2 )
xlim([-400 1000])
set(gca, 'fontsize',30,'ytick',[])
legend({'x', 'y'})
title('Raw amplitude over time')


% Only 2D plotting here but could also plot time-frequency plot of PTE output
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