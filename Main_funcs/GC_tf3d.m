
function tf_granger = GPDC_tf3(data1, data2, times, times2save, freqs, srate, order,twin)

% twin is timewindow in ms
% data1 and data 2 should be 3D data all chans by times by trials
% freqs should be vector of frequencies of interest
% order is model order 

% this is time vector with middle number determining the amount if temporal
% smoothing - can be adjusted to whatever 
times2saveidx = dsearchn(times',times2save');

% twin and order are specified in ms so here we get there indices in terms
% of data points
timewin_points = round(twin/(1000/srate));

trials = size(data1,3);

% initialize
tf_granger=zeros(size(data1,1),2,length(freqs),length(times2save));

eeg = zeros(2, size(data1,1), length(times), size(data1,3));


%% 

for chani = 1:size(data1,1)
    
    % subtraction of ERP
eeg(1,chani,:,:) = bsxfun(@minus,data1(chani,:,:),mean(data1(chani,:,:),3));
eeg(2,chani,:,:) = bsxfun(@minus,data2(chani,:,:),mean(data2(chani,:,:),3));


for timei=1:length(times2save)
    
    % data from all trials in this time window 
    % (note that the ERP-subtracted data are used)
    timeidx = times2saveidx(timei)-floor(timewin_points/2):times2saveidx(timei)+floor(timewin_points/2)-mod(timewin_points+1,2);
    dataX = squeeze(eeg(:,chani,timeidx,:));
    
%     for triali=1:size(dataX,3)
%         dataX(1,:,triali) = zscore(detrend(squeeze(dataX(1,:,triali))));
%         dataX(2,:,triali) = zscore(detrend(squeeze(dataX(2,:,triali))));
% 
%         % At this point with real data, you might want to check for stationarity
%         % and possibly discard or mark data epochs that are non-stationary.
%     end
 
    % reshape tempdata for armorf function
    tempdat = reshape(dataX,2,length(timeidx)*trials);
    
    % fit AR model
    [Axy,E] = armorf(tempdat ,trials,timewin_points,order);
    
    % code below is adapted from bsmart toolbox function pwcausal.m  corrected covariance
    eyx = E(2,2) - E(1,2)^2/E(1,1); 
    exy = E(1,1) - E(2,1)^2/E(2,2); 
    N = size(E,1);
    
    for fi=1:length(freqs)
        
        % transfer matrix (note the similarity to Fourier transform)
        H = eye(N);
        
        for m = 1:order
            H = H + Axy(:,(m-1)*N+1:m*N)*exp(-1i*m*2*pi*freqs(fi)/srate);
            
        end
        
        Hi = inv(H); 
        S  = H\E*Hi'/srate;
        
        % granger prediction per frequency
        tf_granger(chani,1,fi,timei) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/srate) );
        tf_granger(chani,2,fi,timei) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/srate) );
    end
    
end % end time loop

end% end channel loop
