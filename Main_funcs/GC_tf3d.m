
%% This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk



%% function for time varying granger

% Note that this function assumes that appropriate steps have been taken to
% ensure data is stationary and that an appropriate model order is used for
% the data


% Call as tv_gc = TdGc_3d(data1, data2, winsize, morder, srate)

% INPUTS:
% Data1 and data2 2d matricies of single channel EEG data e.g., in samples by
% trails format. To compute this for all possible channel pairing simply loop throught
% channel combinations calling this fucniton each time.
% Twin is the size of the sliding window used to compute the GC estimate
% - given in ms.
% Morder is the model order used in the autoreggressive model fit.
% Strate is the sampling frequency of the data.

% times is time vector with middle number determining the amount if temporal
% smoothing - can be adjusted to whatever


% OUTPUT:
% is matrix of time varying GC inlfuence in format 2 by samples by
% Where tv_gc(1,:,:) is the influecne of channel 1 to 2
%  and tv_gc(2,:,:) is the influecne of channel 2 to 1


% This function uses code from
% Jie Cui, Lei Xu, Steven L. Bressler, Mingzhou Ding, Hualou Liang,
% BSMART: a Matlab/C toolbox for analysis of multichannel neural time series, Neural Networks, 21:1094 - 1104, 2008.

% https://brain-smart.org/


%%


function tf_granger = GC_tf3d(data1, data2, freqs, srate, order,twin)


% twin and order are specified in ms so here we get there indices in terms
% of data points
twin_pnts = round(twin/(1000/srate));

pnts = size(data1,1);
trials = size(data1,2);


% initialize
tf_granger=zeros(2,length(freqs),pnts);


%%


for ti=1:pnts % loop over time points
    
    % figure out time points, considering edges
    tidx = max(1,ti-twin_pnts) : min(pnts,ti+twin_pnts);
    
    
    tempdata = permute(cat(3, data1,data2),[3 1 2]);
%     tempdata = cat(1,data1,data2);
    
    % reshape tempdata for armorf function
    tempdat = reshape(tempdata(:,tidx,:),2,length(tidx)*trials);
    
    
    % fit AR model
    [Axy,E] = armorf(tempdat ,trials,twin_pnts,order);
    
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
        tf_granger(1,fi,ti) = log( abs(S(2,2))/abs(S(2,2)-(Hi(2,1)*exy*conj(Hi(2,1)))/srate) );
        tf_granger(2,fi,ti) = log( abs(S(1,1))/abs(S(1,1)-(Hi(1,2)*eyx*conj(Hi(1,2)))/srate) );
    end
    
end % end time loop

% end% end channel loop
