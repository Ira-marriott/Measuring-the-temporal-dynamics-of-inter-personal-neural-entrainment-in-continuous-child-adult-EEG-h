%% This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk



%% function for time varying granger

% Note that this function assumes that appropriate steps have been taken to
% ensure data is stationary and that an appropriate model order is used for
% the data


% Call as tv_gc = TdGc_3d(data1, data2, winsize, morder, srate)

% INPUTS:
% Where data1 and data2 2d matricies of single channel EEG data e.g., in samples by
% trails format. To compute this for all possible channel pairing simply loop throught
% channel combinations calling this fucniton each time.
% Winsize is the size of the sliding window used to compute the GC estimate - given in ms.
% Morder is the model order used in the autoreggressive model fit.
% Strate is the sampling frequency of the data.


% OUTPUT:
% is matrix of time varying GC inlfuence in format 2 by samples by
% Where tv_gc(1,:,:) is the influecne of channel 1 to 2
%  and tv_gc(2,:,:) is the influecne of channel 2 to 1


% This function uses code from
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014

% https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html

%% 


function tv_gc = GC_td3d(data1, data2, winsize, morder, srate, times2saveidx)


regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
acmaxlags = [];

nbtrials = size(data1,2);
twin_pnts = ceil((srate/1000)*winsize);

tv_gc = zeros(2, length(times2saveidx));


eegdata = [];

eegdata(1,:,:,:) = data1;
eegdata(2,:,:,:) = data2;


for ti=1:length(times2saveidx)
    
    % figure out time points, considering edges
    tempdata = squeeze(eegdata(:,times2saveidx(ti)-floor(twin_pnts/2):times2saveidx(ti)+floor(twin_pnts/2)-mod(twin_pnts+1,2),:));

    % reshape tempdata for armorf
    tempdata = reshape(tempdata,2,twin_pnts*nbtrials);

    [A,SIG] = tsdata_to_var(tempdata,morder,regmode);
    
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    var_info(info,true);
    
    F = autocov_to_pwcgc(G);
    
    tv_gc(1,ti) = F(1,2);
    tv_gc(2,ti) = F(2,1);
    
    
end


end