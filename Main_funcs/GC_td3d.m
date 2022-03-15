% Measuring the temporal dynamics of inter-personal neural entrainment in continuous child-adult EEG hyperscanning data

% https://doi.org/10.1016/j.dcn.2022.101093
% contact u1434978@uel.ac.uk


%% function for time varying granger

% Note that this function assumes that appropriate steps have been taken to
% ensure data is stationary and that an appropriate model order is used for
% the data


% Call as tv_gc = GC_td3d(data1, data2, winsize, morder, srate)

 
% INPUTS:
% Where data1 and data2 2d matrices of single channel EEG data e.g., in samples by trials format. To compute this for all possible channel pairing simply loop through channel combinations calling this function each time.

% Winsize is the size of the sliding window used to compute the GC estimate - given in ms.

% Morder is the model order used in the autoreggressive model fit.

% Strate is the sampling frequency of the data.
  
% OUTPUT:
% Matrix of time varying GC inlfuence in format 2 by samples by where tv_gc(1,:) is the influence of channel 1 to 2 and tv_gc(2,:) is the influence of channel 2 to 1.

% This function uses code from
% [1] L. Barnett and A. K. Seth <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%  Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014

% https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html


%% 

function tv_gc = GC_td3d(data1, data2, winsize, morder, srate)


regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
acmaxlags = [];

pnts = size(data1,1);
trials = size(data1,2);

twin_pnts = ceil((srate/1000)*winsize);

tv_gc = zeros(2, pnts);


for ti=1:pnts % loop over time points
    %
    %     % figure out time points, considering edges
    tidx = max(1,ti-twin_pnts) : min(pnts,ti+twin_pnts);
    
    tempdata = permute(cat(3, data1,data2),[3 1 2]);
    
    tempdat = reshape(tempdata(:,tidx,:),2,length(tidx)*trials);

    [A,SIG] = tsdata_to_var(tempdat,morder,regmode);
    
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    var_info(info,true);
    
    F = autocov_to_pwcgc(G);
    
    tv_gc(1,ti) = F(1,2);
    tv_gc(2,ti) = F(2,1);

end


end