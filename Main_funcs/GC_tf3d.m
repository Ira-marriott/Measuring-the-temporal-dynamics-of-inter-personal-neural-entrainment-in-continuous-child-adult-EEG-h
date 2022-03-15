
% Measuring the temporal dynamics of inter-personal neural entrainment in continuous child-adult EEG hyperscanning data

% https://doi.org/10.1016/j.dcn.2022.101093
% contact u1434978@uel.ac.uk

%% function for time varying granger

% Note that this function assumes that appropriate steps have been taken to
% ensure data is stationary and that an appropriate model order is used for
% the data


% Call as tf_granger = GC_tf3d(data1, data2, freqs, srate, order,twin)

% INPUTS:
% Data1 and data2 2d matrices of single channel EEG data e.g., in samples by trails format. To compute this for all possible channel pairing simply loop through channel combinations calling this function each time.

% Freqs is a vector of frequencies

% Order is the model order used in the autoreggressive model fit.

% Strate is the sampling frequency of the data.

% Twin is the size of the sliding window used to compute the GC estimate- given in ms.
 
 
% OUTPUT:
% Matrix of time-frequency varying GC influence in format 2 by frequencies by samples
% where tf_granger(1,:,:) is the influence of channel 1 to 2 and tf_granger(2,:,:) is the influence of channel 2 to 1



% This function uses code from
% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
% Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014

% https://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html


function tf_granger = GC_tf3d(data1, data2, freqs, srate, order,twin)

icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

% twin and order are specified in ms so here we get there indices in terms
% of data points
pnts = size(data1,1);
trials = size(data1,2);


% initialize
tf_granger=zeros(2,length(freqs),pnts);
twin_pnts = round(twin/(1000/srate));


%% compute time frequency granger causality 

for ti=1:pnts % loop over time points
    
    %     % figure out time points, considering edges
    tidx = max(1,ti-twin_pnts) : min(pnts,ti+twin_pnts);
    
    tempdata = permute(cat(3, data1,data2),[3 1 2]);
    
    % reshape tempdata for armorf function
    tempdat = reshape(tempdata(:,tidx,:),2,length(tidx)*trials);
    
    [A,SIG,~] = tsdata_to_var(tempdat,order,icregmode);
    [G,~] = var_to_autocov(A,SIG,[]);
    f = autocov_to_spwcgc(G,[]);
    
    tf_granger(1,:,ti) = squeeze(f(1,2,freqs));
    tf_granger(2,:,ti) = squeeze(f(2,1,freqs));
    
    
end % end time loop


end
