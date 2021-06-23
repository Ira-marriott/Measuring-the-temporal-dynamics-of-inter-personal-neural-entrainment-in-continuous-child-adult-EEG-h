% function for time varying granger

% Adapted from

% [1] L. Barnett and A. K. Seth,
% <http://www.sciencedirect.com/science/article/pii/S0165027013003701 The MVGC
%     Multivariate Granger Causality Toolbox: A New Approach to Granger-causal
% Inference>, _J. Neurosci. Methods_ 223, 2014



function tv_gc = TdGc_3d(data1, data2, winsize, morder, srate)


regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
acmaxlags = [];


pnts = size(data1,1);
nbtrials = size(data1,2);

twin_pnts = ceil((srate/1000)*winsize);

tv_gc = zeros(2, pnts);


for traili = 1:nbtrials
    
    for ti=1:pnts % loop over time points
        
        % figure out time points, considering edges
        tidx = max(1,ti-twin_pnts) : min(pnts,ti+twin_pnts);
        
        tmpdat(1,:,:) = data1(tidx,:);
        tmpdat(2,:,:) = data2(tidx,:);
        
        [A,SIG] = tsdata_to_var(tmpdat,morder,regmode);
        
        [G,~] = var_to_autocov(A,SIG,acmaxlags);
        
%         var_info(info,true);
        
        F = autocov_to_pwcgc(G);
        
        tv_gc(1,ti) = F(1,2);
        tv_gc(2,ti) = F(2,1);
        
        clear tmpdat
        
    end
end





