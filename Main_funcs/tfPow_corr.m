% power_corr

function corr_ts = tfPow_corr(a,b)

% assumes a and b and frequencies by time by trails data

freqs = size(a,1);
pnts =size(a,2);
ntrials = size(a,3);


% correlation at original sample rate
corr_ts = zeros(length(freqs),pnts,ntrials);

for i = 1:freqs
    
    for ti=1:pnts
        
        x = corrcoef(a(i,ti,:), b(i,ti,:));
        corr_ts(i,ti,:) = x(2,1);
        
    end
    
end


