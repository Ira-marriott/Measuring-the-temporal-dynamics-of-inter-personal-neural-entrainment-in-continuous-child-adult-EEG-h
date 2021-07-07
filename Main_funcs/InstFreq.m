

%  This code accompanies the paper "Using dual EEG to analyse event-locked changes in child-adult neural connectivity"

% ref Marriott Haresign, Phillips, Whitehorn, Goupil, & Wass, 2021
% contact u1434978@uel.ac.uk

%% Function for narrow band filering of EEG.

% call as inst_freq = InstFreq(data, srate)

% % Where data is vector or a matrix (either 2d or 3d) of hilbert transformed data - but
% where one dimension is time.

% Strate is the sampling frequency of the data.

% OUTPUT is time-varying instantaneous frequency of input


function inst_freq = InstFreq(data, srate)


% function for outputting instantaneous frequency from hilbert transformed data

inst_freq = srate*diff(unwrap(angle(data)))/(2*pi);

