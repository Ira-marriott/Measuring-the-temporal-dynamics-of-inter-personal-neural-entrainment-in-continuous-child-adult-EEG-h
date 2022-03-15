

% Measuring the temporal dynamics of inter-personal neural entrainment in continuous child-adult EEG hyperscanning data

% https://doi.org/10.1016/j.dcn.2022.101093
% contact u1434978@uel.ac.uk

%% Function for narrow band filering of EEG.

% call as inst_freq = InstFreq(data, srate)
 
% INPUT:
% data is vector or a matrix (either 2d or 3d) of hilbert transformed data – but where one dimension is time.
 
% Strate is the sampling frequency of the data.
 
% OUTPUT:
% time-varying instantaneous frequency of input



function inst_freq = InstFreq(data, srate)


% function for outputting instantaneous frequency from hilbert transformed data
inst_freq = srate*diff(unwrap(angle(data)))/(2*pi);

