function inst_freq = InstFreq(data, fs)


% fucntion for outputting instantaneous frequency from hilbert transformed data

inst_freq = fs*diff(unwrap(angle(data)))/(2*pi);

