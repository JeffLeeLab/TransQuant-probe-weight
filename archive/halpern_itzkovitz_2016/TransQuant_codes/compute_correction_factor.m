function [W,L,N]=compute_correction_factor()

% Computes the probe weight factor, W, gene length, L and probe
% localization profile, N

% default values
W=0.5;
L=6000;
N=ones(1,6000);
[filename, path]=uigetfile([pwd '\*.*'],'Select UCSC sequence file');
if filename==0
    return;
end
UCSC_full_sequence_file=[path filename];

[filename, path]=uigetfile([path '\*.*'],'Select probe file');
if filename==0
    return;
end
probe_file=[path filename];
track_file='track_temp.txt';

output_UCSC_probe_track(UCSC_full_sequence_file,probe_file, track_file);
[W,L,N]=compute_weight_factor(UCSC_full_sequence_file,track_file);