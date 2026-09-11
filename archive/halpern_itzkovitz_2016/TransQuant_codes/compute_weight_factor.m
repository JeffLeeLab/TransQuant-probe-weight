function [W,L,N]=compute_weight_factor(UCSC_full_sequence_file,UCSC_track_file)

% [W,L,N]=compute_weight_factor(UCSC_full_sequence_file,UCSC_track_file)
% The function computes 1/L * integral_0^L(N(x)) where N(x) is the number
% of probes found until coordinate x
% Input - UCSC_full_sequence_file - the fasta file of the full UCSC gene
%         (including exons and introns and UTRs
%         UCSC_track_file - the UCSC track file produced by
%         output_UCSC_probe_track.m
% Output - W - the correction factor for dividing the TS intensity. An
%            average Pol2 will have an intensity of W*(intensity of mature mRNA)
%          L - the length of the gene transcribed.
%          N - the number of probes up to coordinate i along the gene.
% Note: The function considers the strand transcribed when computing W and
%       N


if nargin<1,
    [filename, path]=uigetfile([pwd '\*.*'],'Select UCSC sequence file');
    UCSC_full_sequence_file=[path filename];
end

if nargin<2,
    [filename, path]=uigetfile([pwd '\*.*'],'Select track');
    UCSC_track_file=[path filename];
end

[a,b]=fastaread(UCSC_full_sequence_file);
index=findstr(a,':');
tmp=a((index+1):end);
index=findstr(tmp,'-');
strt=str2num(tmp(1:(index-1)));
tmp=tmp((index+1):end);
index=findstr(tmp,' ');
stp=str2num(tmp(1:(index-1)));
index=findstr(tmp,'strand=');
if strcmp(tmp(index+7),'+')
    pos=1;
else
    pos=0;
end

[chr probe_start probe_end]=textread(UCSC_track_file,'%s\t%d\t%d\n','headerlines',1);

probe_mid=(probe_start+probe_end)/2;

N=zeros(1,stp-strt+1);
figure;
for i=1:length(N);
    if pos
        N(i)=length(find(probe_mid<(strt+i-1)));
    else
        N(i)=length(find(probe_mid>(stp-i+1)));
    end
end
plot(1:length(N),N,'bo')

xlabel('Coordinate along gene (bps)','fontsize',21);
ylabel('Number of probes bound on Pol2 at the x coordinate','fontsize',18);
L=abs(stp-strt+1);
W=(1/L)*sum(N)/N(end);
title(['Gene length=' num2str(L) ', W=' num2str(W)],'fontsize',25);
set(gca,'fontsize',18);


        
