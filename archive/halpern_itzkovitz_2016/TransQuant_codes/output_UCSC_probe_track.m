function output_UCSC_probe_track(UCSC_fasta_file,probe_file, track_file)

% produces a track which can be loaded into UCSC

s=textread(probe_file,'%s');
[a,b]=fastaread(UCSC_fasta_file);
% Extract the chromosome and location from the fasta header
ind1=findstr(a,'chr');
ind2=findstr(a,':');
ind3=findstr(a,'-');
ind4=findstr(a,'5''pad');
chr=a(ind1:(ind2-1));
pos_start=str2num(a((ind2+1):(ind3-1)));
pos_end=str2num(a((ind3+1):(ind4-2)));



% find probe locations
[missmatch,pos]=find_probe_fit(b,s);


% output track file
fid=fopen(track_file,'w');
fprintf(fid,'track name=%s\n',track_file);

% Determine if the strand is + or -
ind5=findstr(a,'strand=');
if a(ind5+7)=='+',
    for i=1:length(pos),
        fprintf(fid,'%s\t%s\t%s\n',chr, num2str(pos_start-1+pos(i)), num2str(pos_start-1+pos(i)+19));
    end;
else % negative strand
    for i=1:length(pos),
        fprintf(fid,'%s\t%s\t%s\n',chr, num2str(pos_end-pos(i)-19), num2str(pos_end-pos(i)+1));
    end;
end
%close all
fclose all
