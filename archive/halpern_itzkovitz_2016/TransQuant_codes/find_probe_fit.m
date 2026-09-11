function [missmatch,pos]=find_probe_fit(seq,probe_set)
 
% Return the mismatch distribution of a probe library and a given
% sequence.
% [probe2,nm2,ind2,pos2,gc]=textread('mouse_ephB2_probes.txt','%s%s%d%d%f','delimiter',',')
 
pos=zeros(length(probe_set),1);
missmatch=zeros(length(probe_set),1);
for i=1:length(probe_set),
    str=probe_set{i};
    str(str==' ')='';
    %temp=lower(seqrcomplement(probe_set{i}));
    temp=lower(seqrcomplement(str));
    fit=zeros(length(seq)-length(temp)+1,1);
    for j=1:length(seq)-length(temp)+1,
       fit(j)=length(find(lower(seq(j:j+length(temp)-1))==temp));
    end
    [yy,ii]=max(fit);
    missmatch(i)=length(temp)-yy;
    pos(i)=ii;
end
