function [new_pixel_list]=expand_pixel_list_to_minimal_size(pixel_list,pixel_values)

% [new_pixel_list]=expand_pixel_list_to_minimal_size(pixel_list)
% Makes sure that in the stack in which there are most pixels there is a
% least 3*3 square

new_pixel_list=pixel_list;
% find the Z section with the maximal intensity
%pu=unique(pixel_list(:,3));
pu=min(pixel_list(:,3)):max(pixel_list(:,3));
%pu2=min(pixel_list(:,3)):max(pixel_list(:,3));
% if ~isempty(setdiff(pu,pu)) | ~isempty(setdiff(pu,pu))
%     type('ERROR');
% end
s=zeros(length(pu),1);
L=zeros(length(pu),1);
for i=1:length(pu),
    index=find(pixel_list(:,3)==pu(i));
    s(i)=sum(pixel_values(index));
    L(i)=length(pixel_values(index));
end
[yy,ii]=max(s);
if L(ii)<9,
    index=find(pixel_list(:,3)==pu(ii));
    new_pixel_list(index,:)=[];
    if length(index)==1
        center=round((pixel_list(index,:)));
    else
        center=round(mean(pixel_list(index,:)));
    end
    % make sure we are not right at the edge
    center(1)=max(center(1),2);
    center(1)=min(1023,center(1));
    center(2)=max(center(2),2);
    center(2)=min(1023,center(2));
    [x,y]=meshgrid(center(1)-1:center(1)+1,center(2)-1:center(2)+1);
    new_pixel_list=[new_pixel_list;x(:) y(:) repmat(pu(ii),length(x(:)),1)];
end