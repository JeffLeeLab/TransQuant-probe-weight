function [uid]=find_background_dots(dots,pixel_list,DOT_SIZE_THRESHOLD,fiji_xy)

% the function returns the indices (row numbers) of dots that are not
% removed (above DOT_SIZE_THRESHOLD)

% if no threshold coefficient specified
THRESH_DIST=3;

indout=[];
indin=1:size(dots,1);

% Add cytoplasmic dots that are too large
LL=length(pixel_list{1});
md=zeros(LL,3);
sz=zeros(LL,1);
counter=1;

for j=1:length(pixel_list{1}),
    mat=pixel_list{1}(j).PixelList;
    if size(mat,1)==1,
        md(counter,:)=pixel_list{1}(j).PixelList;
    else
        md(counter,:)=mean(pixel_list{1}(j).PixelList);
    end
    sz(counter)=size(pixel_list{1}(j).PixelList,1);
    counter=counter+1;
end


% To make sure we don't remove TS we will now setdiff the fiji validated
% dots
TS_dot_indices=[];
for i=1:size(fiji_xy,1)
    [min_dist,pt_ind]=min(sqrt((dots(:,1)-fiji_xy(i,1)).^2+(dots(:,2)-fiji_xy(i,2)).^2));
    if min_dist<THRESH_DIST
        TS_dot_indices=[TS_dot_indices pt_ind];
    end
end

uid=find(sz>DOT_SIZE_THRESHOLD);
uid=setdiff(uid,TS_dot_indices);


