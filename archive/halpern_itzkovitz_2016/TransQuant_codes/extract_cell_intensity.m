function vec=extract_cell_intensity(UserData,ims)

% iterates over all cells and for each cell outputs the average pixel
% intensity of the relevant image

if ~isfield(UserData,'cell')
    vec=[];
    return;
end
if isempty(UserData.cell)
    vec=[];
    return;
end

% for i=1:length(UserData.cell)
%     mask=poly2mask(UserData.cell(i).edge(:,1),UserData.cell(i).edge(:,2),size(ims,1),size(ims,2));
%     index=find(mask);
%     m=zeros(size(ims,3),1);
%     for j=1:size(ims,3)
%         im_temp=ims(:,:,j);
%         % find the background
%         [yy,ord]=sort(im_temp(:));
%         BG=yy(round(length(yy)/10));
%         m(j)=mean(im_temp(index))-BG;
%     end
%     vec(i)=mean(m);
% end


% Changed on 12/8/2015
% We will obtain the background for each Z-stack as the median intensity of
% all nuclear voxels. Then we will average all background-subtracted CYTOPLASMIC pixels for each
% cell

BG_mat=NaN*ones(length(UserData.cell),size(ims,3)); % The BG based on the nucleus of each cell
mask2=zeros(size(ims,1),size(ims,2));
for i=1:length(UserData.cell)
    for j=1:size(ims,3)
        %mask1=poly2mask(UserData.cell(i).edge(:,1),UserData.cell(i).edge(:,2),size(ims,1),size(ims,2)); % the cell
        if ~isempty(UserData.nuc_contour{i}{1})
            mask2=mask2+poly2mask(UserData.nuc_contour{i}{1}(:,1),UserData.nuc_contour{i}{1}(:,2),size(ims,1),size(ims,2)); % the nucleus
            bg_mask=poly2mask(UserData.nuc_contour{i}{1}(:,1),UserData.nuc_contour{i}{1}(:,2),size(ims,1),size(ims,2)); % the nucleus
            index=find(bg_mask);
            im_temp=ims(:,:,j);
            BG_mat(i,j)=median(im_temp(index));
        end
    end
end
mask2=mask2>0;

% BG_mat has NaNs for every cell that did not have a nucleus, replace with
% the median nuclear bg intensity in that Z-stack
vec_nan_median=NaN*ones(1,size(BG_mat,2));
for j=1:size(BG_mat,2)
   vec_nan_median(j)=nanmedian(BG_mat(:,j));
end
for i=1:size(BG_mat,1)
    for j=1:size(BG_mat,2)
        if isnan(BG_mat(i,j))
            BG_mat(i,j)=vec_nan_median(j);
        end
    end
end
        

% now to make sure the background does not include "nuclei" which are
% actually not nuclei we will take for each Z-stack the median of the lower
% 50% nuclei
BG_vec=zeros(1,size(ims,3));
for i=1:size(ims,3)
    [yy,ord]=sort(BG_mat(:,i));
    indin=ord(1:floor(length(ord)/2)); % I take the bottom 50%
    BG_vec(i)=median(BG_mat(indin,i));
end

% now subtract the appropriate background from each stack
m=zeros(length(UserData.cell),size(ims,3)); % The BG corrected intensities for each cell
%index_BG=find(mask2);
for j=1:size(ims,3)
    im_temp=ims(:,:,j);
    %BG=median(im_temp(index_BG));
    for k=1:length(UserData.cell)
        mask=poly2mask(UserData.cell(k).edge(:,1),UserData.cell(k).edge(:,2),size(ims,1),size(ims,2));
        mask=mask-mask2;
        index=find(mask==1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% NOTE: ON 20/9/2015 we removed the BG since in zonation profiles
        % if half of the picture has lower expression this artificially
        % increases the intentisty of the other half by substracting a
        % lower BG.
        %m(k,j)=median(im_temp(index))-BG_vec(j);
        m(k,j)=median(im_temp(index));
    end    
end
vec=mean(m,2);
vec(vec<0)=0;



