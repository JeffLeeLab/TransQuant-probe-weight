function [ volume_im, THRESH, pixel_list pixel_value pixel_value_raw,SN] = get_dot_3d_image( pixel_list,BORDER,ims )

% for each dot defined by its pixel_list returns the 3D image (cur_ims),
% the THRESH, which is the graythresh level and I, the integrated summed
% intensity of the background subtracted image. Iraw is the sum of the
% non-subtracted values.

DEBUG_PLOTS=0;
I=[];
Iraw=[];

%returns a sub_image with the neighborhood of a dot.
[M,N,P]=size(ims);

x_min_pix=round(max(1,min(pixel_list(:,1))-BORDER));
x_max_pix=round(min(M,max(pixel_list(:,1))+BORDER));
y_min_pix=round(max(1,min(pixel_list(:,2))-BORDER));
y_max_pix=round(min(N,max(pixel_list(:,2))+BORDER));
zmin=min(pixel_list(:,3));
zmax=max(pixel_list(:,3));
volume_im=ims(y_min_pix:y_max_pix,x_min_pix:x_max_pix,zmin:zmax); % IMPORTANT - note the inversion of x and y when we work with images


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Create the thresholded image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

THRESH_vec=zeros(size(volume_im,3),1);
noise_vec=zeros(size(volume_im,3),1);
signal_vec=zeros(size(volume_im,3),1);

newim=zeros(size(volume_im));
bwim=zeros(size(volume_im));

for i=1:size(volume_im,3),
    vec=volume_im(:,:,i);
    vec=vec(:);
    [THRESH,EM] = graythresh(vec/max(vec)); % Choose the threshold for this dot stack
    THRESH=max(vec) * THRESH;
    vec2=vec(vec>THRESH);
    [THRESH,EM] = graythresh(vec2/max(vec2)); % Choose the threshold for the already thresholded image
    THRESH_vec(i)=max(vec2) * THRESH;
    BG=median(vec(vec<THRESH_vec(i)));
    signal_vec(i)=median(vec(vec>THRESH_vec(i)))-BG;
    noise_vec(i)=std(vec(vec<THRESH_vec(i)));
    newim(:,:,i)=max(0,volume_im(:,:,i)-BG);
    
    bwim(:,:,i)=volume_im(:,:,i)>THRESH_vec(i);    
    
end
SN=[median(signal_vec) median(noise_vec)];
bwim=bwlabeln(bwim,6);


% find the relevant connected component
p=pixel_list-repmat([x_min_pix y_min_pix zmin],[size(pixel_list,1) 1])+1;
%index=sub2ind([size(volume_im,1),size(volume_im,2),size(volume_im,3)],p(:,1),p(:,2),p(:,3));
index=sub2ind([size(volume_im,1),size(volume_im,2),size(volume_im,3)],p(:,2),p(:,1),p(:,3)); % note change from previous line
vector=bwim(index);
comp_rel=median(vector(vector>0)); % the relevant connected component


index2=find(bwim==comp_rel);
[i,j,k]=ind2sub(size(volume_im),index2);
%p2=[i j k];
p2=[j i k]; % Note that empirically this is the way it fit the bwim, still need to understand the inversion of axes


% testing - actually use the original p from the LOG counting method (rather than p2), later
% we will use the watershed to refine the dot regions
p2=p;

if ~isempty(p2)
    %pout=refine_by_watershed(p2,p,volume_im); % intersect with the watershed region that includes the original pixels in each plane
    pout=p2; % watershed simply takes too long
else
    pout=p2;
end


% DEBUG plots
if DEBUG_PLOTS
    figure(3);clf
    temp=newim;
    for j=zmin:zmax,
        index=find(pixel_list(:,3)==j);
        index2=find(p2(:,3)==j-zmin+1);
        index3=find(pout(:,3)==j-zmin+1);
        center=[mean(pixel_list(index,1:2),1)]; % the pixel in j z-plane and x,y coordinates same as this dot's centroid
        %get the neighborhood of the centroid in the current z (i.e. the centorid x and y coordinates in the current stack)
        %     x_shift=(BORDER+1)-center(1);
        %     y_shift=(BORDER+1)-center(2);
        pt=pixel_list(index,1:2)-repmat([x_min_pix y_min_pix],[length(index) 1])+1;
        pt=round(pt);
        figure(j-zmin+1)%subplot(3,3,j-zmin+1);
        imagesc(temp(:,:,j-zmin+1));colormap gray
        hold on;
        
        %plot(p2(index2,1),p2(index2,2),'b.');
        plot(pout(index3,1),pout(index3,2),'cs');
        plot(pt(:,1),pt(:,2),'ro');
        title(num2str(j));pause(1);
        axis square;
        hold off;
    end
end

% The final coordinates will be the union of p and p2
%temp=[p;p2];
temp=[pout];
pixel_list=unique(temp,'rows');
sz=[size(volume_im,1) size(volume_im,2) size(volume_im,3)];
if length(sz)<3,
    sz=[size(volume_im) 1];
end
%index=sub2ind(sz,pixel_list(:,1),pixel_list(:,2),pixel_list(:,3));
index=sub2ind(sz,pixel_list(:,2),pixel_list(:,1),pixel_list(:,3)); % CHECK THIS!!!
pixel_value=newim(index);
pixel_value_raw=volume_im(index);
%plot3(pixel_list(:,1),pixel_list(:,2),pixel_list(:,3),'r.');
%title(num2str(sum(pixel_value)));
%pause(1);
pixel_list=pixel_list+repmat([x_min_pix y_min_pix zmin],[size(pixel_list,1) 1])-1;


function [pout,watershed_pts]=refine_by_watershed(pts,root_pts,im)

if size(pts,1)==1,
    pout=pts;
    watershed_pts=pts;
    return;
end
zmin=min([pts]);zmin=zmin(3);
zmax=max([pts]);zmax=zmax(3);

%include=zeros(size(pts,1),1); % 1 if include, 0 otherwise
watershed_pts=[];
h = fspecial('gaussian', 4, 3); % low pass filter to smooth image
for i=zmin:zmax
    cur_ims=im(:,:,i);
    % watershed
    
    %B=imfilter(cur_ims,h);
    B=cur_ims;
    %D = cur_ims;
    D=B;
    D = -D;
    D(D==0) = -Inf;
    L = watershed(D);
    index=find(root_pts(:,3)==i);
%    index=sub2ind([size(cur_ims,1) size(cur_ims,2)],root_pts(index,1),root_pts(index,2)); % CHECK!!!
    index=sub2ind([size(cur_ims,1) size(cur_ims,2)],root_pts(index,2),root_pts(index,1)); % CHECK!!!
    vector=L(index);
        comp_rel=median(double(vector(vector>0))); % the relevant connected component
    index2=find(L==comp_rel);
    [ii,jj]=ind2sub(size(cur_ims),index2);
    watershed_pts=[watershed_pts;jj ii repmat(i,length(ii),1)]; % note the required axes inversion, empirically this is what works
    %figure(44); imagesc(L);
    %temp = cur_ims;
    %temp(L==0) = 10000;
    %figure(55); imagesc(temp);
    
end

pout=intersect(pts,watershed_pts,'rows');



