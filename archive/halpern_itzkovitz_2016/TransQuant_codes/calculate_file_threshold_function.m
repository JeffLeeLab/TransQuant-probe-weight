function [thresholdfn,thresholds,ims2,ims3]=calculate_file_threshold_function(handles,ims,BW,N,sigma,n_thresholds,ZTHRESH)

% [thresholdfn,thresholds]=calculate_threshold_function(ims,mat_name,N,sigma)
% The function applies Arjun's method for dot counting
% The function first applies a laplacian filter with size N and std sigma
% to the original image. It then normalizes and counts connected components
% on a binary image obtained with different thresholds.The function returns
% the threshold function, and saves the function as well as filtered and
% non-filtered images in a matlab with a specified name.



if nargin<4,
    sigma=1.5;
end
if nargin<3,
    N=15;
end

set(handles.text1,'String','Applying Laplacian filter ...');
ims2 = LOG_filter(ims,N,sigma);
ims3=ims2;

% October 6, 2014 - TRY THIS AND THEN REMOVE IF NECESSARY
% This part saturates all really bright pixels to improve the dynamic range
% in the threshold selection
if ZTHRESH>0
    %ZTHRESH=2;
    Zim=(ims3-mean(ims3(:)))/std(ims3(:));
    ind=find(Zim>ZTHRESH);
    if ~isempty(ind),
        index=min(find(Zim>ZTHRESH));
        ims3(ind)=ims3(index);
    end
    ims2=ims3;
end
%%%%%% REMOVE UNTIL HERE

% Normalize ims2
ims2 = ims2/max(ims2(:));

set(handles.text1,'String','Start computing intensity threshold ...');
% This function call will find the number of mRNAs for all thresholds
% Note BW is used to mask the outline, where abberant dots appear
[thresholdfn,thresholds] = BW_multithreshstack(handles,ims2,BW,n_thresholds);

% These are the thresholds
%thresholds = (1:100)/100;

%mat_name(findstr(mat_name,' '))='_';
%save(mat_name,'thresholds','thresholdfn','ims','ims2');
