function ims=parse_stack(filename,first,last,FILTER_OUTLIER)

% ims=parse_stack(filename)
% 12/11/09
% Parses metamorph stack file with S 2D images of size n*m and outputs a
% n*m*S double 3D array

if nargin<4,
    FILTER_OUTLIER=0;
end
ZTHRESH=5;

info=imfinfo(filename);
if last>length(info),
    last=length(info);
end

if nargin < 2
    stack = imread(filename);
else
    stack=zeros(1024,1024,last-first+1);
    %stack=uint16(stack);
    for i=first:last,
        temp = imread(filename,i,'info',info);
        stack(1:size(temp,1),1:size(temp,2),(i-first+1))=temp;
        
    end
end

ims=stack;

% [n,m]=size(stack(1).data);
% ims=zeros(n,m,length(stack));
%
% for i=1:length(stack),
%     %ims(:,:,i)=double(stack(i).data);
%     ims(:,:,i)=stack(i).data;
% end
%
% if FILTER_OUTLIER,
%     ind=find((ims-mean(ims(:)))/std(ims(:))>ZTHRESH);
%     if ~isempty(ind),
%         index=min(find((ims-mean(ims(:)))/std(ims(:))>ZTHRESH));
%         ims(ind)=ims(index);
%     end
% end