function channel=all_channel_names(dirname,file_index)
% find files with specified file_index in the directory dirname.

d=dir([dirname '*.tif']);
if length(d)==0,
    d=dir([dirname '*.stk']);
end

if length(d)==0,
    d=dir([dirname '*.TIF']);
end

channel={};
for i=1:length(d),
    %index = d(i).name((end-6):(end-4));
    index = findstr(d(i).name,'_');
    %if strcmp(file_index,index)
    %    channel=[channel d(i).name(1:end-7)];
    %end
    channel=[channel d(i).name(1:(index-1))];
end

channel=unique(channel);