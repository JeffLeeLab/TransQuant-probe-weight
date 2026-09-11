function UserData2=refine_dots(UserData,bgid)

L=[];
for i=1:length(UserData.dot_pixel_list)
    L(i)=length(UserData.dot_pixel_list{i});
end
indin=zeros(sum(L),1);
indin(bgid)=1; % remove these

limits=[1 L(1)];
for i=2:length(L),
    limits=[limits;limits(end,2)+1 limits(end,2)+L(i)];
end

for i=1:length(UserData.dot_pixel_list)
    index=intersect(bgid,limits(i,1):limits(i,2));
    index=find(indin(limits(i,1):limits(i,2)));
    UserData.dot_pixel_list{i}(index)=[];
    UserData.dot_pixel_values{i}(index)=[];
    UserData.dot_intensities{i}(index)=[];
end



% pixel_list=UserData.dot_pixel_list;
% for i=1:length(pixel_list),
%     centers_pixel_list{i}=zeros(length(pixel_list{i}),3);
%     for j=1:length(pixel_list{i}),
%         mat=pixel_list{i}(j).PixelList;
%         centers_pixel_list{i}(j,:)=mean(mat);
%     end
% end
% 
% for j=1:length(UserData.dot_pixel_list)
%     indin{j}=[];
% end
% 
% for i=1:size(UserData.dots,1),
%     for j=1:length(UserData.dot_pixel_list)
%         temp=repmat(UserData.dots(i,1:3),size(centers_pixel_list{j},1),1);
%         temp2=centers_pixel_list{j}-temp;
%         [yy,ii]=min(sqrt(sum((temp2'.^2))));
%         if yy<4 & UserData.dots(i,4)==j, % the dot has to be same channel and close enough to the pixel_list mean
%             indin{j}=[indin{j};ii];
%         end
%     end
% end
% 
% for j=1:length(UserData.dot_pixel_list)
%     UserData.dot_pixel_list{j}=UserData.dot_pixel_list{j}(indin{j});
%     UserData.dot_pixel_values{j}=UserData.dot_pixel_values{j}(indin{j});
% end





UserData2=UserData;

