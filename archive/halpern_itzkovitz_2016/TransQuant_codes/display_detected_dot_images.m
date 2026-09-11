function display_detected_dot_images(UserData,channel,num_per_plot)

% display_detected_dot_images(UserData,channel,num_per_plot)
% Displays the dot intensity profiles on a Z max-projection 



if nargin<2,
    channel=1;
end
if nargin<3,
    num_per_plot=5;
end

figure;clf;
counter=1;
for j=1:length(UserData.dot_pixel_list{channel}),
    mat=UserData.dot_pixel_list{channel}(j).PixelList;
    val=UserData.dot_pixel_values{channel}(j).PixelValues;
    temp=zeros(max(mat(:,1)),max(mat(:,2)),max(mat(:,3)));
    for i=1:size(mat,1), temp(mat(i,1),mat(i,2),mat(i,3))=val(i);end
    mat1=max(temp,[],3);
    [indexi,indexj]=find(mat1>0);
    mat1=mat1(min(indexi):max(indexi),min(indexj):max(indexj));
    mat2=max(temp,[],2);
    [indexi,indexj]=find(mat2>0);
    mat2=mat2(min(indexi):max(indexi),min(indexj):max(indexj));
    mat3=max(temp,[],1);
    [indexi,indexj]=find(mat3>0);
    mat3=mat3(min(indexi):max(indexi),min(indexj):max(indexj));
    subplot(num_per_plot,num_per_plot,counter);imagesc(mat1);axis square;colormap gray;
    counter=counter+1;
    title(['Z projection of ' num2str(j) ', Size=' num2str(length(val))]);
%     subplot(10,2,2*mod(j,10));imagesc(mat2);axis square;colormap gray;
%     title(['X projection of ' num2str(j)]);
%     %     subplot(1,3,3);imagesc(mat3);axis square;colormap gray;
    %     title(['Y projection of ' num2str(j)]);
    if mod(j,num_per_plot^2)==0
        pause;
        counter=1;
    end
end