function [num_dots,locations,threshold,pixel_list,pixel_values,SN]=threshold_selection(handles,thresholds,thresholdfn,ims,ims2,BW,width,offset,proj_method,enhance_method,ims3,USE_LOG,DOT_SIZE_THRESHOLD,SHARP_FACTOR,FAST_MODE)

% [num_dots,locations]=process_file_threshold_selection(thresholds,thresholdfn,ims,ims2,BW)
% 20/11/09 Shalev
% Allows the user to choose the appropriate threshold. Receives the cropped
% original image - ims, and the LOG filtered image - ims2.


if nargin<15
    FAST_MODE=0;
end

if ~isempty(BW)
    USE_ERODE_MASK=1;
    n_stack=size(ims2,3);
    ERODE_MASK_SIZE=20;
    
    if USE_ERODE_MASK,
        BW2=imerode(BW,strel('square',ERODE_MASK_SIZE));
        for kk=1:n_stack,
            ims2(:,:,kk)=ims2(:,:,kk).*BW2;
        end
    end
    for kk=1:n_stack,
        temp=ims2(:,:,kk);
        temp(temp==0)=min(min(temp(temp~=0)));
        ims2(:,:,kk)=temp;
    end
end

% Now ask the user to select a threshold and animate the resulting image
[t nc cv]= auto_thresholding(thresholds,thresholdfn,width,offset);
%t2=graythresh(ims2);
y = nc;
x=thresholds(t);
bwl = ims2 > x;
[lab,num_dots] = bwlabeln(bwl);
locations=regionprops(lab,'centroid'); % The really important matlab function for getting info on connected components
for i=1:length(locations),coor(i,:)=locations(i).Centroid;end
locations=coor;
pixel_list=regionprops(lab,'PixelList');
%pixel_values=regionprops(lab,ims,'PixelValues');
%Aviezer-option to get the values by the log filtered image

if USE_LOG==1
    pixel_values=regionprops(lab,ims3,'PixelValues');
else
    pixel_values=regionprops(lab,ims,'PixelValues');
end

screen_size = get(0, 'ScreenSize');
figure;
set(gcf, 'Position', [0 0 screen_size(3) screen_size(4) ] );

%datacursormode on

h1=subplot('Position',[0.1 0.55 0.25 0.35]);
plot(thresholds,thresholdfn);ylabel('Number of dots');xlabel('Threshold');
%axis([0 inf 0 1000])
title('automatic thresholding')
line([t/100 t/100],ylim,'Color','r');
hold on;
hxy = plot(x,y,'r+');
htip = text(x,y,[' (' num2str(x,'%10.2f') ',' num2str(y) ')']);
%plot(t2,interp1(thresholds,cv,t2),'go');
hold off;

h2=subplot('Position',[0.1 0.1 0.25 0.35]);
plot(thresholds,cv);
hold on;
plot(t/100,cv(t),'r*');
%plot(t2,interp1(thresholds,cv,t2),'go');
xlabel(['default threshold=' num2str(t/100) ', ' 'num. of dots=' num2str(nc)])
ylabel('inverse of coefficient variation')
hold off

im0=zprojection(ims,proj_method);
im_log0=zprojection(ims3,proj_method);
%im0=enhance_image(im0,enhance_method);

h3=subplot('Position',[0.4 0.1 0.5 0.8]);
%imagesc((uint16(SHARP_FACTOR*im_log0+im0)));axis square;colormap gray;axis off;
imagesc(imadjust(uint16(SHARP_FACTOR*im_log0+im0)));axis square;colormap gray;axis off;
title('click at appropriate x/threshold value at top left figure, or press any key to exit');drawnow;
hold on;
for i=1:2,
    h4=plot(locations(:,1),locations(:,2),'o','MarkerSize',4);
    pause(0.2);
    delete(h4);
    pause(0.2);
end
h4=plot(locations(:,1),locations(:,2),'o','MarkerSize',4);
hold off

go_on=1; % As the long as the user is not satisfied
while go_on
    %clear locations bwl SelectionType CurrentPoint
    UserInput=2;
    while UserInput == 2
        UserInput = waitforbuttonpress;      % Wait for click
    end
    
    if UserInput == 1
        
        %        selection = questdlg('Happy with this?','Happy with this?','Yes','No','Yes');
        %        if strcmp(selection,'Yes')
        go_on=0;
        title('closing this window ...');drawnow;
        close
        figure(handles.figure1);
        
        %title('Start to process next image...','FontSize',18,'Color','r');drawnow;
        threshold=thresholds(round(x*100));
        bwl = ims2 > threshold;
        [lab,num_dots] = bwlabeln(bwl);
        locations=regionprops(lab,'centroid');
    else
        subplot(h1);
        %         UserInput=1;
        %         while UserInput > 0
        %             UserInput = waitforbuttonpress;      % Wait for click
        %         end
        CurrentPoint = get(gca,'CurrentPoint');
        x = CurrentPoint(1);
        ndots = thresholdfn(round(x*100));
        subplot(h1)
        hold on;
        delete(hxy)
        hxy = plot(x,ndots,'r+');
        hold off;
        delete(htip)
        htip=text(x,ndots,[' (' num2str(x,'%10.2f') ',' num2str(ndots) ')']);
        
        x=thresholds(round(x*100));
        bwl = ims2 > x;
        
        [lab,num_dots] = bwlabeln(bwl);
        locations=regionprops(lab,'centroid');
        pixel_list=regionprops(lab,'PixelList');
        %Also here = option to change to ims3
        pixel_values=regionprops(lab,ims,'PixelValues');
        clear coor;
        coor=[];
        for i=1:length(locations),coor(i,:)=locations(i).Centroid;end
        locations=coor;
        
        subplot(h3);
        
        %imagesc((uint16(SHARP_FACTOR*im_log0+im0)));axis square;colormap gray;axis off;
        imagesc(imadjust(uint16(SHARP_FACTOR*im_log0+im0)));axis square;colormap gray;axis off;
        title('click at appropriate x/threshold value at top left figure, or press any key to exit');drawnow;
        hold on;
        %   for i=1:2,
        subplot(h3);
        h4=plot(locations(:,1),locations(:,2),'o','MarkerSize',4);
        pause(0.2);
        delete(h4);
        pause(1);
        %   end
        h4=plot(locations(:,1),locations(:,2),'o','MarkerSize',4);
        hold off
        
    end
    %SelectionType = get(FigureHandle,'SelectionType');         % Get information about the last button press
end

locations_non_filtered=locations;

% filter out large dots
%DOT_SIZE_THRESHOLD=200;
for i=1:length(pixel_list), L(i)=size(pixel_list(i).PixelList,1);end
% Note we removed this because we lose the large exon dots in the nucleus
%indin=find(L<DOT_SIZE_THRESHOLD);   
indin=1:size(locations,1);
num_dots=length(indin);
pixel_list=pixel_list(indin);
pixel_values=pixel_values(indin);
locations=locations(indin);

clear coor;
threshold=x;
if num_dots
    for i=1:length(locations),coor(i,:)=locations(i).Centroid;end
    locations=coor;
end


SN=zeros(length(pixel_list),2); % will hold the signal and the noise
if ~FAST_MODE    
    % Get refined pixel list    
    BORDER=10;
    for k=1:length(pixel_list),
        pixel_list(k).PixelList=expand_pixel_list_to_minimal_size(pixel_list(k).PixelList,pixel_values(k).PixelValues);
        k;
        if mod(k,length(pixel_list)/10)==0
            [k length(pixel_list)]
        end
        
        [ volume_im, THRESH, pixel_list2(k).PixelList pixel_values2(k).PixelValues pixel_values_raw,SN(k,:)]  = get_dot_3d_image( pixel_list(k).PixelList,BORDER,ims );
    end
    
    pixel_list=pixel_list2;
    pixel_values=pixel_values2;
    
    
end