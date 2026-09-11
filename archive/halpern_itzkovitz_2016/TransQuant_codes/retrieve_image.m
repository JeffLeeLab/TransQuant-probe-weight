function [im1a im2a im3a]=retrieve_image(im1,im2,im3,coor,stack_range,RANGE);

% [im1a im2a im3a]=retrieve_image(im1,im2,im3,coor,RANGE);
% Retrieve a subimage from each which is centered on the [x,y,z] in coor.
% RANGE is [XRANGE YRANGE ZRANGE]

if nargin<6,
    RANGE=[10 10 2];
end

im1=im1(:,:,stack_range(1):stack_range(2));
im2=im2(:,:,stack_range(1):stack_range(2));
im3=im3(:,:,stack_range(1):stack_range(2));

coor=round(coor);
YLOW=max(1,coor(1)-RANGE(1));
YHIGH=min(size(im1,1),coor(1)+RANGE(1));

XLOW=max(1,coor(2)-RANGE(2));
XHIGH=min(size(im1,2),coor(2)+RANGE(2));

ZLOW=max(1,coor(3)-RANGE(3));
ZHIGH=min(size(im1,3),coor(3)+RANGE(3));

im1a=im1(XLOW:XHIGH,YLOW:YHIGH,ZLOW:ZHIGH);im1aF=im1(XLOW:XHIGH,YLOW:YHIGH,coor(3));
im2a=im2(XLOW:XHIGH,YLOW:YHIGH,ZLOW:ZHIGH);im2aF=im2(XLOW:XHIGH,YLOW:YHIGH,coor(3));
im3a=im3(XLOW:XHIGH,YLOW:YHIGH,ZLOW:ZHIGH);im3aF=im3(XLOW:XHIGH,YLOW:YHIGH,coor(3));

figure;
subplot(1,3,1);imagesc(im1aF');axis square; colormap gray;axis off;
subplot(1,3,2);imagesc(im2aF');axis square; colormap gray;axis off;
subplot(1,3,3);imagesc(im3aF');axis square; colormap gray;axis off;

figure;
% Now add the relevant plots
% Automatically detect background - the 5 percentile low of the histograms
%[p,x]=hist(im1(:),1000);p=p/sum(p);CP=cumsum(p);plot(x,CP);BG1=x(min(find(CP>0.05)));
BG1=mean(im1(:));
%[p,x]=hist(im2(:),1000);p=p/sum(p);CP=cumsum(p);plot(x,CP);BG2=x(min(find(CP>0.05)));
BG2=mean(im2(:));
%[p,x]=hist(im3(:),1000);p=p/sum(p);CP=cumsum(p);plot(x,CP);BG3=x(min(find(CP>0.05)));
BG3=mean(im3(:));
%BG1=0;BG2=0;BG3=0;
vec1=[im1(XLOW:XHIGH,YLOW+round((YHIGH-YLOW)/2),ZLOW:ZHIGH)-BG1];
vec2=[im2(XLOW:XHIGH,YLOW+round((YHIGH-YLOW)/2),ZLOW:ZHIGH)-BG2];
vec3=[im3(XLOW:XHIGH,YLOW+round((YHIGH-YLOW)/2),ZLOW:ZHIGH)-BG3];
mny=min([vec1(:);vec2(:);vec3(:)]);
mxy=max([vec1(:);vec2(:);vec3(:)]);
for i=1:(ZHIGH-ZLOW+1),
    subplot(5,3,3*(i-1)+1);
    colr='b';
    %if i>(ZHIGH-ZLOW)/2 & i<(ZHIGH-ZLOW)/2+1,
    if i==RANGE(end)+1,
    %if i==coor(3),
        colr='g';
    end
    plot(im1(XLOW:XHIGH,YLOW+round((YHIGH-YLOW)/2),ZLOW+i-1)-BG1,colr,'LineWidth',2);axis tight;set(gca,'XTick',[]);axis([xlim mny mxy]);      
end

for i=1:(ZHIGH-ZLOW+1),
    subplot(5,3,3*(i-1)+2);
    colr='b';
    %if i>(ZHIGH-ZLOW)/2 & i<(ZHIGH-ZLOW)/2+1,
    if i==RANGE(end)+1,
    %if i==coor(3),
        colr='g';
    end
    plot(im2(XLOW:XHIGH,YLOW+round((YHIGH-YLOW)/2),ZLOW+i-1)-BG2,colr,'LineWidth',2); axis tight;set(gca,'XTick',[]);axis([xlim mny mxy]);        
end

for i=1:(ZHIGH-ZLOW+1),
    subplot(5,3,3*(i-1)+3);
    colr='b';
    %if i>(ZHIGH-ZLOW)/2 & i<(ZHIGH-ZLOW)/2+1,
    if i==RANGE(end)+1,
    %if i==coor(3),
        colr='g';
    end
    plot(im3(XLOW:XHIGH,YLOW+round((YHIGH-YLOW)/2),ZLOW+i-1)-BG3,colr,'LineWidth',2);axis tight;set(gca,'XTick',[]);axis([xlim mny mxy]);  
end



