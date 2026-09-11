function [dot_area,dot_intensity]=compute_dot_features(imt,mid)

% [dot_area,dot_intensity]=compute dot_features(imt);
% Receives a 2D image which contains the vicinity of a dot. starts at the
% dot center and dilates the image recomputing the intensity of sorrounding
% rings. Stops when ring intensity starts to increas

%DOT_AREA=floor(length(imt)/2);

% compute background
[y,ord]=sort(imt(:));
BG=y(floor(length(y)/10));

imt=-imt;
% save_pixel=imt(mid(1),mid(2));
imt(1,:)=-inf;imt(end,:)=-inf;imt(:,1)=-inf;imt(:,end)=-inf;
% imt(mid(1),mid(2))=save_pixel;
L=watershed(imt);
[x,y]=find(L==L(mid(1),mid(2)));
ind=find(L==L(mid(1),mid(2)));
dot_area=length(ind);
dot_intensity=sum(-imt(ind))-dot_area*BG;



M=-median(-imt(:));
imt(1,:)=M;imt(end,:)=M;imt(:,1)=M;imt(:,end)=M;
imagesc(-imt);
hold on;
plot(y,x,'rx');
hold on;
plot(mid(1),mid(2),'go');
title(['Area=' num2str(dot_area) ', Intensity=' num2str(dot_intensity)]);

