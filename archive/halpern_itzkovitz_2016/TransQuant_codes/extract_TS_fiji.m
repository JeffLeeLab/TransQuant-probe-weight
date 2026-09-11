function [UserData]=extract_TS_fiji(UserData)

% For each cell adds a cell that includes the position and intensity of its
% transcription sites
% The added info is in cell.TS_list - its a n*11 matrix with n TS per cell
% with the following column - x_intron y_intron z_intron I_intron A_intron x_exon
% y_exon z_exon I_exon A_exon nuc_number (1/2)
% The main thing we need to still correct here is to refine the pixel_list
% of each dot and the pixel value by using graythresh


THRESH_DIST=5; % the fiji_xy must be within 3 pixels of the xy of the TS dot to be considered

VELOCITY=34; % We assume 34 bp/sec as measured in Bahar Halpern et al.

if ~isfield(UserData,'L')
    L=6000;
else
    L=UserData.L;
end
answer = inputdlg('Enter gene length','Gene length dialog',1,{num2str(L)});

if isempty(answer)
    return;
end
GENE_LENGTH=str2num(answer{1});

if ~isfield(UserData,'W')
    W=0.5'
else
    W=UserData.W;
end
answer = inputdlg('Enter probe weighting factpr','Probe weight dialog',1,{num2str(W)});
if isempty(answer)
    return;
end
PROBE_WEIGHT_FACTOR=str2num(answer{1});


if ~isfield(UserData,'dot_pixel_values'),
    return;
end


num_exons=length(UserData.dot_pixel_list{1});

c=[];

c=[];
for i=1:num_exons,
    temp=UserData.dot_pixel_list{1}(i).PixelList;
    if size(temp,1)==1, temp=[temp; temp];end
    c=[c;mean(temp,1)];
    A_exons(i)=size(temp,1);
end
c_exons=c;


% compute exon integrated intensity
exon_I=zeros(num_exons,1);

for i=1:num_exons,
    exon_I(i)=sum_max_stack(UserData.dot_pixel_list{1}(i).PixelList,UserData.dot_pixel_values{1}(i).PixelValues);
end


% Extract the dots that are in segmented regions
if ~isfield(UserData,'cell')
    indin=1:size(UserData.dots,1);
elseif isempty(UserData.cell) 
    indin=1:size(UserData.dots,1);
else
    indin=find(~isnan(UserData.dots(:,end)));
end
exon_I=exon_I(indin);
c_exons=c_exons(indin,:);
dots=UserData.dots(indin,:);
Xst=size(c_exons,1);

%UserData.fiji_xy=c_exons(:,1:2); % This was a debugging line to ensure
%median occupancy of a mRNA dot is 1/PROBE_WEIGHT_FACTOR

UserData.TS=NaN*ones(size(UserData.fiji_xy,1),6); % for every TS we have x,y,z, intensity, M and mu
for i=1:size(UserData.fiji_xy,1)
    [min_dist,pt_ind]=min(sqrt((dots(:,1)-UserData.fiji_xy(i,1)).^2+(dots(:,2)-UserData.fiji_xy(i,2)).^2));
    if min_dist<THRESH_DIST
        % compute the estimated M
        indZ=find((c_exons(:,3)<(ceil(dots(pt_ind,3)+1)) & (c_exons(:,3)>=(floor(dots(pt_ind,3))-1))));
        exon_I_Z=median(exon_I(indZ));
        M=exon_I(pt_ind)./(PROBE_WEIGHT_FACTOR*exon_I_Z);
        UserData.TS(i,:)=[dots(pt_ind,1:3) exon_I(pt_ind) M M*3600*VELOCITY/GENE_LENGTH];
    end
end
indout=find(any(isnan(UserData.TS')));
indin=setdiff(1:size(UserData.TS,1),indout);
UserData.TS=UserData.TS(indin,:);

beta=sum(UserData.TS(:,end));
Mtot=sum(UserData.TS(:,end-1));
delta=beta/Xst;    % units of hr^-1
UserData.Mtot=Mtot;
UserData.delta=delta;
UserData.beta=beta;
UserData.Xst=Xst;
str=['Total Transcription rate =' num2str(beta,'%.2f') ' [mRNA/hr], Degradation rate=' num2str(delta,'%.2f') ', mRNA lifetime=' num2str(log(2)/delta,'%.2f') ' hr, Total #mRNA=',num2str(Xst,'%d')];
msgbox(str,'modal');




