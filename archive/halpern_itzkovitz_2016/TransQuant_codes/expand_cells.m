function UserData=expand_cells(UserData,MARGIN,plot_all)

% cells=expand_cells(UserData)
% Expands the cells by a certain MARGIN and makes sure they do not
% intersect
%load('Z:\users\shalevi\Langerhans_Islets\Donatello_November_30_2010\D11_Ins2b_GCG_20_percent_FD\D11_Ins2b_GCG_20_percent_FD.001.ima','-mat');


if nargin<2,
    %MARGIN=10; % how many pixels to exapand along the periphery
    MARGIN=UserData.MARGIN; % how many pixels to exapand along the periphery
end

if nargin<3,
    plot_all=0;
end

% Plot the cells

if plot_all,
hold on;
for i=1:length(UserData.cell),
    
    plot(UserData.cell(i).edge(:,1),UserData.cell(i).edge(:,2));
    text(UserData.cell(i).center(1),UserData.cell(i).center(2),num2str(i));
end
end

for i=1:length(UserData.cell),
    num_edges=size(UserData.cell(i).edge,1);
    vec=UserData.cell(i).edge; % Original edges
    vec_shift=vec-repmat(UserData.cell(i).center,num_edges,1);
    D=sqrt(vec_shift(:,1).^2+vec_shift(:,2).^2);
    FACTOR=(D+MARGIN)./D;
    FACTOR_IN=(D-MARGIN)./D;
    if FACTOR_IN<0,
        FACTOR_IN=0.1*ones(length(FACTOR_IN),1);
    end
    % create an expanded rim version for the calculation of the rim intensities
    UserData.cell(i).extended_edge=repmat(UserData.cell(i).center,num_edges,1)+vec_shift.*repmat(FACTOR,1,size(vec_shift,2));
    % also create a trimmed in version for the calculation of the nuclear intensities
    UserData.cell(i).internal_edge=repmat(UserData.cell(i).center,num_edges,1)+vec_shift.*repmat(FACTOR_IN,1,size(vec_shift,2));
    if plot_all,
    plot(UserData.cell(i).extended_edge(:,1),UserData.cell(i).extended_edge(:,2),'r');
    end
end

% Smooth the polygons
% INTERP_FACTOR=4;
% for i=1:length(UserData.cell),
%     x1=UserData.cell(i).extended_edge(:,1);
%     y1=UserData.cell(i).extended_edge(:,2);
%     q = 0:length(x1)-1;
%     qq = 0:1/INTERP_FACTOR:length(x1)-1;
%     xpts_spline = spline(q,x1,qq);
%     ypts_spline = spline(q,y1,qq);
%     %SplineCurve = plot(xpts_spline,ypts_spline,'r');
%     UserData.cell(i).extended_edge=[xpts_spline' ypts_spline'];
% end

% Now correct each cell by subtracting from it each of the other
% extended cells
display('Subtracting overlapping cells');
for i=1:length(UserData.cell),
    UserData.cell(i).extended_edge_new=UserData.cell(i).extended_edge;
end
for i=1:length(UserData.cell),
    %display(['Cell #' num2str(i) 'of ' num2str(length(UserData.cell))]);
    UserData.cell(i).extended_edge_new=UserData.cell(i).extended_edge;
    for j=setdiff(1:length(UserData.cell),i),
        x1=UserData.cell(i).extended_edge_new(:,1);
        y1=UserData.cell(i).extended_edge_new(:,2);
        [x1,y1]=poly2cw(x1,y1);
        x2=UserData.cell(j).extended_edge_new(:,1);
        y2=UserData.cell(j).extended_edge_new(:,2);
        [x2,y2]=poly2cw(x2,y2);
        x3=UserData.cell(j).edge(:,1);
        y3=UserData.cell(j).edge(:,2);
        [x3,y3]=poly2cw(x3,y3);
        [x4 y4]=polybool('union', x2, y2, x3, y3);
        [x4,y4]=poly2cw(x4,y4);
        [subtx subty]=polybool('subtraction', x1, y1, x4, y4);
        if ~isempty(subtx)
            UserData.cell(i).extended_edge_new=[subtx subty];
        end
    end
    ind=find(isnan(UserData.cell(i).extended_edge_new(:,1)));
    UserData.cell(i).extended_edge_new=UserData.cell(i).extended_edge_new(setdiff(1:size(UserData.cell(i).extended_edge_new,1),ind),:);
    UserData.cell(i).extended_edge=UserData.cell(i).extended_edge_new;
    if plot_all,
    plot(UserData.cell(i).extended_edge(:,1),UserData.cell(i).extended_edge(:,2),'g');
    end
end



% Subtract the nucleus curve from the subtracted extended polygons
%figure;
for i=1:length(UserData.cell),
    x1=UserData.cell(i).extended_edge_new(:,1);
    y1=UserData.cell(i).extended_edge_new(:,2);
    [x1,y1]=poly2cw(x1,y1);
    x2=UserData.cell(i).edge(:,1);
    y2=UserData.cell(i).edge(:,2);
    [x2,y2]=poly2cw(x2,y2);
    [subtx subty]=polybool('subtraction', x1, y1, x2, y2);
    if ~isempty(subtx)
        UserData.cell(i).rim_edge=[subtx subty];
    end
    ind=find(isnan(UserData.cell(i).rim_edge(:,1)));
    UserData.cell(i).rim_edge=UserData.cell(i).rim_edge(setdiff(1:size(UserData.cell(i).rim_edge,1),ind),:);
    %plot(UserData.cell(i).extended_edge(:,1),UserData.cell(i).extended_edge(:,2),'k');
    %plot(UserData.cell(i).internal_edge(:,1),UserData.cell(i).internal_edge(:,2),'c');
end

% Plot all cells
if plot_all,
    cmap=gray(10);
    figure;
    hold on;
    for i=1:length(UserData.cell),
        ind=randperm(10);ind=ind(1);
        patch(UserData.cell(i).extended_edge(:,1),UserData.cell(i).extended_edge(:,2),cmap(ind));
        title(num2str(i));
        %pause;
    end
    axis square;
end


