
function [distance,Dt]=find_dot_cell_dist(UserData,MAX_MARGIN,MAX_DIST_TO_CONSIDER)

% expands the marked cells by MAX_MARGIN and find the distance to the
% nucleus for all dots within cells
% MAX_DIST_TO_CONSIDER - the maximal distance from the nuclear periphery
% which we would want to analyze
% Inputs: UserData - the UserData struct
%         MAX_MARGIN - How much to expand the nuclear contours
%         MAX_DIST_TO_CONSIDER - beyond this distance a distance of a dot to its
%                nucleus is not considered.
% Outputs: distance - N*1 vector for the N dots for their distance from the
%          closest nucleus.
%          Dt - distances to nucleus stratified by channel          

if nargin<3,
    MAX_DIST_TO_CONSIDER=10; % 10 microns
end
if nargin<2,
    MAX_MARGIN=5; % 5 microns
end

MAX_MARGIN=round(MAX_MARGIN/0.13);
MAX_DIST_TO_CONSIDER=round(MAX_DIST_TO_CONSIDER/0.13);

UserData1=expand_cells(UserData,MAX_MARGIN,0);

% reassign dots in UserData1
n_cell=length(UserData1.cell);
for i=1:n_cell % for each cell
    %if isfield(UserData.cell(1),'rim_edge'),
    %   edges=UserData.cell(i).rim_edge;
    if isfield(UserData1.cell(i),'extended_edge_new') ...
            & ~isnan(UserData1.cell(i).extended_edge_new),
        edges=UserData1.cell(i).extended_edge_new;
    else
        edges=UserData1.cell(i).edge;
    end
    %in = inpoly(UserData.dots(:,1:2),UserData.cell(i).edge);
    in = inpoly(UserData1.dots(:,1:2),edges);
    UserData1.dots(in,5) = i;
    % if a dot is in multiple cells, it will be assigned to the one
    % with larger label
end

dist_dot_cell1=1000*ones(size(UserData.dots,1),1);
dist_dot_cell2=1000*ones(size(UserData.dots,1),1);
index=find(~isnan(UserData1.dots(:,5)));

for kk=1:length(index)
    i=index(kk);
%     if mod(i,100)==0
%         i
%     end
    z=round(UserData1.dots(i,3));
    j=UserData1.dots(i,5);
    
    contour1=UserData1.nuc_contour_stack{j}{1}{z};
    if ~isempty(contour1)
        dist_dot_cell1(i) = p_poly_dist(UserData1.dots(i,1), UserData1.dots(i,2), contour1(:,1), contour1(:,2));
    end
    if length(UserData.nuc_contour_stack{j})>1,
        contour2=UserData.nuc_contour_stack{j}{2}{z};
        if ~isempty(contour2)
            dist_dot_cell2(i) = p_poly_dist(UserData1.dots(i,1), UserData1.dots(i,2), contour2(:,1), contour2(:,2));
        end
    end    
end

dist_dot_cell1(dist_dot_cell1<0)=1000;
dist_dot_cell2(dist_dot_cell2<0)=1000;

dist_dot_cell1(dist_dot_cell1>MAX_DIST_TO_CONSIDER)=1000;
dist_dot_cell2(dist_dot_cell2>MAX_DIST_TO_CONSIDER)=1000;

D=min(dist_dot_cell1,dist_dot_cell2);
distance=D;

% stratify by channels
num_channels=max(UserData.dots(:,4))
for i=1:num_channels
    Dt{i}=D(find(UserData1.dots(:,4)==i));
    ind=find(Dt{i}<1000);
    Dt{i}=Dt{i}(ind);
end
D=D(D<1000);
figure;
clrs={'r','g','b'};
for i=1:num_channels
    %subplot(num_channels,1,i);
    subplot(2,2,i);
    [p,x]=hist(Dt{i},linspace(0,max(D),15));
    x=x*0.13;
    bar(x,p);
    axis square;
    xlabel('Distance from nucleus (um)','FontSize',18);
    ylabel('Number of dots','FontSize',18);
    title(UserData.channel_name{i},'FontSize',24);
    m=median(Dt{i}*0.13);
    if i>1
        p1=ranksum(Dt{i},Dt{1});
        title([UserData.channel_name{i} ', pval=' num2str(p1)],'FontSize',24);
    end
    line([m m],ylim,'linestyle','--','color','k');
    axis tight;
    set(gca,'FontSize',18);
    subplot(2,2,4);hold on;
    vec=p/sum(p);
    plot(x,cumsum(vec),'linewidth',2,'color',clrs{i});    
end
subplot(2,2,4);
axis square;
set(gca,'FontSize',18);
box on;
legend(UserData.channel_name,'Location','northwest');
xlabel('Distance from nucleus (um)','FontSize',18);
ylabel('Fraction of dots with smaller distance','FontSize',18);
axis tight

