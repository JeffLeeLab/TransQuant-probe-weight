function ncell = plot_all_cell(cell,scale,edge_color,label_color,labelfontsize)

ncell = length(cell);

if ncell
    if nargin <2
        scale = [1 1];
    end
    if nargin < 3
        edge_color = 'w';
    end
    if nargin <4
        label_color = 'g';
    end
    if nargin <5
        labelfontsize = 14;
    end
    hold on

    for i=1:ncell
        if isfield(cell(i),'extended_edge_new') & ~isnan(cell(i).extended_edge_new)
            x = cell(i).extended_edge_new(:,1)*scale(1);
            y = cell(i).extended_edge_new(:,2)*scale(2);
        else
            x = cell(i).edge(:,1)*scale(1);
            y = cell(i).edge(:,2)*scale(2);
        end
        xc = cell(i).center(1)*scale(1);
        yc = cell(i).center(2)*scale(2);
        plot(x,y,'Color',edge_color(1),'Tag',['plot_cell_edge' num2str(i)],'LineWidth',2);
        
        if isfield(cell(i),'rim_edge') & ~isnan(cell(i).rim_edge),
            x2 = cell(i).rim_edge(:,1)*scale(1);
            y2 = cell(i).rim_edge(:,2)*scale(2);            
            plot(x2,y2,'Color',edge_color(1),'Tag',['plot_cell_edge' num2str(i)],'LineWidth',2,'linestyle',':');
        end    
        
        text(xc,yc,cell(i).label,'FontSize',labelfontsize,'Color',label_color(1),'HorizontalAlignment','center','Tag',['plot_cell_label' cell(i).label])
    end
else
    disp('No cells!')
end