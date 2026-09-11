function ndot = plot_all_dot(dots,marker,scale)

ndot = length(dots)
if ndot
    hold on;
    
    if nargin <3
        scale = [1 1];
    end
    if nargin <2
        marker = 'o';
    end
    
    if marker{1}=='x',
        col='mcgw';
        marker_size=[9 9 9];
    else
        col='rgbm';
        marker_size=[5 5 5];
    end
    
    %marker={'s','^','o'};
    %marker_size=[5 6 3];
    scale
    dots(:,1) = dots(:,1)*scale(1);
    dots(:,2) = dots(:,2)*scale(2);
    n_channel=numel(unique(dots(:,4)));
    n_channel
    for i=1:n_channel
        %i
        p = dots(dots(:,4)==i,:);
        %p
        %plot(p(:,1),p(:,2),marker,'Color',col(i),'Tag','plot_dot','MarkerSize',3,'MarkerFaceColor',col(i));
        plot(p(:,1),p(:,2),marker{1},'Color',col(i),'Tag','plot_dot','MarkerSize',marker_size(i),'MarkerFaceColor',col(i));
    end
else
    disp('No dots!')
end