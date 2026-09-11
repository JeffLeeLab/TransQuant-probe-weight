function I=sum_max_stack(pixel_list,pixel_values)

z=unique(pixel_list(:,3));
for i=1:length(z),
    index=find(pixel_list(:,3)==z(i));
    I(i)=sum(pixel_values(index));
end
I=max(I);

% Add 2d gaussian estimation
