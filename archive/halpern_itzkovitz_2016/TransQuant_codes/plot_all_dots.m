function [ndot,ndot_nuc] = plot_all_dots(UserData)

% Plots both nuclear and cytoplasmic dots

ndot=0;
ndot_nuc=0;

if ~isempty(UserData.dots),
    ndot=plot_all_dot(UserData.dots,{'o'});
end
if ~isempty(UserData.nuclear_dots),
    ndot_nuc=plot_all_dot(UserData.nuclear_dots,{'x'});
end

