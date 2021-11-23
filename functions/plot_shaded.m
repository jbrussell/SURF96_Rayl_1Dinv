function [ h ] = plot_shaded( ax, x_l, x_u, y, opt, clr, alph )
% Plot a shaded region defined by two vectors of upper and lower bounds and 
% a third vector along the second axis.
%
% x_u : vector containing upper bound (size:1xN)
% x_l : vector containing lower bound (size:1xN)
% y   : axis along which to plot (size:1xN)
% opt : 'x' = upper and lower bounds along x-direction
%       'y' = upper and lower bounds along y-direction
% clr : [R,G,B] color
% alph: transparency from 0-1
%

if size(x_u,1)~=1
    x_u = x_u';
end
if size(x_l,1)~=1
    x_l = x_l';
end
if size(y,1)~=1
    y = y';
end

if strcmp(opt,'x')==1
    h = patch(ax,[x_l fliplr(x_u)],[y fliplr(y)],clr,'linestyle','none'); hold on;
    set(h, 'FaceAlpha', alph);
elseif strcmp(opt,'y')==1
    h = patch(ax,[y fliplr(y)],[x_l fliplr(x_u)],clr,'linestyle','none'); hold on;
    set(h, 'FaceAlpha', alph);
end


end

