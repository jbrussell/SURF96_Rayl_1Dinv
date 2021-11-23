function [h,mod,dep] = plotlayermods2(THICK,MODELS,lintype,offs)
%PLOTLAYERMODS Plot layered models
%
%H = PLOTLAYERMODS(THICK,MODELS) plots the layered model arrays in
% MODELS of layer thickness given in THICK.  At least one
% dimension of the model array must be
% equal to length(THICK).  If both are, then PLOTLAYERMODS
% assumes that each column is a different model.
% H is a handle to the plot
%PLOTLAYERMODS(THICK,MODELS,'linetype') plots the model(s) with
% the given linetype or symbol.
%PLOTLAYERMODS(THICK,MODELS,'linetype',offs) plots the model(s) with
% the given linetype or symbol, and a y-offset.

if (nargin<2)
    error('plotlayermods() needs at least two arguments')
end
[M,N] = size(MODELS);
L = length(THICK);
if (M==L) % the default way, we're happy
    L = M; %i.e., do nothing
elseif (N==L)
    MODELS = MODELS';
    [M,N] = size(MODELS);
else
    error('No dimension of model array equal to dim of depth array');
end

if (nargin ~= 4), offs = 0; end;
for i=1:L			%for all layers
    if i > 1
        dep(2*i-1,1) = dep(2*i-2,1);
    else
        dep(2*i-1,1) = offs;
    end
    dep(2*i,1) = dep(2*i-1,1)+THICK(i);
    for j = 1:N
        mod(2*i-1,j) = MODELS(i,j);
        mod(2*i,j) = MODELS(i,j);
    end
end
if (nargin==2)
    h = plot(mod,dep);
else
    h = plot(mod,dep,lintype);
end
%axis('ij');

