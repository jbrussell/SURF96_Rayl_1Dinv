% Plot 1-D layered isotropic or transversely isotropic (TI) models
%
% plot_model_TI(models)
% plot_model_TI(models, labels)
% plot_model_TI(models, labels, zmax)
%
% models: a model matrix or a cell array of model matrices, each either
%         isotropic [H Vp Vs rho (Qp Qs)] or TI [H Vpv Vsv rho Vph Vsh eta (Qp Qs)]
% labels: cell array of legend labels (OPTIONAL)
% zmax:   maximum depth to plot (km) (OPTIONAL, default: 20 km below the
%         top of the half-space)
%
% Panels (in the current figure): Vp (Vpv solid, Vph dashed), Vs (Vsv
% solid, Vsh dashed), density, anisotropy (xi = (Vsh/Vsv)^2, phi =
% (Vpv/Vph)^2, eta) and, if any model has Q, Qp and Qs. Layers are drawn as
% constant-property steps; the half-space is drawn down to zmax.
%
% jbrussell 11/2026
%
function plot_model_TI(models, varargin)

if ~iscell(models)
    models = {models};
end
nmod = length(models);
labels = {};
if nargin >= 2 && ~isempty(varargin{1})
    labels = varargin{1};
end
% convert to TI format [H Vpv Vsv rho Vph Vsh eta Qp Qs] (Q NaN if absent)
hasQ = false;
for im = 1:nmod
    m = models{im};
    switch size(m,2)
        case 4
            m = [convert_iso2TI(m), nan(size(m,1),2)];
        case 6
            m = convert_iso2TI(m);
        case 7
            m = [m, nan(size(m,1),2)];
        case 9
        otherwise
            error('plot_model_TI: models must have 4, 6, 7 or 9 columns');
    end
    hasQ = hasQ || any(~isnan(m(:,8)));
    models{im} = m;
end
zhs = max(cellfun(@(m) sum(m(1:end-1,1)), models));
zmax = zhs + 20;
if nargin >= 3 && ~isempty(varargin{2})
    zmax = varargin{2};
end

npan = 4 + hasQ;
clr = lines(nmod);
ax = zeros(1,npan);
for ip = 1:npan
    ax(ip) = subplot(1,npan,ip); hold on; box on;
end
hleg = zeros(1,nmod);
for im = 1:nmod
    m = models{im};
    [z, v] = steps(m, zmax);
    c = clr(im,:);
    vsv = v(:,3); vsh = v(:,6);
    xi = vsh.^2 ./ vsv.^2; xi(vsv == 0) = NaN;
    phi = v(:,2).^2 ./ v(:,5).^2;
    % Vp
    axes(ax(1));
    hleg(im) = plot(v(:,2), z, '-', 'color', c, 'linewidth', 2);
    plot(v(:,5), z, '--', 'color', c, 'linewidth', 1.5);
    % Vs
    axes(ax(2));
    plot(vsv, z, '-', 'color', c, 'linewidth', 2);
    plot(vsh, z, '--', 'color', c, 'linewidth', 1.5);
    % density
    axes(ax(3));
    plot(v(:,4), z, '-', 'color', c, 'linewidth', 2);
    % anisotropy
    axes(ax(4));
    plot(xi, z, '-', 'color', c, 'linewidth', 2);
    plot(phi, z, '--', 'color', c, 'linewidth', 1.5);
    plot(v(:,7), z, ':', 'color', c, 'linewidth', 2);
    % Q
    if hasQ
        axes(ax(5));
        plot(v(:,9), z, '-', 'color', c, 'linewidth', 2);
        plot(v(:,8), z, '--', 'color', c, 'linewidth', 1.5);
    end
end
xl = {'Vp (km/s)  [- Vpv, -- Vph]', 'Vs (km/s)  [- Vsv, -- Vsh]', '\rho (g/cm^3)', ...
    'anisotropy  [- \xi, -- \phi, : \eta]', 'Q  [- Qs, -- Qp]'};
for ip = 1:npan
    axes(ax(ip));
    set(gca, 'ydir', 'reverse');
    ylim([0 zmax]);
    xlabel(xl{ip});
    if ip == 1; ylabel('Depth (km)'); end
end
if hasQ
    axes(ax(5)); set(gca,'xscale','log');
end
if ~isempty(labels)
    axes(ax(1));
    legend(hleg, labels, 'location', 'southwest');
end

end

function [z, v] = steps(m, zmax)
% depth and property arrays for a step (layered) plot
ztop = [0; cumsum(m(1:end-1,1))];
zbot = [ztop(2:end); max(zmax, ztop(end)+1)];
n = size(m,1);
z = zeros(2*n,1);
v = zeros(2*n,size(m,2));
z(1:2:end) = ztop;
z(2:2:end) = zbot;
v(1:2:end,:) = m;
v(2:2:end,:) = m;
end
