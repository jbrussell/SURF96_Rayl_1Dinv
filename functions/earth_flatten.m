function [ mod_fl ] = earth_flatten( mod )
% Apply earth flattening transformation from Muller 1971
%
% jbrussell 6/4/2020

R = 6371; % radius of earth
z = cumsum(mod(:,1));

zfac = R./(R - z);

% Calculate z flat
zflat = R * log(zfac);

% Calculate v flat
vp = mod(:,2);
vp_flat = zfac .* vp;

vs = mod(:,3);
vs_flat = zfac .* vs;

rho = mod(:,4);
% p = 4;
p = -2.275;
rho_flat = (zfac).^p .* rho;


hflat = diff([0; zflat]);
mod_fl = [hflat, vp_flat, vs_flat, rho_flat];

if 1
    figure(87); clf; hold on;
    plotlayermods(mod(:,1),mod(:,3),'-k');
    plotlayermods(mod_fl(:,1),mod_fl(:,3),'-r');
    legend({'input','flattened'},'location','southwest');
end

end

