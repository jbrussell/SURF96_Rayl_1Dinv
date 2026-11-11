function [mod, discs, zlays] = card2mod_TI(cardfile,zmax)
% Convert MINEOS card to a transversely isotropic (TI) layered model for
% input to the *_TI surf96 wrappers (TI analog of card2mod.m)
%
% INPUT
% cardfile: path to card file
% zmax: maximum depth of model
%
% OUTPUT
% mod: TI layered model [H Vpv Vsv rho Vph Vsh eta] (km, km/s, g/cm^3)
%      (see writemod_surf96_TI.m). card2mod.m returns the Voigt-average
%      isotropic model with the same layering.
% discs: depth to discontinuities in card file
% zlays: depth of the bottom of each layer
%
% jbrussell 11/2026

card=read_model_card(cardfile);

% cut card at max depth and flip so Earth's surface at index=1
ind = find(card.z <= zmax);
z = flip(card.z(ind));
vpv = flip(card.vpv(ind))/1000;
vph = flip(card.vph(ind))/1000;
vsv = flip(card.vsv(ind))/1000;
vsh = flip(card.vsh(ind))/1000;
rho = flip(card.rho(ind))/1000;
eta = flip(card.eta(ind));

% Remove redundant card knots, if any (same test as card2mod.m)
vp = sqrt(1/5*vpv.^2+4/5*vph.^2);
Idup = find(diff(vp)==0 & diff(z)==0);
z(Idup) = []; vpv(Idup) = []; vph(Idup) = []; vsv(Idup) = []; vsh(Idup) = [];
rho(Idup) = []; eta(Idup) = [];

% MINEOS cards define velocity at a node, but surf96 takes layers of
% constant velocity. Assign layer values as the average between card nodes.
avg = @(x) [0.5*(x(1:end-1)+x(2:end)); x(end)];
vpv = avg(vpv); vph = avg(vph); vsv = avg(vsv); vsh = avg(vsh);
rho = avg(rho); eta = avg(eta);

% Save discontinuities
discs = unique(z(diff(z)==0));

% Layerize model by removing discontinuities
ind_nodis = find(diff(z)~=0);
z = z(ind_nodis);
vpv = vpv(ind_nodis); vph = vph(ind_nodis); vsv = vsv(ind_nodis);
vsh = vsh(ind_nodis); rho = rho(ind_nodis); eta = eta(ind_nodis);

% Build model
mod(:,1) = diff(z(:)); % layer thicknesses
mod(end,1) = 0;
mod(:,2) = vpv(1:end-1);
mod(:,3) = vsv(1:end-1);
mod(:,4) = rho(1:end-1);
mod(:,5) = vph(1:end-1);
mod(:,6) = vsh(1:end-1);
mod(:,7) = eta(1:end-1);

zlays = z(2:end);

end
