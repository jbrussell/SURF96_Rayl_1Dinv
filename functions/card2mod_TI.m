function [mod, discs, zlays] = card2mod_TI(cardfile,zmax,varargin)
% Convert MINEOS card to a transversely isotropic (TI) layered model for
% input to the *_TI surf96 wrappers (TI analog of card2mod.m)
%
% INPUT
% cardfile: path to card file
% zmax: maximum depth of model
% withQ: (OPTIONAL, default 0) if 1, append Qp and Qs columns from the
%        card's Qkappa and Qmu (see below)
%
% OUTPUT
% mod: TI layered model [H Vpv Vsv rho Vph Vsh eta] (km, km/s, g/cm^3)
%      (see writemod_surf96_TI.m). card2mod.m returns the Voigt-average
%      isotropic model with the same layering.
%      With withQ = 1: [H Vpv Vsv rho Vph Vsh eta Qp Qs], where
%        Qs = Qmu  and
%        1/Qp = L/Qmu + (1-L)/Qkappa,  L = 4/3 (Vs/Vp)^2
%      (Voigt-average Vs and Vp; Qs = 0 in fluid layers).
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
qkap = flip(card.qkap(ind));
qmu = flip(card.qmu(ind));

% Remove redundant card knots, if any (same test as card2mod.m)
vp = sqrt(1/5*vpv.^2+4/5*vph.^2);
Idup = find(diff(vp)==0 & diff(z)==0);
z(Idup) = []; vpv(Idup) = []; vph(Idup) = []; vsv(Idup) = []; vsh(Idup) = [];
rho(Idup) = []; eta(Idup) = []; qkap(Idup) = []; qmu(Idup) = [];

% MINEOS cards define velocity at a node, but surf96 takes layers of
% constant velocity. Assign layer values as the average between card nodes.
avg = @(x) [0.5*(x(1:end-1)+x(2:end)); x(end)];
vpv = avg(vpv); vph = avg(vph); vsv = avg(vsv); vsh = avg(vsh);
rho = avg(rho); eta = avg(eta); qkap = avg(qkap); qmu = avg(qmu);

% Save discontinuities
discs = unique(z(diff(z)==0));

% Layerize model by removing discontinuities
ind_nodis = find(diff(z)~=0);
z = z(ind_nodis);
vpv = vpv(ind_nodis); vph = vph(ind_nodis); vsv = vsv(ind_nodis);
vsh = vsh(ind_nodis); rho = rho(ind_nodis); eta = eta(ind_nodis);
qkap = qkap(ind_nodis); qmu = qmu(ind_nodis);

% Build model
mod(:,1) = diff(z(:)); % layer thicknesses
mod(end,1) = 0;
mod(:,2) = vpv(1:end-1);
mod(:,3) = vsv(1:end-1);
mod(:,4) = rho(1:end-1);
mod(:,5) = vph(1:end-1);
mod(:,6) = vsh(1:end-1);
mod(:,7) = eta(1:end-1);

if nargin >= 3 && ~isempty(varargin{1}) && varargin{1}
    vs = sqrt((2*mod(:,3).^2 + mod(:,6).^2)/3);
    vp = sqrt((mod(:,2).^2 + 4*mod(:,5).^2)/5);
    [mod(:,8), mod(:,9)] = q_kapmu2pqs(qkap(1:end-1), qmu(1:end-1), vp, vs);
end

zlays = z(2:end);

end

function [qp, qs] = q_kapmu2pqs(qkap, qmu, vp, vs)
% Qp and Qs from Qkappa and Qmu (Q values of 0 are treated as elastic)
iqk = zeros(size(qkap)); iqk(qkap > 0) = 1./qkap(qkap > 0);
iqm = zeros(size(qmu));  iqm(qmu > 0)  = 1./qmu(qmu > 0);
L = 4/3 * (vs./vp).^2;
iqp = L.*iqm + (1-L).*iqk;
qp = zeros(size(iqp)); qp(iqp > 0) = 1./iqp(iqp > 0);
qs = zeros(size(iqm)); qs(iqm > 0) = 1./iqm(iqm > 0);
qs(vs <= 0.01) = 0; % fluid
end
