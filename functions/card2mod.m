function [mod, discs, zlays] = card2mod(cardfile,zmax,varargin)
% Convert MINEOS card to layered model for input to surf96
%
% INPUT
% cardfile: path to card file
% zmax: maximum depth of model
% withQ: (OPTIONAL, default 0) if 1, append Qp and Qs columns:
%        mod = [H Vp Vs rho Qp Qs] with Qs = Qmu and
%        1/Qp = L/Qmu + (1-L)/Qkappa, L = 4/3 (Vs/Vp)^2 (Qs = 0 in fluid)
%        (jbrussell 11/2026)
%
% OUTPUT
% mod: layered model file for input to surf96
% discs: depth to discontinuities in card file
%
% jbrussell 5/28/2020

card=read_model_card(cardfile);

% cut card at max depth and flip so Earth's surface at index=1
ind = find(card.z <= zmax);
z = flip(card.z(ind));
% vp = flip(0.5*(card.vpv(ind)+card.vph(ind)))/1000;
% vs = flip(0.5*(card.vsv(ind)+card.vsh(ind)))/1000;
vp = flip(sqrt(1/5*card.vpv(ind).^2+4/5*card.vph(ind).^2))/1000;
vs = flip(sqrt(2/3*card.vsv(ind).^2+1/3*card.vsh(ind).^2))/1000;
rho = flip(card.rho(ind))/1000;
qkap = flip(card.qkap(ind));
qmu = flip(card.qmu(ind));

% Remove redundant card knots, if any
Idup = find(diff(vp)==0 & diff(z)==0);
z(Idup) = [];
vp(Idup) = [];
vs(Idup) = [];
rho(Idup) = [];
qkap(Idup) = [];
qmu(Idup) = [];

% MINEOS cards define velocity at a node, but surf96 takes layers of
% constant velocity. We can can approximate this by assigning surf96 layer
% velocities as average between card nodes.
vp(1:end-1) = 0.5*(vp(1:end-1)+vp(2:end));
vs(1:end-1) = 0.5*(vs(1:end-1)+vs(2:end));
rho(1:end-1) = 0.5*(rho(1:end-1)+rho(2:end));
qkap(1:end-1) = 0.5*(qkap(1:end-1)+qkap(2:end));
qmu(1:end-1) = 0.5*(qmu(1:end-1)+qmu(2:end));

% Save discontinuities
discs = unique(z(diff(z)==0));

% Layerize model by removing discontinuities
ind_nodis = find(diff(z)~=0);
z = z(ind_nodis);
vp = vp(ind_nodis);
vs = vs(ind_nodis);
rho = rho(ind_nodis);
qkap = qkap(ind_nodis);
qmu = qmu(ind_nodis);

% Build model
mod(:,1) = diff(z(:)); % layer thicknesses
mod(end,1) = 0;
mod(:,2) = vp(1:end-1); % vp
mod(:,3) = vs(1:end-1); % vs
mod(:,4) = rho(1:end-1); % density
if nargin >= 3 && ~isempty(varargin{1}) && varargin{1}
    iqk = zeros(size(qkap)); iqk(qkap > 0) = 1./qkap(qkap > 0);
    iqm = zeros(size(qmu));  iqm(qmu > 0)  = 1./qmu(qmu > 0);
    iqk = iqk(1:end-1); iqm = iqm(1:end-1);
    L = 4/3 * (mod(:,3)./mod(:,2)).^2;
    iqp = L.*iqm + (1-L).*iqk;
    mod(:,5) = 0; mod(iqp > 0,5) = 1./iqp(iqp > 0); % Qp
    mod(:,6) = 0; mod(iqm > 0,6) = 1./iqm(iqm > 0); % Qs
    mod(mod(:,3) <= 0.01,6) = 0; % fluid
end
% mod(:,2) = 0.5*(vp(1:end-1) + vp(2:end)); % vp
% mod(:,3) = 0.5*(vs(1:end-1) + vs(2:end)); % vs
% mod(:,4) = 0.5*(rho(1:end-1) + rho(2:end)); % density

zlays = z(2:end);
    
if 0
    figure(99); clf; hold on;
    plot(card.vsv(ind)/1000,card.z(ind)*-1,'-k','linewidth',2)
    h = plotlayermods(mod(:,1),mod(:,3));
    set(h,'linewidth',2,'color',[1 0 0]);
    legend('card','mod','location','southwest');
    xlabel('Vs');
    ylabel('Depth (km)');
    set(gca,'fontsize',15,'linewidth',1.5,'box','on'); 
end

end