function [mod, discs, zlays] = card2mod(cardfile,zmax)
% Convert MINEOS card to layered model for input to surf96
%
% INPUT
% cardfile: path to card file
% zmax: maximum depth of model
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

% Remove redundant card knots, if any
Idup = find(diff(vp)==0 & diff(z)==0);
z(Idup) = [];
vp(Idup) = [];
vs(Idup) = [];
rho(Idup) = [];

% MINEOS cards define velocity at a node, but surf96 takes layers of
% constant velocity. We can can approximate this by assigning surf96 layer
% velocities as average between card nodes.
vp(1:end-1) = 0.5*(vp(1:end-1)+vp(2:end));
vs(1:end-1) = 0.5*(vs(1:end-1)+vs(2:end));
rho(1:end-1) = 0.5*(rho(1:end-1)+rho(2:end));

% Save discontinuities
discs = unique(z(diff(z)==0));

% Layerize model by removing discontinuities
ind_nodis = find(diff(z)~=0);
z = z(ind_nodis);
vp = vp(ind_nodis);
vs = vs(ind_nodis);
rho = rho(ind_nodis);

% Build model
mod(:,1) = diff(z(:)); % layer thicknesses
mod(end,1) = 0;
mod(:,2) = vp(1:end-1); % vp
mod(:,3) = vs(1:end-1); % vs
mod(:,4) = rho(1:end-1); % density
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