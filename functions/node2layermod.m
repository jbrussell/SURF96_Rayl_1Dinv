function [layermod, discs, zlays] = node2layermod(z_node,vp_node,vs_node,rho_node)
% Convert model with velocity at depth nodes (i.e. MINEOS) to
% model with constant velocity layers.
%
% INPUT
% z_node: depth of node (earth's surface at z_node(1))
% vp_node: vp of node
% vs_node: vs of node
% rho_node: density of node
% zmax: maximum depth of desired model
%
% OUTPUT
% layermod: layered model file for input to surf96
% discs: depth to discontinuities in card file
%
% jbrussell 6/15/2020


z = z_node;
vp = vp_node;
vs = vs_node;
rho = rho_node;

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
layermod(:,1) = diff(z(:)); % layer thicknesses
layermod(end,1) = 0;
layermod(:,2) = vp(1:end-1); % vp
layermod(:,3) = vs(1:end-1); % vs
layermod(:,4) = rho(1:end-1); % density

zlays = z(2:end);
    
if 0
    figure(99); clf; hold on;
    plot(vs_node/1000,z_node*-1,'-k','linewidth',2)
    h = plotlayermods(layermod(:,1),layermod(:,3));
    set(h,'linewidth',2,'color',[1 0 0]);
    legend('card','mod','location','southwest');
    xlabel('Vs');
    ylabel('Depth (km)');
    set(gca,'fontsize',15,'linewidth',1.5,'box','on'); 
end

end