function [mod_int] = layermod2node_interp(modn,zinterp)

% % Layerize model
% [mod] = layerizemod(modn);

% Build structure defined by velocities at each node at center of layer
z = [0; cumsum(modn(1:end-1,1))];
vp = [0; modn(1:end-1,2)];
vs = [0; modn(1:end-1,3)];
rho = [0; modn(1:end-1,4)];
ih2o = find(vp==1.5);
zh2o = z(ih2o(end));
vp(1) = modn(1,2);
vs(1) = modn(1,3);
rho(1) = modn(1,4);
z = z - [0; modn(1:end-1,1)]/2;
% shift seafloor back to where it should be
ilays = find(z < zh2o);
z(ilays(end)) = zh2o;

% add another node at seafloor to make discontinuity
z = [z(ilays); z(ilays(end)); z(ilays(end)+1:end)];
vp = [vp(ilays); vp(ilays(end)+1); vp(ilays(end)+1:end)];
vs = [vs(ilays); vs(ilays(end)+1); vs(ilays(end)+1:end)];
rho = [rho(ilays); rho(ilays(end)+1); rho(ilays(end)+1:end)];

mod.z = z;
mod.vp = vp;
mod.vs = vs;
mod.rho = rho;

% z = [0; cumsum(modn(1:end-1,1))];
% vp = [0; 0.5*(modn(1:end-1,2) + modn(2:end,2))];
% vs = [0; 0.5*(modn(1:end-1,3) + modn(2:end,3))];
% rho = [0; 0.5*(modn(1:end-1,4) + modn(2:end,4))];
% vp(1) = modn(1,2);
% vs(1) = modn(1,3);
% rho(1) = modn(1,4);

inoh2o = find(mod.vp~=1.5);
ih2o = find(mod.vp==1.5);

% Deal with discontinuities
mod.z(inoh2o(2:end)) = mod.z(inoh2o(2:end))+(1:length(mod.z(inoh2o(2:end))))'*1e-13;
% mod.z(inoh2o(1)) = 0;

Ih2o = find(zinterp==zh2o);
zinterp_noh2o = zinterp(Ih2o(end):end);

mod_int.vs = interp1(mod.z(inoh2o),mod.vs(inoh2o),zinterp_noh2o);
mod_int.vp = interp1(mod.z(inoh2o),mod.vp(inoh2o),zinterp_noh2o);
mod_int.rho = interp1(mod.z(inoh2o),mod.rho(inoh2o),zinterp_noh2o);
mod_int.z = zinterp_noh2o;

% Add h2o back
mod_int.vs = [mod.vs(ih2o); mod_int.vs]; 
mod_int.vp = [mod.vp(ih2o); mod_int.vp]; 
mod_int.rho = [mod.rho(ih2o); mod_int.rho];
mod_int.z = [mod.z(ih2o); mod_int.z];
% mod_int.z = mod_int.z - [0; diff(mod_int.z)/2];

if 0
    figure(9999); clf; hold on;
    h = plotlayermods(modn(:,1),modn(:,3),'o-k');
    h.LineWidth = 2;
    plot(mod.vs,mod.z,'o--b');
    plot(mod_int.vs,mod_int.z,'s--r');
    plot(vs,z,'^--c');
end

end
