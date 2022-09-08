function [mod_int] = layerizemod_interp(modn,zinterp)

% Layerize model
[mod] = layerizemod(modn);

inoh2o = find(mod.vs~=0);
ih2o = find(mod.vs==0);
zh2o = mod.z(ih2o(end));

% Deal with discontinuities
mod.z(inoh2o(2:end)) = mod.z(inoh2o(2:end))+(1:length(mod.z(inoh2o(2:end))))'*1e-10;
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

if 0
    figure(9999); clf; hold on;
    h = plotlayermods(modn(:,1),modn(:,3),'o-k');
    h.LineWidth = 2;
    plot(mod.vs,mod.z,'o--b');
    plot(mod_int.vs,mod_int.z,'s--r');
end

end

