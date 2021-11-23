function [ cardo ] = mod2card( mod,cardn )
% Convert surf96 layered model to MINEOS card file
%

[MODEL] = layerizemod(mod);

zmax = max(MODEL.z);

ind_replace = cardn.z <= zmax;
cardo = cardn;
cardo.z = [cardo.z(~ind_replace); flip(MODEL.z)];
cardo.rad = [cardo.rad(~ind_replace); flip((6371-MODEL.z)*1000)];
cardo.vpv = [cardo.vpv(~ind_replace); flip(MODEL.vp*1000)];
cardo.vph = [cardo.vph(~ind_replace); flip(MODEL.vp*1000)];
cardo.vsv = [cardo.vsv(~ind_replace); flip(MODEL.vs*1000)];
cardo.vsh = [cardo.vsh(~ind_replace); flip(MODEL.vs*1000)];
cardo.rho = [cardo.rho(~ind_replace); flip(MODEL.rho*1000)];
cardo.eta = [cardo.eta(~ind_replace); flip(ones(length(MODEL.z),1))];
cardo.qmu = [cardo.qmu(~ind_replace); flip(ones(length(MODEL.z),1) * mean(cardo.qmu(ind_replace)))];
cardo.qkap = [cardo.qkap(~ind_replace); flip(ones(length(MODEL.z),1) * mean(cardo.qkap(ind_replace)))];




end

