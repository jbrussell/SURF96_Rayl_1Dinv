function [mod_out] = spline2mod(modsp,vs_sp,vp_vs,rho_vs)

dz = [diff(modsp.z(:)); 0];

vp = vp_vs * vs_sp; vp(vs_sp<0.01)=1.5;
% rho = vp / 2.5; rho(1)=1.03;
rho = rho_vs * vs_sp; rho(vs_sp<0.01)=1.03;

idisc = find(dz(1:end-1) == 0);
dz(idisc) = [];
vp(idisc) = [];
vs_sp(idisc) = [];
rho(idisc) = [];

I = find(vs_sp>0);
dzshift = dz(I(1))/2;
dz(I(1)) = dz(I(1)) - dzshift;
dz(end-1) = dz(end-1) + dzshift;

% refmod_sp = [dz mod_ref.vp vs_ref mod_ref.rho];
mod_out = [dz vp vs_sp rho];

end

