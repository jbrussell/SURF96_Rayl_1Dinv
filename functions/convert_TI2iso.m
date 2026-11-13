% Convert a transversely isotropic (TI) model to an equivalent isotropic
% model using Voigt-average velocities
%
% [model, xi, phi, eta] = convert_TI2iso(model_TI)
%
% model_TI: [H Vpv Vsv rho Vph Vsh eta] or [H Vpv Vsv rho Vph Vsh eta Qp Qs]
% model:    [H Vp Vs rho] or [H Vp Vs rho Qp Qs]
%     Vs^2 = (2*Vsv^2 + Vsh^2)/3
%     Vp^2 = (Vpv^2 + 4*Vph^2)/5
% xi  = (Vsh/Vsv)^2   (NaN in fluid layers)
% phi = (Vpv/Vph)^2
% eta = F/(A-2L)
%
% This is the inverse of convert_iso2TI.m
%
% jbrussell 11/2026
%
function [model, xi, phi, eta] = convert_TI2iso(model_TI)

vpv = model_TI(:,2);
vsv = model_TI(:,3);
rho = model_TI(:,4);
vph = model_TI(:,5);
vsh = model_TI(:,6);
eta = model_TI(:,7);

vs = sqrt((2*vsv.^2 + vsh.^2)/3);
vp = sqrt((vpv.^2 + 4*vph.^2)/5);
xi = vsh.^2 ./ vsv.^2;
xi(vsv==0) = NaN;
phi = vpv.^2 ./ vph.^2;

model = [model_TI(:,1), vp, vs, rho];
if size(model_TI,2) >= 9
    model = [model, model_TI(:,8), model_TI(:,9)];
end

end
