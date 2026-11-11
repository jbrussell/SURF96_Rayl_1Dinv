% Convert an isotropic model to a transversely isotropic (TI) model
%
% model_TI = convert_iso2TI(model)
% model_TI = convert_iso2TI(model, xi, phi, eta)
%
% model:  isotropic model [H Vp Vs rho] or [H Vp Vs rho Qp Qs]
% xi:     (Vsh/Vsv)^2   (OPTIONAL, scalar or nlayer vector, default 1)
% phi:    (Vpv/Vph)^2   (OPTIONAL, scalar or nlayer vector, default 1)
% eta:    F/(A-2L)      (OPTIONAL, scalar or nlayer vector, default 1)
%
% model_TI: [H Vpv Vsv rho Vph Vsh eta] or [H Vpv Vsv rho Vph Vsh eta Qp Qs]
%           (the format used by the *_TI functions, see writemod_surf96_TI.m)
%
% With the default xi = phi = eta = 1 the TI model is exactly equivalent to
% the isotropic model (Vpv=Vph=Vp, Vsv=Vsh=Vs). Otherwise, the input Vp and
% Vs are treated as Voigt-average velocities, which are preserved:
%     Vs^2 = (2*Vsv^2 + Vsh^2)/3   ->  Vsv^2 = 3*Vs^2/(2+xi),  Vsh^2 = xi*Vsv^2
%     Vp^2 = (Vpv^2 + 4*Vph^2)/5   ->  Vph^2 = 5*Vp^2/(4+phi), Vpv^2 = phi*Vph^2
%
% jbrussell 11/2026
%
function model_TI = convert_iso2TI(model,varargin)

nlay = size(model,1);
xi = ones(nlay,1);
phi = ones(nlay,1);
eta = ones(nlay,1);
if nargin >= 2 && ~isempty(varargin{1}); xi  = varargin{1}(:).*ones(nlay,1); end
if nargin >= 3 && ~isempty(varargin{2}); phi = varargin{2}(:).*ones(nlay,1); end
if nargin >= 4 && ~isempty(varargin{3}); eta = varargin{3}(:).*ones(nlay,1); end

vp = model(:,2);
vs = model(:,3);
rho = model(:,4);

vsv = sqrt(3*vs.^2 ./ (2+xi));
vsh = sqrt(xi) .* vsv;
vph = sqrt(5*vp.^2 ./ (4+phi));
vpv = sqrt(phi) .* vph;

% keep exact isotropic values where no anisotropy is requested
iso_s = (xi==1);
vsv(iso_s) = vs(iso_s); vsh(iso_s) = vs(iso_s);
iso_p = (phi==1);
vpv(iso_p) = vp(iso_p); vph(iso_p) = vp(iso_p);

model_TI = [model(:,1), vpv, vsv, rho, vph, vsh, eta];
if size(model,2) >= 6
    model_TI = [model_TI, model(:,5), model(:,6)];
end

end
