function [kernel, dispersion] = calc_anelastic_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot, varargin)
% Calculate anelastic (Q-corrected) surface wave sensitivity kernels for a
% transversely isotropic (TI) model (TI analog of calc_anelastic_kernel96.m)
%
% [kernel, dispersion] = calc_anelastic_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot);
% [kernel, dispersion] = calc_anelastic_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot, nmode, fref);
%
% model: TI model [H Vpv Vsv rho Vph Vsh eta Qp Qs]  (see writemod_surf96_TI.m)
% vec_T: period in seconds
% wavetype: 'L' or 'R'
% ifnorm: 1: normalize kernel by layer thickness; 0: no normalization
% ifplot: whether to plot the kernel
% nmode: 0 for fund. mode; 1 for 1st overtone, etc. (OPTIONAL, default 0)
% fref: reference frequency (Hz) of the constant-Q model (OPTIONAL, default 1)
%
% METHOD
% CPS has no TI version of srfker96, so the anelastic correction applied by
% srfker96 (subroutines getgam and cslmod in VOLIV/src/srfker96.f) is
% reproduced here, applied to the elastic TI kernels from calc_kernel96_TI.m:
%   gamma = w/(2 c0^2) * sum[ (dc/dVsv Vsv + dc/dVsh Vsh) Qs^-1
%                            + (dc/dVpv Vpv + dc/dVph Vph) Qp^-1 ]
%   c     = c0 + 2 gamma c0^2 ln(f/fref) / (pi w)
%   U     = U0 (1 + (2 - U0/c0)(c-c0)/c0 + 2 gamma U0/(pi w))
% and the velocity kernels are converted to their causal values. The S
% velocity kernels (Vsv, Vsh) use Qs and the P velocity kernels (Vpv, Vph)
% use Qp. As in srfker96, the thickness kernels are not modified; the eta
% and density kernels are also left as their elastic values.
%
% OUTPUT
% kernel: structure of [nlayer x nperiod] kernels with the same fields as
%   calc_anelastic_kernel96.m (dcdvs, dcdvp, dudvs, dudvp, dcdQs, dcdQp,
%   dudQs, dudQp, dgdQs, dgdQp, dQdQs, dQdQp, dgdQmu, dgdQkap, dQdQmu,
%   dQdQkap, z) plus the TI kernels of calc_kernel96_TI.m (dcdvpv, dcdvph,
%   dcdvsv, dcdvsh, dcdeta, dcdrho, dcdh and the dU equivalents).
%   NOTE: here dcdrho/dudrho are true density kernels and dcdh/dudh are
%   thickness kernels (in calc_anelastic_kernel96.m the fields named
%   dcdrho/dudrho hold srfker96's dc/dh and dU/dh).
%   Qs and Qp kernels are with respect to Qs^-1 and Qp^-1. The Qmu/Qkappa
%   conversion uses Voigt-average Vs and Vp.
% dispersion: structure with fields period, mode, phv, grv, gamma, qinv,
%   phv_elastic, grv_elastic
%
% jbrussell 11/2026
%

if size(model,2)<9
    error('Must include Qp and Qs in TI model [H Vpv Vsv rho Vph Vsh eta Qp Qs] to generate anelastic kernels');
end

if nargin == 5
    nmode = 0; % default fundamental mode
    fref = 1; % default reference frequency
elseif nargin == 6
    nmode = varargin{1};
    fref = 1; % default reference frequency
elseif nargin == 7
    nmode = varargin{1};
    fref = varargin{2};
end
wavetype = upper(wavetype);

vec_T = vec_T(:)';
Nper = length(vec_T);
nlayer = size(model,1);

% Elastic kernels (not normalized)
[kel, disp_el] = calc_kernel96_TI(model(:,1:7), vec_T, wavetype, 0, 0, nmode);
if isempty(kel)
    disp('Something went wrong calculating kernels');
    kernel = []; dispersion = [];
    return
end

% Model velocities and Q^-1 (CPS convention: Q>1 is converted to 1/Q)
vpv = model(:,2); vsv = model(:,3); vph = model(:,5); vsh = model(:,6);
qa = model(:,8); qb = model(:,9);
qa(qa>1) = 1./qa(qa>1);
qb(qb>1) = 1./qb(qb>1);
qb(vsv<=0.01) = 0; % fluid

Sparams = {'vsv','vsh'}; Svel = [vsv, vsh];
Pparams = {'vpv','vph'}; Pvel = [vpv, vph];
other = {'eta','rho','h'};

for ii = 1:2
    for pre = {'dcd','dud'}
        kernel.([pre{1},Sparams{ii}]) = nan(nlayer,Nper);
        kernel.([pre{1},Pparams{ii}]) = nan(nlayer,Nper);
    end
end
for ii = 1:length(other)
    kernel.(['dcd',other{ii}]) = kel.(['dcd',other{ii}]);
    kernel.(['dud',other{ii}]) = kel.(['dud',other{ii}]);
end
qflds = {'dcdQs','dcdQp','dudQs','dudQp','dgdQs','dgdQp'};
for ii = 1:length(qflds)
    kernel.(qflds{ii}) = nan(nlayer,Nper);
end

dispersion.period = vec_T(:);
dispersion.mode = nmode;
dispersion.phv = nan(Nper,1);
dispersion.grv = nan(Nper,1);
dispersion.gamma = nan(Nper,1);
dispersion.phv_elastic = disp_el.phv(:);
dispersion.grv_elastic = disp_el.grv(:);

for ip = 1:Nper
    c0 = disp_el.phv(ip);
    u0 = disp_el.grv(ip);
    if isnan(c0)
        continue
    end
    f = 1/vec_T(ip);
    omega = 2*pi*f;
    faclog = log(f/fref)/pi;

    % sum of dc/dv * v for S and P velocities (elastic kernels)
    sumS = zeros(nlayer,1); sumP = zeros(nlayer,1);
    for ii = 1:2
        sumS = sumS + kel.(['dcd',Sparams{ii}])(:,ip) .* Svel(:,ii);
        sumP = sumP + kel.(['dcd',Pparams{ii}])(:,ip) .* Pvel(:,ii);
    end
    sumS(isnan(sumS)) = 0; sumP(isnan(sumP)) = 0;

    % attenuation coefficient (getgam)
    gam = 0.5*omega/(c0*c0) * sum(sumS.*qb + sumP.*qa);

    % causal phase and group velocities (cslmod)
    facg = 2*gam*u0/(pi*omega);
    c = c0 + 2*gam*c0*c0*faclog/omega;
    cmc0c0 = (c-c0)/c0;
    u0c0 = u0/c0;
    u = u0*(1 + (2-u0c0)*cmc0c0 + facg);
    uu0 = u/u0;

    dispersion.phv(ip) = c;
    dispersion.grv(ip) = u;
    dispersion.gamma(ip) = gam;

    % causal velocity kernels (cslmod)
    for ii = 1:2
        for jj = 1:2
            if jj==1
                p = Sparams{ii}; q = qb;
            else
                p = Pparams{ii}; q = qa;
            end
            dcdp = kel.(['dcd',p])(:,ip);
            dudp = kel.(['dud',p])(:,ip);
            kernel.(['dcd',p])(:,ip) = dcdp.*(1 + faclog*q);
            kernel.(['dud',p])(:,ip) = dudp.*(uu0 - u0c0*cmc0c0 + facg) ...
                + dcdp.*u0c0.*(-2*facg + u0c0*q/pi + (2-u0c0)*(faclog*q - cmc0c0) + u0c0*cmc0c0);
        end
    end

    % Q^-1 kernels
    kernel.dcdQs(:,ip) = faclog*sumS;
    kernel.dcdQp(:,ip) = faclog*sumP;
    kernel.dgdQs(:,ip) = 0.5*omega*sumS/(c0*c0);
    kernel.dgdQp(:,ip) = 0.5*omega*sumP/(c0*c0);
    kernel.dudQs(:,ip) = u0c0*(2-u0c0)*kernel.dcdQs(:,ip) + u0c0*u0c0*sumS/pi;
    kernel.dudQp(:,ip) = u0c0*(2-u0c0)*kernel.dcdQp(:,ip) + u0c0*u0c0*sumP/pi;
end
if strcmp(wavetype,'L')
    kernel.dcdQp(:) = 0; kernel.dudQp(:) = 0; kernel.dgdQp(:) = 0;
end

% Isotropic-equivalent kernels
kernel.dcdvp = kernel.dcdvpv + kernel.dcdvph;
kernel.dcdvs = kernel.dcdvsv + kernel.dcdvsh;
kernel.dudvp = kernel.dudvpv + kernel.dudvph;
kernel.dudvs = kernel.dudvsv + kernel.dudvsh;

% Choose kernel normalization
flds_kern = fieldnames(kernel);
if(ifnorm)
    for ii = 1:length(flds_kern)
        kernel.(flds_kern{ii}) = kernel.(flds_kern{ii}) ./ repmat(model(:,1),1,Nper);
    end
end
% Set inf values to --> 0
for ii = 1:length(flds_kern)
    kernel.(flds_kern{ii})(isinf(kernel.(flds_kern{ii}))) = 0;
end

% Convert gamma(Qs,Qp) kernels to --> Q(Qs,Qp) kernels
% gamma = omega / 2 / grv / Q
omega = 2*pi./dispersion.period(:)';
fac = repmat(2 .* dispersion.grv(:)' ./ omega, nlayer, 1);
kernel.dQdQs = fac .* kernel.dgdQs;
kernel.dQdQp = fac .* kernel.dgdQp;

% Convert gamma(Qs,Qp) kernels to --> gamma(Qmu, Qkappa) kernels
% Derived following eq. 8 of Wang and Cai (2017)
% Qs^-1 = Qmu^-1
% Qp^-1 = (1-4/3 * vs^2/vp^2)*Qkap^-1 + 4/3*(vs^2/vp^2)*Qmu^-1
% (Voigt-average velocities are used for the TI model)
model_iso = convert_TI2iso(model);
vp = repmat(model_iso(:,2),1,Nper);
vs = repmat(model_iso(:,3),1,Nper);
kernel.dgdQmu = kernel.dgdQs + kernel.dgdQp.*(4/3 .* vs.^2 ./ vp.^2);
kernel.dgdQkap = kernel.dgdQp .* (1 - 4/3 .* vs.^2 ./ vp.^2);

% Convert gamma(Qmu, Qkappa) --> Q(Qmu, Qkappa) kernel
kernel.dQdQmu = fac .* kernel.dgdQmu;
kernel.dQdQkap = fac .* kernel.dgdQkap;

% Convert alpha (gamma) --> Q
dispersion.qinv = 2 .* dispersion.grv .* dispersion.gamma ./ (2*pi./dispersion.period);

% Get mid point of layer
vec_h=model(:,1);
z1=cumsum(vec_h);z1=z1(:);
z2=[0;z1(1:end-1)];
z_mid=(z1+z2)/2;
kernel.z = z_mid;

if(ifplot)
    for ip = 1:Nper
        period = vec_T(ip);
        if(ifnorm)
            t=['Depth sensitivity kernel',num2str(period),' S'];
        else
            t=['Layer sensitivity kernel',num2str(period),' S'];
        end
        figure;hold on
        if(wavetype=='L')
            subplot(1,3,1);
            plot(kernel.dcdvsv(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVsv'); set(gca,'fontsize',20);
            subplot(1,3,2);
            plot(kernel.dcdvsh(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVsh'); set(gca,'fontsize',20);
            subplot(1,3,3);
            plot(kernel.dQdQmu(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('d{Q^{-1}}/dQ^{-1}_{\mu}'); set(gca,'fontsize',20);
        elseif(wavetype=='R')
            xl = {'dC/dVpv','dC/dVph','dC/dVsv','d{Q^{-1}}/dQ^{-1}_{\kappa}','d{Q^{-1}}/dQ^{-1}_{\mu}'};
            fl = {'dcdvpv','dcdvph','dcdvsv','dQdQkap','dQdQmu'};
            for jj = 1:5
                subplot(1,5,jj);
                plot(kernel.(fl{jj})(:,ip),z_mid,'-b'); axis ij
                title(t); xlabel(xl{jj}); set(gca,'fontsize',20);
            end
        end
    end
end

return
end
