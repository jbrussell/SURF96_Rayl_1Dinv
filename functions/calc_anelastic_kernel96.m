function [kernel, dispersion] = calc_anelastic_kernel96(model, vec_T, wavetype, ifnorm, ifplot, varargin)
% calculate anelastic (Q-corrected) surface wave sensitivity kernels
%
% [kernel, dispersion] = calc_anelastic_kernel96(model, period, wavetype, ifNorm, ifplot);
% [kernel, dispersion] = calc_anelastic_kernel96(model, period, wavetype, ifNorm, ifplot, nmode, fref);
% model: surf96 format model: [thickness, vp, vs, rho, Qp, Qs];
% period : in seconds
% wavetype: 'L' or 'R'
% ifnorm: 1: normalize kernel by layer thickness; 0: no normalization
% ifplot: whether to plot the kernel
% nmode: 0 for fund. mode; 1 for 1st overtone, etc. (OPTIONAL, default 0)
% fref: reference frequency (Hz) of the constant-Q model (OPTIONAL, default 1)
%
% OUTPUT
% kernel: structure of [nlayer x nperiod] kernels
%   dcdvs, dcdvp, dudvs, dudvp     causal velocity kernels
%   dcdrho, dudrho                 density kernels (elastic)
%   dcdh, dudh                     srfker96 layer-thickness kernels (elastic)
%   dcdQs, dcdQp, dudQs, dudQp, dgdQs, dgdQp   kernels w.r.t. Qs^-1, Qp^-1
%   dQdQs, dQdQp, dgdQmu, dgdQkap, dQdQmu, dQdQkap
%   z                              layer mid-point depth
% dispersion: structure with fields period, mode, phv, grv, gamma, qinv,
%   phv_elastic, grv_elastic
%
% Yang Zha 11/06/2013
%
% jbrussell 6/5/2020: Note, if using the kernels for velocity inversion,
% set ifnorm=0. Kernels are scaled such that input data
% is in units of km/s and the output model is velocity perurbation in
% units of km/s (not percent!).
%
% jbrussell 11/20/2022: Include density (rho) kernels in addition to Vs and Vp.
% Also update to allow for overtone kernels to be output
% nmode=0 for fund. mode; nmode=1 for 1st overtone, etc...
%
% jbrussell 6/12/2023: Significant rewrite to include anelastic kernels. This version
% reads in a new structure called "dispersion" that includes attenuation coefficient
% gamma as well as phv and grv with the influence of physical dispersion.
%
% jbrussell 11/2026: The anelastic quantities used to be read from the
% "Anelastic" section of the srfker96 output. srfker96 derives them from
% the phase-velocity kernels, which srfdrr96/srfdrl96 output without the
% spherical-earth mapping when only group velocity is requested (see
% calc_kernel96.m), so gamma, the Q kernels and the causal velocity kernels
% inherited that error. They are now computed here from the corrected
% elastic kernels of calc_kernel96.m with the same formulas as srfker96
% (subroutines getgam and cslmod):
%   gamma = w/(2 c0^2) * sum[ dc/dVs Vs Qs^-1 + dc/dVp Vp Qp^-1 ]
%   c     = c0 + 2 gamma c0^2 ln(f/fref) / (pi w)
%   U     = U0 (1 + (2 - U0/c0)(c-c0)/c0 + 2 gamma U0/(pi w))
% Other changes: fref is now used as the reference frequency of the Q model
% (srfker96 used the surf96 default); dcdrho/dudrho are true density
% kernels (previously srfker96's dc/dh, dU/dh, now dcdh/dudh); Love waves
% are supported. As in srfker96, the density and thickness kernels are the
% elastic values.
%

if size(model,2)<6
    error('Must include Qp and Qs in model to generate anelastic kernels');
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
if strcmp(wavetype,'J')
    error('calc_anelastic_kernel96: joint (J) kernels are not supported; call separately for R and L');
end

vec_T = vec_T(:)';
Nper = length(vec_T);
nlayer = size(model,1);

% Elastic kernels (not normalized), with the spherical-earth correction
[dcdvs, dcdvp, dudvs, dudvp, z_mid, dcdrho, dudrho, dcdh, dudh, phv0, grv0] = ...
    calc_kernel96(model(:,1:4), vec_T, wavetype, 0, 0, nmode);
if isempty(dcdvs)
    disp('Something went wrong calculating kernels');
    kernel = []; dispersion = [];
    return
end

% Q^-1 (CPS convention: Q>1 is converted to 1/Q)
vp = model(:,2); vs = model(:,3);
qa = model(:,5); qb = model(:,6);
qa(qa>1) = 1./qa(qa>1);
qb(qb>1) = 1./qb(qb>1);
qb(vs<=0.01) = 0; % fluid

kernel.dcdvs = nan(nlayer,Nper); kernel.dcdvp = nan(nlayer,Nper);
kernel.dudvs = nan(nlayer,Nper); kernel.dudvp = nan(nlayer,Nper);
kernel.dcdrho = dcdrho; kernel.dudrho = dudrho;
kernel.dcdh = dcdh; kernel.dudh = dudh;
qflds = {'dcdQs','dcdQp','dudQs','dudQp','dgdQs','dgdQp'};
for ii = 1:length(qflds)
    kernel.(qflds{ii}) = nan(nlayer,Nper);
end

dispersion.period = vec_T(:);
dispersion.mode = nmode;
dispersion.phv = nan(Nper,1);
dispersion.grv = nan(Nper,1);
dispersion.gamma = nan(Nper,1);
dispersion.phv_elastic = phv0(:);
dispersion.grv_elastic = grv0(:);

for ip = 1:Nper
    c0 = phv0(ip);
    u0 = grv0(ip);
    if isnan(c0) || isnan(u0)
        continue
    end
    f = 1/vec_T(ip);
    omega = 2*pi*f;
    faclog = log(f/fref)/pi;

    dcds = dcdvs(:,ip); dcdp = dcdvp(:,ip);
    duds = dudvs(:,ip); dudp = dudvp(:,ip);
    if strcmp(wavetype,'L')
        dcdp = zeros(nlayer,1); dudp = zeros(nlayer,1);
    end
    sumS = dcds.*vs; sumP = dcdp.*vp;
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
    kernel.dcdvs(:,ip) = dcds.*(1 + faclog*qb);
    kernel.dcdvp(:,ip) = dcdp.*(1 + faclog*qa);
    kernel.dudvs(:,ip) = duds.*(uu0 - u0c0*cmc0c0 + facg) ...
        + dcds.*u0c0.*(-2*facg + u0c0*qb/pi + (2-u0c0)*(faclog*qb - cmc0c0) + u0c0*cmc0c0);
    kernel.dudvp(:,ip) = dudp.*(uu0 - u0c0*cmc0c0 + facg) ...
        + dcdp.*u0c0.*(-2*facg + u0c0*qa/pi + (2-u0c0)*(faclog*qa - cmc0c0) + u0c0*cmc0c0);

    % Q^-1 kernels
    kernel.dcdQs(:,ip) = faclog*sumS;
    kernel.dcdQp(:,ip) = faclog*sumP;
    kernel.dgdQs(:,ip) = 0.5*omega*sumS/(c0*c0);
    kernel.dgdQp(:,ip) = 0.5*omega*sumP/(c0*c0);
    kernel.dudQs(:,ip) = u0c0*(2-u0c0)*kernel.dcdQs(:,ip) + u0c0*u0c0*sumS/pi;
    kernel.dudQp(:,ip) = u0c0*(2-u0c0)*kernel.dcdQp(:,ip) + u0c0*u0c0*sumP/pi;
end

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
% Derived following eq. 8 of Wang and Cai (2017) "A Method to Estimate Shear Quality Factor of Hard Rocks"
% Qs^-1 = Qmu^-1
% Qp^-1 = (1-4/3 * vs^2/vp^2)*Qkap^-1 + 4/3*(vs^2/vp^2)*Qmu^-1
VP = repmat(vp,1,Nper);
VS = repmat(vs,1,Nper);
kernel.dgdQmu = kernel.dgdQs + kernel.dgdQp.*(4/3 .* VS.^2 ./ VP.^2);
kernel.dgdQkap = kernel.dgdQp .* (1 - 4/3 .* VS.^2 ./ VP.^2);

% Convert gamma(Qmu, Qkappa) --> Q(Qmu, Qkappa) kernel
kernel.dQdQmu = fac .* kernel.dgdQmu;
kernel.dQdQkap = fac .* kernel.dgdQkap;

% Convert alpha (gamma) --> Q
dispersion.qinv = 2 .* dispersion.grv .* dispersion.gamma ./ (2*pi./dispersion.period);

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
            subplot(1,2,1);
            plot(kernel.dcdvs(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVs'); set(gca,'fontsize',20);
            subplot(1,2,2);
            plot(kernel.dQdQmu(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('d{Q^{-1}}/dQ^{-1}_{\mu}'); set(gca,'fontsize',20);
        elseif(wavetype=='R')
            subplot(1,4,1);
            plot(kernel.dcdvp(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVp'); set(gca,'fontsize',20);
            subplot(1,4,2);
            plot(kernel.dcdvs(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVs'); set(gca,'fontsize',20);
            subplot(1,4,3);
            plot(kernel.dQdQkap(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('d{Q^{-1}}/dQ^{-1}_{\kappa}'); set(gca,'fontsize',20);
            subplot(1,4,4);
            plot(kernel.dQdQmu(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('d{Q^{-1}}/dQ^{-1}_{\mu}'); set(gca,'fontsize',20);
        end
    end
end

return
end
