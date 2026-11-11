function [kernel, dispersion] = calc_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot, varargin)
% Calculate surface wave sensitivity kernels for a transversely isotropic
% (TI) model using the CPS TI programs (TI analog of calc_kernel96.m)
%
% [kernel, dispersion] = calc_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot);
% [kernel, dispersion] = calc_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot, nmode);
% [kernel, dispersion] = calc_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot, nmode, fref);
%
% model: TI model [H Vpv Vsv rho Vph Vsh eta]  (see writemod_surf96_TI.m;
%        use convert_iso2TI.m to convert an isotropic [H Vp Vs rho] model)
%        If Qp and Qs are also provided ([... Qp Qs], 9 columns), the
%        Q-corrected kernels are calculated with calc_anelastic_kernel96_TI.m
%        (same behavior as calc_kernel96.m).
% vec_T: period in seconds
% wavetype: 'L' or 'R'  (joint 'J' is not supported)
% ifnorm: 1: normalize kernel by layer thickness; 0: no normalization
% ifplot: whether to plot the kernel
% nmode: 0 for fund. mode; 1 for 1st overtone, etc. (OPTIONAL, default 0)
% fref: reference frequency (Hz) for Q (OPTIONAL, default 1; only used for anelastic)
%
% OUTPUT
% kernel: structure of [nlayer x nperiod] layer kernels. Partial
%   derivatives are taken with respect to Vpv, Vph, Vsv, Vsh, eta and rho,
%   each holding the others fixed (this is the CPS parameterization).
%     kernel.dcdvpv, kernel.dcdvph, kernel.dcdvsv, kernel.dcdvsh,
%     kernel.dcdeta, kernel.dcdrho, kernel.dcdh
%     kernel.dudvpv, kernel.dudvph, kernel.dudvsv, kernel.dudvsh,
%     kernel.dudeta, kernel.dudrho, kernel.dudh
%   Isotropic-equivalent kernels (perturbing Vpv & Vph, or Vsv & Vsh,
%   together at fixed eta); for an isotropic model these are the
%   dc/dVp, dc/dVs, dU/dVp, dU/dVs kernels of calc_kernel96.m:
%     kernel.dcdvp = dcdvpv + dcdvph,  kernel.dcdvs = dcdvsv + dcdvsh
%     kernel.dudvp = dudvpv + dudvph,  kernel.dudvs = dudvsv + dudvsh
%   kernel.z: layer mid-point depth
% dispersion: structure with fields period, mode, phv, grv
%
% NOTES
%  - Phase-velocity kernels come directly from tregn96/tlegn96 -DER.
%  - CPS does not compute group-velocity kernels for TI media (there is no
%    TI version of srfker96), so they are computed from the phase-velocity
%    kernels following Rodi et al. (1975):
%        dU/dm = (U/c)(2 - U/c) dc/dm + (U^2 w / c^2) d(dc/dm)/dw
%    where d/dw is evaluated by centered finite difference using two
%    additional frequencies w(1 +/- dw_frac) with dw_frac = 2e-2
%    (see calc_der96_TI.m, which also rejects spurious tdisp96 roots).
%  - kernel.dcdh is the CPS TI thickness partial (DC/DH in tdpder96); it
%    agrees with the srfker96 dc/dh returned by calc_kernel96.m (outputs 8-9).
%  - kernel.dcdrho is a density kernel, comparable with outputs 6-7 of
%    calc_kernel96.m.
%  - Love waves have no sensitivity to Vpv, Vph, or eta (returned as zeros).
%
% jbrussell 11/2026
%

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
    error('calc_kernel96_TI: joint (J) kernels are not supported; call separately for R and L');
end

if size(model,2)>7 % if Qp and Qs provided, calculate the Q-corrected kernels

    [kernel, dispersion] = calc_anelastic_kernel96_TI(model, vec_T, wavetype, ifnorm, ifplot, nmode, fref);
    return

end

vec_T = vec_T(:)';
Nper = length(vec_T);
nlayer = size(model,1);

% Elastic calculation: Q^-1 = 0 in every layer. calc_der96_TI also returns
% the same mode branch at w(1+dw_frac) and w(1-dw_frac)
model_e = [model(:,1:7), zeros(nlayer,2)];
try
    [eig, eig_p, eig_m, info] = calc_der96_TI(vec_T, model_e, wavetype, nmode, fref);
catch
    eig = [];
end
if isempty(eig)
    disp('Something went wrong calculating kernels');
    kernel = []; dispersion = [];
    return
end
dw_frac = info.dw_frac;

params = {'vpv','vph','vsv','vsh','eta','rho','h'};
dispersion.period = vec_T(:);
dispersion.mode = nmode;
dispersion.phv = eig.phv(:);
dispersion.grv = eig.grv(:);
for ii = 1:length(params)
    kernel.(['dcd',params{ii}]) = eig.(['dcd',params{ii}]);
    kernel.(['dud',params{ii}]) = nan(nlayer,Nper);
end

for ip = 1:Nper
    c = eig.phv(ip);
    U = eig.grv(ip);
    if isnan(c)
        continue
    end
    for ii = 1:length(params)
        dcdm = eig.(['dcd',params{ii}])(:,ip);
        Kp = eig_p.(['dcd',params{ii}])(:,ip);
        Km = eig_m.(['dcd',params{ii}])(:,ip);
        if ~any(isnan(Kp)) && ~any(isnan(Km))
            wdKdw = (Kp - Km) / (2*dw_frac);   % w * d(dc/dm)/dw, centered
        elseif ~any(isnan(Kp))
            wdKdw = (Kp - dcdm) / dw_frac;     % one-sided (e.g. near a cutoff)
        elseif ~any(isnan(Km))
            wdKdw = (dcdm - Km) / dw_frac;
        else
            continue
        end
        kernel.(['dud',params{ii}])(:,ip) = (U/c)*(2-U/c)*dcdm + (U^2/c^2)*wdKdw;
    end
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
        % plotting sensitivity
        figure;hold on
        if(wavetype=='L')
            subplot(1,2,1);
            plot(kernel.dcdvsv(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVsv'); set(gca,'fontsize',20);
            subplot(1,2,2);
            plot(kernel.dcdvsh(:,ip),z_mid,'-b'); axis ij
            title(t); xlabel('dC/dVsh'); set(gca,'fontsize',20);
        elseif(wavetype=='R')
            xl = {'dC/dVpv','dC/dVph','dC/dVsv','dC/d\eta'};
            fl = {'dcdvpv','dcdvph','dcdvsv','dcdeta'};
            for jj = 1:4
                subplot(1,4,jj);
                plot(kernel.(fl{jj})(:,ip),z_mid,'-b'); axis ij
                title(t); xlabel(xl{jj}); set(gca,'fontsize',20);
            end
        end
    end
end

return
end
