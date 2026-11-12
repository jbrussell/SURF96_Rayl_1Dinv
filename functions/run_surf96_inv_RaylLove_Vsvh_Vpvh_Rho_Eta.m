function [finalmod,cpre,unc] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,wavetype,startmod,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,nit,nit_recalc_c,nmode,varargin)
% [finalmod,cpre,unc] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(..., nmode, isfigure, vsv_thresh_xi, fref, nit_recalc_G, nmodes)
% Do linearized inversion of Rayleigh and/or Love wave phase velocities for
% a transversely isotropic (radially anisotropic) model: Vsv, Vsh, Vpv,
% Vph, density (rho) and eta, using the CPS TI programs (via the *_TI
% wrappers) to calculate the kernels and phase velocities.
%
% The inversion solves for 6 parameters in each layer: Vsv, Vsh, Vpv, Vph,
% Rho, and Eta. Constraint equations include:
% - Norm damping of Vsv and Vsh toward the starting model (increasing
%   below z_dampbot)
% - First derivative (flatness) smoothing of Vsv and Vsh
% - Second derivative (roughness) smoothing of Vsv and Vsh
% - Damping Vpv/Vsv, Vph/Vsh and Rho/Vsv toward the starting-model ratios
% - Damping xi = (Vsh/Vsv)^2 toward the starting model (weight eps_xi):
%   use a large eps_xi for Rayleigh-only or Love-only data, which cannot
%   separate Vsv and Vsh
% - Damping eta toward the starting model (weight eps_eta)
% - (optional) xi held at the starting model in the crust and sediments,
%   i.e. layers with 0 < Vsv < vsv_thresh_xi (otherwise the strong
%   sensitivity of the shallow layers tends to absorb the Rayleigh-Love
%   discrepancy as unrealistically large shallow anisotropy)
% - Water layer(s) fixed to the starting model
%
% Thus, although we technically solve for all 6 parameters, we solve
% primarily for Vsv and Vsh and scale Vpv, Vph and Rho accordingly. The
% kernels for all 6 parameters are included, so that the scaled parameters
% contribute to the predicted phase velocities.
%
% INPUTS (N=number of data; M=number of layers in model)
% cobs - phase velocity observed km/s [N x 1]
% cstd - phase velocity uncertainty km/s [N x 1]
% periods - seconds [N x 1]
% wavetype - 'R' or 'L' for each datum [N x 1 char], e.g. ['R';'R';'L'],
%            or a single character if all data are the same type.
%            Rayleigh and Love data may be mixed in any order.
% startmod - starting TI model [H Vpv Vsv rho Vph Vsh eta] [M x 7]
%            (km, km/s, g/cm^3; see writemod_surf96_TI.m), or
%            [H Vpv Vsv rho Vph Vsh eta Qp Qs] [M x 9] to include
%            attenuation in the forward calculation (see below)
% discs - depth of desired discontinuities in smoothing (km) [any length x 1]
% eps_data - weight for data fit [scalar]
% eps_H - weight for norm damping as fraction of G norm [scalar]
% eps_J - weight for 1st derivative smoothing as fraction of G norm [scalar]
% eps_F - weight for 2nd derivative smoothing as fraction of G norm [scalar]
% eps_vpvs - weight for enforcing Vpv/Vsv and Vph/Vsh as fraction of G norm [scalar]
% eps_rhovs - weight for enforcing rho/Vsv as fraction of G norm [scalar]
% eps_xi - weight for damping xi to the starting model as fraction of G norm [scalar]
% eps_eta - weight for damping eta to the starting model as fraction of G norm [scalar]
% z_dampbot - depth below which to damp to starting model (km) [scalar]
% nit - number of iterations [scalar]
% nit_recalc_c - number of iterations after which to recalculate the phase
%                velocities (otherwise the linearized prediction G*dm is
%                used) [scalar]
% nmode - mode number (0 = fundamental) [scalar]
% isfigure - plot the progress of the inversion (OPTIONAL, default 1)
% vsv_thresh_xi - hold xi at the starting model where 0 < Vsv < vsv_thresh_xi
%                 (km/s) in the starting model (OPTIONAL, default 0 = off)
% fref - reference frequency (Hz) of the Q model (OPTIONAL, default 1)
% nit_recalc_G - number of iterations after which to recalculate the
%                kernels (OPTIONAL, default nit_recalc_c)
% nmodes - number of modes computed by tdisp96 at each period (OPTIONAL,
%          default nmode+1; see calc_der96_TI.m)
%
% COST: the phase velocities and kernels are computed with the CPS TI
% programs; tdisp96 (Rayleigh waves especially) takes nearly all the time.
% - When the kernels are recalculated, the phase velocities come from the
%   same calculation (calc_kernel96_TI with opts.phase_only: 2 frequencies
%   per period, no group-velocity kernels).
% - When only the phase velocities are recalculated, the linearized
%   prediction is used to pick the root at each period (calc_der96_TI with
%   opts.c_guess: 1 frequency per period, with a fallback to the full root
%   validation where no root is within 2% of the prediction).
%
% OUTPUTS
% finalmod - final TI model [H Vpv Vsv rho Vph Vsh eta] [M x 7]
% cpre - phase velocity predicted for the final model (km/s) [N x 1]
% unc - formal uncertainties (structure): vsv_std, vsh_std, vpv_std,
%       vph_std, rho_std, eta_std [M x 1], the model covariance Cm and
%       model resolution R [6M x 6M] (parameter order Vsv, Vsh, Vpv, Vph,
%       Rho, Eta) and data resolution Ndat [N x N]
%
% ATTENUATION: if startmod has Qp and Qs columns, the phase velocities and
% kernels include the causal (physical dispersion) correction for that Q
% model (dispR/dispL_surf96_TI and calc_kernel96_TI with Q, reference
% frequency fref). Q is not inverted for: Qp and Qs are held at their
% starting values. Without Q columns the forward problem is elastic
% (Q^-1 = 0).
%
% jbrussell - 11/2026

isfigure = 1;
if ~isempty(varargin) && ~isempty(varargin{1})
    isfigure = varargin{1};
end
vsv_thresh_xi = 0;
if length(varargin) >= 2 && ~isempty(varargin{2})
    vsv_thresh_xi = varargin{2};
end
fref = 1;
if length(varargin) >= 3 && ~isempty(varargin{3})
    fref = varargin{3};
end
nit_recalc_G = nit_recalc_c;
if length(varargin) >= 4 && ~isempty(varargin{4})
    nit_recalc_G = varargin{4};
end
nmodes = nmode + 1; % default: only the modes needed (see calc_der96_TI.m)
if length(varargin) >= 5 && ~isempty(varargin{5})
    nmodes = varargin{5};
end
if ~any(size(startmod,2) == [7 9])
    error('startmod must have 7 columns [H Vpv Vsv rho Vph Vsh eta] or 9 columns [... Qp Qs]');
end

cobs = cobs(:); cstd = cstd(:); periods = periods(:);
N = length(cobs);
wavetype = upper(wavetype(:));
if length(wavetype) == 1
    wavetype = repmat(wavetype, N, 1);
end
if length(wavetype) ~= N
    error('wavetype must have one entry per datum (or be a single character)');
end

eps_large = 1e4; % large weight (relative to the G norm) to force constraint equations

nlayer = size(startmod,1);
z = cumsum(startmod(:,1)); % bottom of each layer (km)
zmid = z - startmod(:,1)/2; zmid(end) = z(end-1);

% Kernels for the G matrix
% kernels and phase velocities of the starting model (one calculation)
[G, cstart] = build_G(startmod, periods, wavetype, nmode, fref, nmodes);

% Data weighting
W = diag(1./cstd);

% Starting-model ratios
vsv0 = startmod(:,3); vsh0 = startmod(:,6);
vpv_vsv = startmod(:,2) ./ vsv0;
vph_vsh = startmod(:,5) ./ vsh0;
rho_vsv = startmod(:,4) ./ vsv0;
xi0 = vsh0.^2 ./ vsv0.^2;
ind_h2o = find(vsv0 == 0);
ind_solid = find(vsv0 > 0);
xi0(ind_h2o) = 1;

% Set up smoothing and damping matrices (applied to Vsv and Vsh)
dz = gradient(zmid);
dz_mat = repmat(dz,1,length(dz));
% Damping matrix
H00 = eye(nlayer);
% first derivative flatness
J00 = build_flatness(nlayer) ./ dz_mat;
% second derivative smoothing
F00 = build_smooth( nlayer ) ./ dz_mat.^2;

% Break constraints at discontinuities
z_brks = sort(discs);
J00 = break_constraint(J00, z, z_brks);
F00 = break_constraint(F00, z, z_brks);

% Rescale the constraint weights relative to the data kernels
NA = norm(W*G,1);
eps_H0 = eps_H*NA/norm(H00,1);
eps_J0 = eps_J*NA/norm(J00,1);
eps_F0 = eps_F*NA/norm(F00,1);

% Damp towards starting model below z_dampbot
hv = vsv0; hh = vsh0;
ind_dampstart = find(z > z_dampbot);
ramp = linspace(1,1000,length(ind_dampstart));
H00(ind_dampstart,ind_dampstart) = H00(ind_dampstart,ind_dampstart).*ramp;
hv(ind_dampstart) = hv(ind_dampstart).*ramp';
hh(ind_dampstart) = hh(ind_dampstart).*ramp';

% Column blocks of the model vector m = [Vsv; Vsh; Vpv; Vph; Rho; Eta]
Z = zeros(nlayer);
I = eye(nlayer);
blk = @(k) [repmat(Z,1,k-1), I, repmat(Z,1,6-k)];

% Damping and smoothing on Vsv (block 1) and Vsh (block 2)
Hd = [H00*blk(1); H00*blk(2)];  hd = [hv; hh];
Jd = [J00*blk(1); J00*blk(2)];  jd = zeros(2*size(J00,1),1);
Fd = [F00*blk(1); F00*blk(2)];  fd = zeros(2*size(F00,1),1);

% Ratio constraints (solid layers only)
%       Vpv - (Vpv/Vsv)0 Vsv = 0,  Vph - (Vph/Vsh)0 Vsh = 0,  Rho - (Rho/Vsv)0 Vsv = 0
S = I(ind_solid,:);
VPV_mat = S*blk(3) - S*diag(fixnan(vpv_vsv))*blk(1);
VPH_mat = S*blk(4) - S*diag(fixnan(vph_vsh))*blk(2);
RHO_mat = S*blk(5) - S*diag(fixnan(rho_vsv))*blk(1);
eps_vpvs0 = eps_vpvs*NA/norm(VPV_mat,1);
eps_rhovs0 = eps_rhovs*NA/norm(RHO_mat,1);

% xi constraint (solid layers):  Vsh - sqrt(xi0) Vsv = 0
XI_mat = S*blk(2) - S*diag(sqrt(xi0))*blk(1);
eps_xi0 = eps_xi*NA/norm(XI_mat,1);
% hold xi at the starting model in the crust/sediments (Vsv < vsv_thresh_xi)
w_xi = ones(length(ind_solid),1);
w_xi(vsv0(ind_solid) < vsv_thresh_xi) = eps_large*NA/eps_xi0;
XI_mat = diag(w_xi)*XI_mat;

% eta damping toward the starting model
ETA_mat = blk(6);
eta_vec = startmod(:,7);
eps_eta0 = eps_eta*NA/norm(ETA_mat,1);

% Water layer(s): all parameters fixed to the starting model
W_h2o = zeros(6*length(ind_h2o), 6*nlayer); w_h2o = zeros(6*length(ind_h2o),1);
m0 = model2vec(startmod);
irow = 0;
for k = 1:6
    for ii = ind_h2o(:)'
        irow = irow + 1;
        W_h2o(irow, (k-1)*nlayer + ii) = 1;
        w_h2o(irow) = m0((k-1)*nlayer + ii);
    end
end

% combine all constraints
H = [Hd*eps_H0; Jd*eps_J0; Fd*eps_F0; VPV_mat*eps_vpvs0; VPH_mat*eps_vpvs0; RHO_mat*eps_rhovs0; ...
     XI_mat*eps_xi0; ETA_mat*eps_eta0; W_h2o*eps_large*NA];
h = [hd*eps_H0; jd*eps_J0; fd*eps_F0; zeros(3*length(ind_solid),1); ...
     zeros(length(ind_solid),1); eta_vec*eps_eta0; w_h2o*eps_large*NA];

% Data vector
if any(isnan(cstart))
    error('run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta: %d phase velocities of the starting model could not be calculated', sum(isnan(cstart)));
end
cpre = cstart;

% Least squares inversion
premod = startmod;
m_pre = model2vec(premod);
clrs = jet(nit);
for ii = 1:nit

    % reformulate inverse problem so that constraints apply directly to model
    % G(m-m0) = dobs-d0  -->  Gm = G*m0 + dobs-d0
    dc = cobs - cpre;
    d = dc + G*m_pre;

    % least squares
    F = [W*G*eps_data; H];
    f = [W*d*eps_data; h];
    m = (F'*F)\(F'*f);
    m(1:2*nlayer) = max(m(1:2*nlayer), 0); % no negative shear velocities

    % update data vector
    dm = m - m_pre;
    cpre = cpre + G*dm;

    % update model
    m_pre = m;
    premod = vec2model(m_pre, premod);

    % Model uncertainties and resolution
    FtF = F'*F;
    Ginv = eps_data^2 * (FtF \ (G'*(W'*W))); % generalized inverse [6M x N]
    unc.R = Ginv * G;                        % model resolution [6M x 6M]
    unc.Ndat = G * Ginv;                     % data resolution [N x N]
    unc.Cm = Ginv * diag(cstd.^2) * Ginv';   % model covariance from data errors
    m_std = sqrt(diag(unc.Cm));
    unc.vsv_std = m_std(1:nlayer);
    unc.vsh_std = m_std(nlayer+1:2*nlayer);
    unc.vpv_std = m_std(2*nlayer+1:3*nlayer);
    unc.vph_std = m_std(3*nlayer+1:4*nlayer);
    unc.rho_std = m_std(4*nlayer+1:5*nlayer);
    unc.eta_std = m_std(5*nlayer+1:6*nlayer);

    do_c = mod(ii,nit_recalc_c)==0 && ii < nit; % recalculate phase velocities
    do_G = mod(ii,nit_recalc_G)==0 && ii < nit; % recalculate kernels
    if do_G
        % kernels (and phase velocities, from the same calculation)
        [G, cnew] = build_G(premod, periods, wavetype, nmode, fref, nmodes);
    elseif do_c
        % phase velocities only, using the linearized prediction to pick roots
        cnew = predict(premod, periods, wavetype, nmode, fref, nmodes, cpre);
    end
    if do_c
        inan = isnan(cnew);
        if any(inan)
            warning('run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta: %d phase velocities could not be calculated; keeping the linearized prediction there', sum(inan));
            cnew(inan) = cpre(inan);
        end
        cpre = cnew;
    end

    if isfigure
        if ii == 1
            figure(2); clf; set(gcf,'position',[370 372 1100 580]);
            subplot(2,3,[1 4]); box on; hold on;
            stairs_mod(startmod(:,1), startmod(:,3), '-k', 4);
            subplot(2,3,[2 5]); box on; hold on;
            stairs_mod(startmod(:,1), startmod(:,6), '-k', 4);
        end
        subplot(2,3,[1 4]);
        stairs_mod(premod(:,1), premod(:,3), clrs(ii,:), 2);
        xlabel('Vsv (km/s)'); ylabel('Depth (km)');
        title('Iterations (black: start)');
        set(gca,'FontSize',15,'linewidth',1.5,'ydir','reverse');
        subplot(2,3,[2 5]);
        stairs_mod(premod(:,1), premod(:,6), clrs(ii,:), 2);
        xlabel('Vsh (km/s)'); ylabel('Depth (km)');
        set(gca,'FontSize',15,'linewidth',1.5,'ydir','reverse');
        for iw = 1:2
            wt = 'RL'; I_w = wavetype == wt(iw);
            if ~any(I_w); continue; end
            subplot(2,3,3*iw); box on; hold on;
            if ii == 1
                plot(periods(I_w),cstart(I_w),'-ok','linewidth',3);
            end
            plot(periods(I_w),cpre(I_w),'-o','color',clrs(ii,:),'linewidth',2);
            he = errorbar(periods(I_w),cobs(I_w),2*cstd(I_w),'or'); set(he,'linewidth',2);
            xlabel('Period (s)'); ylabel([wt(iw),' phase velocity (km/s)']);
            set(gca,'FontSize',15,'linewidth',1.5);
        end
        drawnow
    end
end

finalmod = premod;
cpre = predict(finalmod, periods, wavetype, nmode, fref, nmodes, cpre);

end

%% ------------------------------------------------------------------------
function [G, c] = build_G(model, periods, wavetype, nmode, fref, nmodes)
% phase-velocity kernels [N x 6M] for m = [Vsv; Vsh; Vpv; Vph; Rho; Eta]
% and the phase velocities c [N x 1] from the same calculation
% (Q-corrected if the model has Qp, Qs columns)
N = length(periods);
nlayer = size(model,1);
G = zeros(N, 6*nlayer);
c = nan(N,1);
kopts.phase_only = true; % no group-velocity kernels: 2 frequencies per period
kopts.nmodes = nmodes;
flds = {'dcdvsv','dcdvsh','dcdvpv','dcdvph','dcdrho','dcdeta'};
for wt = 'RL'
    I_w = find(wavetype == wt);
    if isempty(I_w); continue; end
    [kernel, dispersion] = calc_kernel96_TI(model, periods(I_w), wt, 0, 0, nmode, fref, kopts);
    if isempty(kernel)
        error('run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta: kernel calculation failed');
    end
    c(I_w) = dispersion.phv(:);
    for k = 1:6
        K = kernel.(flds{k});
        K(isnan(K)) = 0; % mode does not exist at this period
        G(I_w, (k-1)*nlayer+(1:nlayer)) = K';
    end
end
end

function c = predict(model, periods, wavetype, nmode, fref, nmodes, c_guess)
% phase velocities for the TI model (Q-corrected if the model has Qp, Qs).
% c_guess (expected phase velocities, e.g. the linearized prediction) is
% used to pick the root at each period (see calc_der96_TI.m, opts.c_guess)
c = nan(length(periods),1);
for wt = 'RL'
    I_w = find(wavetype == wt);
    if isempty(I_w); continue; end
    dopts.nmodes = nmodes;
    dopts.c_guess = c_guess(I_w);
    if wt == 'R'
        c(I_w) = dispR_surf96_TI(periods(I_w), model, nmode, 'C', fref, dopts);
    else
        c(I_w) = dispL_surf96_TI(periods(I_w), model, nmode, 'C', fref, dopts);
    end
end
end

function m = model2vec(model)
% [H Vpv Vsv rho Vph Vsh eta] -> [Vsv; Vsh; Vpv; Vph; Rho; Eta]
m = [model(:,3); model(:,6); model(:,2); model(:,5); model(:,4); model(:,7)];
end

function model = vec2model(m, model)
n = size(model,1);
model(:,3) = m(1:n);
model(:,6) = m(n+1:2*n);
model(:,2) = m(2*n+1:3*n);
model(:,5) = m(3*n+1:4*n);
model(:,4) = m(4*n+1:5*n);
model(:,7) = m(5*n+1:6*n);
end

function x = fixnan(x)
x(~isfinite(x)) = 0;
end

function stairs_mod(H, v, clr, lw)
% plot a layered model as steps
ztop = [0; cumsum(H(1:end-1))];
zbot = [ztop(2:end); ztop(end) + max(20, 0.1*ztop(end))];
zz = reshape([ztop zbot]', [], 1);
vv = reshape([v v]', [], 1);
if ischar(clr)
    plot(vv, zz, clr, 'linewidth', lw);
else
    plot(vv, zz, '-', 'color', clr, 'linewidth', lw);
end
end
