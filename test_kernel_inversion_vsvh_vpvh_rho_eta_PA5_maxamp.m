% PA5 inversion with the Love (and Rayleigh) branch that dominates the
% surface records: nmode = 'maxamp' (see calc_der96_TI.m)
%
% In PA5 the slowest (mode-0) Love wave at 10-55 s is trapped in the LVZ
% and barely reaches the sea floor; the Love wave that dominates
% seismograms recorded at the sea floor (from a source near the sea floor)
% is mode 1, guided by the crust and lid. Mode 0 dominates above ~60 s.
% Here the synthetic data AND the forward calculations of the inversion use
% the branch with the largest surface amplitude at each period
% (nmode = 'maxamp'), i.e. what would be measured as the "fundamental
% mode".
%
% CAVEAT: in PA5 the selected branch is unambiguous only at ~10-12 s and
% above ~70 s. At ~13-55 s the crust/lid Love wave is faster than the LVZ
% Vsh and leaks into the LVZ, and its energy is shared by modes 1, 2 and 3:
% their surface amplitudes are within a factor of ~1-2.5 of each other and
% their group velocities within ~0.1 km/s, so a sea-floor record is a
% superposition of them. At 40-80 s mode 0 is also comparable. The
% phase velocity of the strongest branch is then an idealization of what
% would be measured (and it jumps from mode 1 to mode 0 near 60 s). The
% amplitude ratio of the strongest to the next-strongest branch in the true
% model is printed and plotted.
%
% Inversions (same models and regularization as
% test_kernel_inversion_vsvh_vpvh_rho_eta_PA5.m), 10-150 s R + L:
%   A. one stage (R+L) from the smooth starting model
%   B. two stages (Rayleigh only with xi = 1, then R+L)
%
% jbrussell 11/2026
%
clear
path2BIN = getenv('CPS_BIN'); % (optional) override the path to the CPS binaries
if isempty(path2BIN)
    path2BIN = './bin_v3.30/'; % path to surf96 binary
    % Make binary files executable
    system('chmod +x ./bin_v3.30/*');
end
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
    setenv('PATH', [path2BIN,':',PATH]);
end
addpath('./functions/')

figdir = './test_inversion_PA5_figs/'; % output figures ('' = do not save)
if ~isempty(figdir) && ~exist(figdir,'dir'); mkdir(figdir); end


%% Data to invert
data_types = 'RL'; % 'RL' = Rayleigh + Love (needed to constrain xi)
vec_T = logspace(log10(10),log10(150),20); % periods (s)

%% Inversion regularization parameters (as in test_kernel_inversion_vsvh_vpvh_rho_eta.m)
eps_data = 1; % data fit
eps_H = 0.005; % norm damping
eps_J = 0.001; % first derivative smoothing
eps_F = 2; % second derivative smoothing
eps_vpvs = 10; % enforce Vpv/Vsv and Vph/Vsh ratios
eps_rhovs = 10; % enforce rho/Vsv ratio
eps_eta = 10; % damp eta to starting model
eps_xi = 0.01; % damp xi=(Vsh/Vsv)^2 to starting model (weakly: R + L constrain xi)
vsv_thresh_xi = 4; % [km/s] hold xi at the starting model where Vsv < this ("crust")

% 150 s Rayleigh waves are sensitive down to ~400 km, so the model and
% the undamped depth range are deeper than in the 10-40 s test
z_dampbot = 250; % [km] damp to starting model below this depth
zmax = 600; % [km] maximum depth of model space

nmode = 'maxamp'; % branch with the largest surface amplitude among modes 0-2 (see calc_der96_TI.m)

nit = 8; % total number of iterations
nit_recalc_c = 1; % recalculate phase velocities and kernels every iteration

%% Starting and true models
CARD = './MINEOS_CARD/pa5_5km.card';
z_start = 250; % [km] starting model = PA5 below this depth
[pa5, discs_card] = card2mod_TI(CARD,zmax);
pa5(end,7) = 1; % eta ~= 1 in the half-space can give spurious Rayleigh roots in tdisp96
z = [0; cumsum(pa5(:,1))];
zmid = z(1:end-1) + pa5(:,1)/2;
z_moho = max(discs_card(discs_card < 30)); % Moho (12 km below sea level)
imantle = zmid > z_moho;
% smoothing breaks: only discontinuities that are in the starting model
discs = discs_card(discs_card <= z_moho | discs_card > z_dampbot);

vsv = pa5(:,3); vsh = pa5(:,6);
vs_voigt = sqrt((2*vsv.^2 + vsh.^2)/3);
vp_voigt = sqrt((pa5(:,2).^2 + 4*pa5(:,5).^2)/5);
istart = imantle & zmid < z_start;
i1 = find(istart,1); i2 = find(zmid >= z_start,1);
vs0 = vs_voigt;
vs0(istart) = interp1(zmid([i1 i2]), vs_voigt([i1 i2]), zmid(istart));

startmod = pa5;
startmod(istart,3) = vs0(istart); % Vsv
startmod(istart,6) = vs0(istart); % Vsh
startmod(istart,2) = vs0(istart).*pa5(istart,2)./vsv(istart); % Vpv (PA5 Vpv/Vsv)
startmod(istart,5) = vs0(istart).*pa5(istart,5)./vsh(istart); % Vph (PA5 Vph/Vsh)
startmod(istart,4) = vs0(istart).*pa5(istart,4)./vsv(istart); % rho (PA5 rho/Vsv)

% PA5 Voigt-isotropic (LVZ, xi = 1), used for the kernel comparison
truemod_iso = pa5;
truemod_iso(imantle,[3 6]) = repmat(vs_voigt(imantle),1,2);
truemod_iso(imantle,[2 5]) = repmat(vp_voigt(imantle),1,2);
% control: starting model + PA5 xi (Vsh and Vph raised by sqrt(xi))
xi_z = (vsh./vsv).^2; xi_z(~isfinite(xi_z)) = 1;
control = startmod;
control(:,6) = startmod(:,6).*sqrt(xi_z);
control(:,5) = startmod(:,5).*sqrt(xi_z);
truemods = {pa5, control};
truemod = pa5;


%% Data: branch with the largest surface amplitude
periods = [vec_T(:); vec_T(:)];
wavetype = [repmat('R',length(vec_T),1); repmat('L',length(vec_T),1)];
I_R = wavetype=='R'; I_L = wavetype=='L';
cobs = zeros(size(periods)); msel = zeros(size(periods));
eR = calc_eigenfunctions96_TI(vec_T, truemod, 'R', nmode);
eL = calc_eigenfunctions96_TI(vec_T, truemod, 'L', nmode);
cobs(I_R) = eR.phv; cobs(I_L) = eL.phv; msel(I_R) = eR.mode; msel(I_L) = eL.mode;
cstd = cobs * 0.01;
% surface amplitudes of modes 0-2 in the true model -> ambiguity ratio
isol = find(truemod(:,3) > 0, 1);
ampL = nan(3, length(vec_T)); c_modes = nan(3, length(vec_T));
for n = 0:2
    e = calc_eigenfunctions96_TI(vec_T, truemod, 'L', n);
    ampL(n+1,:) = e.amp .* e.ut(isol,:).^2;
    c_modes(n+1,:) = e.phv;
end
a0 = ampL; a0(isnan(a0)) = 0;
as = sort(a0, 1, 'descend');
ratioL = as(1,:)./as(2,:);
fprintf('\nLove (true model): T, selected mode, c, surface amplitude of modes 0-2 / max, ratio best/next\n');
disp([vec_T(:) eL.mode(:) eL.phv(:) (a0./max(a0,[],1))' ratioL(:)]);

%% Inversions
names = {'A: max. amplitude branch, one stage', 'B: max. amplitude branch, two stages'};
ncase = length(names);
res = struct([]);
for ic = 1:ncase
    fprintf('\n=== %s ===\n', names{ic}); tic;
    id = true(size(periods));
    m0 = startmod;
    if ic == 2
        % stage 1: Rayleigh only, xi held at the starting model (1)
        [m0, cpre1] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs(I_R),cstd(I_R),periods(I_R),'R',startmod,discs, ...
            eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,100,eps_eta,z_dampbot,nit,nit_recalc_c,nmode,0,vsv_thresh_xi);
        stage1 = m0;
    end
    [fm, cpre, unc] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs(id),cstd(id),periods(id),wavetype(id),m0,discs, ...
        eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,nit,nit_recalc_c,nmode,0,vsv_thresh_xi);
    toc
    nl = size(fm,1); xf = (fm(:,6)./fm(:,3)).^2; xs = nan(nl,1);
    for k = 1:nl
        if fm(k,3) > 0
            Jx = [-2*xf(k)/fm(k,3), 2*xf(k)/fm(k,6)]; ik = [k, nl+k];
            xs(k) = sqrt(Jx*unc.Cm(ik,ik)*Jx');
        end
    end
    % branch selected in the final model at all periods
    efL = calc_eigenfunctions96_TI(vec_T, fm, 'L', nmode);
    res(ic).finalmod = fm; res(ic).xi = xf; res(ic).xi_std = xs; res(ic).id = id;
    res(ic).cpre = cpre; res(ic).chi2 = mean(((cobs(id)-cpre)./cstd(id)).^2);
    res(ic).modeL = efL.mode; res(ic).cL = efL.phv;
    fprintf('Love branch selected in the final model: %s\n', mat2str(efL.mode));
    fprintf('Love branch selected in the true model:  %s\n', mat2str(eL.mode));
end

%% Summary
xi_of = @(m) (m(:,6)./m(:,3)).^2;
zb = [z_moho 68.5; 68.5 100; 100 135; 135 166; 166 220; 220 260];
xt = xi_of(truemod);
fprintf('\nMean xi in depth bands\n%-36s', 'depth (km)'); fprintf('   %3.0f-%3.0f      ', zb'); fprintf('  chi2\n');
fprintf('%-36s', 'true (PA5)');
for ib = 1:size(zb,1); ii = zmid >= zb(ib,1) & zmid < zb(ib,2); fprintf('   %6.3f         ', mean(xt(ii))); end
fprintf('\n');
for ic = 1:ncase
    fprintf('%-36s', names{ic}(1:min(36,end)));
    for ib = 1:size(zb,1)
        ii = zmid >= zb(ib,1) & zmid < zb(ib,2);
        fprintf('   %6.3f+-%5.3f  ', mean(res(ic).xi(ii)), mean(res(ic).xi_std(ii)));
    end
    fprintf('  %.3f\n', res(ic).chi2);
end

%% Plots
clrc = [0.85 0.1 0.1; 0.2 0.6 0.2];
% data: Love branches and the selected one
figure(120); clf; set(gcf,'position',[50 100 1300 550],'color','w');
subplot(1,2,1); box on; hold on;
hm = plot(vec_T, c_modes', '-', 'color', [1 1 1]*0.6, 'linewidth', 1.5);
h1 = plot(vec_T, eL.phv, 'sk', 'markerfacecolor', 'k', 'markersize', 8);
hc = [];
for ic = 1:ncase
    hc(ic) = plot(vec_T, res(ic).cL, '-o', 'color', clrc(ic,:), 'linewidth', 1.5, 'markersize', 4);
end
set(gca,'xscale','log','fontsize',13,'linewidth',1.5); xlim([9 160]);
xlabel('Period (s)'); ylabel('Love phase velocity (km/s)');
legend([hm(1) h1 hc], [{'PA5 modes 0-2','data (max. amplitude branch)'}, strcat(names, ' (final)')], 'location', 'northwest', 'fontsize', 9);
title('PA5 Love: max-surface-amplitude branch');
subplot(1,2,2); box on; hold on;
plot(vec_T, (a0./max(a0,[],1))', '-o', 'linewidth', 2);
set(gca,'xscale','log','yscale','log','fontsize',13,'linewidth',1.5); xlim([9 160]); ylim([1e-4 2]);
xlabel('Period (s)'); ylabel('surface amplitude / max'); legend({'mode 0','mode 1','mode 2'}, 'location', 'southwest');
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'PA5_maxamp_data.png']); end

figure(121); clf; set(gcf,'position',[50 100 1500 700],'color','w');
pv = {@(m) m(:,3), @(m) m(:,6), xi_of};
pl = {'V_{SV} (km/s)','V_{SH} (km/s)','\xi = (V_{SH}/V_{SV})^2'};
for k = 1:3
    subplot(1,4,k); box on; hold on; hh = [];
    hh(1) = plotlayermods(startmod(:,1), pv{k}(startmod), '-b'); set(hh(1),'linewidth',2);
    hh(2) = plotlayermods(truemod(:,1), pv{k}(truemod), '-k'); set(hh(2),'linewidth',3);
    for ic = 1:ncase
        h = plotlayermods(res(ic).finalmod(:,1), pv{k}(res(ic).finalmod), '-'); set(h,'linewidth',2,'color',clrc(ic,:));
        hh(2+ic) = h;
        if k == 3
            h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi - res(ic).xi_std, ':'); set(h,'color',clrc(ic,:));
            h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi + res(ic).xi_std, ':'); set(h,'color',clrc(ic,:));
        end
    end
    xlabel(pl{k}); ylim([0 300]);
    if k == 1; ylabel('Depth (km)'); title('PA5, max-amplitude branch'); end
    if k < 3; xlim([3.9 5.0]); else; xlim([0.95 1.15]);
        legend(hh, [{'smooth start','true (PA5)'}, names], 'location', 'southeast', 'fontsize', 9); end
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
end
subplot(1,4,4); box on; hold on; plot([0 0],[0 300],'-k');
for ic = 1:ncase
    h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi - xt, '-'); set(h,'linewidth',2,'color',clrc(ic,:));
end
xlabel('\xi_{final} - \xi_{true}'); ylim([0 300]); xlim([-0.1 0.1]);
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'PA5_maxamp_xi.png']); end
