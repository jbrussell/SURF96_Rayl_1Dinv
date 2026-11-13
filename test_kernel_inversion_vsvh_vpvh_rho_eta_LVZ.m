% Test inversion for a radially anisotropic (TI) oceanic model with a low
% velocity zone (LVZ) in the asthenosphere, using 10-150 s Rayleigh + Love
% phase velocities (extension of test_kernel_inversion_vsvh_vpvh_rho_eta.m)
%
% Question: the fundamental-mode Love-wave kernels peak within the LVZ, and
% where they peak depends on the LVZ itself (the LVZ acts as a waveguide, so
% the kernels are a nonlinear function of the model). Does this affect the
% recovery of the radial anisotropy xi = (Vsh/Vsv)^2?
%
% True model: the starting model (MINEOS card, isotropic, weak LVZ) with
%   - a fast lid (+2% at the Moho tapering to 0 at 70 km) and a strong LVZ
%     (Gaussian, -6% at 130 km, half-width ~45 km) in Vsv, Vsh, Vpv, Vph
%     and rho, so that Vsv drops from ~4.45 km/s in the lid to ~4.0 km/s
%   - radial anisotropy with Vsh > Vsv: xi = xi_max from the Moho to
%     z_ani1, tapering linearly to 1 at z_ani2 (Vph raised with Vsh, so
%     phi = (Vpv/Vph)^2 = 1/xi, as in the other TI test)
%
% Three inversions of synthetic data are compared (run_comparisons = 1):
%   1. LVZ true model, kernels recalculated every iteration (the default)
%   2. LVZ true model, kernels held fixed at the (LVZ-free, isotropic)
%      starting model, i.e. a single linearization
%   3. Control: same xi, but no lid/LVZ perturbation, kernels recalculated
% 1 vs 2 shows how much the nonlinearity of the kernels matters; 1 vs 3
% shows whether the LVZ itself degrades the recovery of xi.
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

figdir = './test_inversion_LVZ_figs/'; % output figures ('' = do not save)
if ~isempty(figdir) && ~exist(figdir,'dir'); mkdir(figdir); end

run_comparisons = 1; % 1 = also run inversions 2 and 3 (see above)

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

nmode = 0; % mode branch number (0=fund)

nit = 8; % total number of iterations
nit_recalc_c = 1; % number of iterations after which to recalculate the phase velocities
nit_recalc_G = nit_recalc_c; % number of iterations after which to recalculate the kernels
nmodes = nmode+1; % number of modes computed by tdisp96 (more can help when spurious roots occur; see calc_der96_TI.m)
fref = 1; % [Hz] reference frequency of the Q model (only used if the models have Qp, Qs)

%% True model parameters
dlnv_lid = 0.02;  z_lid = 70;           % fast lid: +2% at the Moho, 0 at z_lid
dlnv_lvz = -0.06; z_lvz = 130; w_lvz = 45; % LVZ: Gaussian, center z_lvz, half-width w_lvz
xi_max = 1.10; z_ani1 = 120; z_ani2 = 220; % xi = xi_max to z_ani1, 1 at z_ani2

%% Starting and true models
cardname = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_Line1_dist100.00km';
CARD = ['./MINEOS_CARD/',cardname,'.card'];
[startmod, discs] = card2mod_TI(CARD,zmax);
% eta ~= 1 in the half-space can give spurious Rayleigh roots in tdisp96
startmod(end,7) = 1;
z = [0; cumsum(startmod(:,1))];
zmid = z(1:end-1) + startmod(:,1)/2;
z_moho = max(discs(discs < 60)); % Moho (discs also has the 410 km etc.)
imantle = zmid > z_moho;

% velocity perturbation (lid + LVZ), applied to Vpv, Vsv, rho, Vph, Vsh
dlnv = zeros(size(zmid));
dlnv(imantle) = dlnv_lid*max(0, 1 - (zmid(imantle)-z_moho)/(z_lid-z_moho)) ...
              + dlnv_lvz*exp(-((zmid(imantle)-z_lvz)/w_lvz).^2);
% radial anisotropy
xi_z = ones(size(zmid));
xi_z(imantle) = 1 + (xi_max-1)*min(1, max(0, (z_ani2-zmid(imantle))/(z_ani2-z_ani1)));

truemod_iso = startmod; % lid + LVZ, no anisotropy
truemod_iso(:,2:6) = startmod(:,2:6).*(1+dlnv);
truemods = {}; % {LVZ + xi, no LVZ + xi}
for imod = 1:2
    if imod == 1; tm = truemod_iso; else; tm = startmod; end
    tm(:,6) = tm(:,6).*sqrt(xi_z); % Vsh
    tm(:,5) = tm(:,5).*sqrt(xi_z); % Vph
    truemods{imod} = tm;
end
truemod = truemods{1};

figure(1); clf; set(gcf,'position',[50 50 1200 600]);
plot_model_TI({startmod, truemod, truemods{2}}, {'start','true (LVZ)','true (no LVZ)'}, 300);
subplot(1,4,1); title('Starting and true models');
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'LVZ_models.png']); end

%% Kernels in the starting and true models: where do they peak?
% The nonlinearity: the kernels of the true model (LVZ + xi) differ from
% those of the starting model used in the first iteration.
T_k = [20 40 60 80 100 125 150];
kmods = {startmod, truemod_iso, truemod};
kname = {'start (no LVZ, \xi=1)','LVZ, \xi=1','true (LVZ, \xi>1)'};
kst = {':','--','-'};
K = cell(3,2); % {model, wave}
for im = 1:3
    for iw = 1:2
        K{im,iw} = calc_kernel96_TI(kmods{im}, vec_T, data_types(iw), 0, 0, nmode);
    end
end
[~, i_vsvmin] = min(truemod(zmid > z_moho & zmid < 300,3));
z_vsvmin = zmid(find(zmid > z_moho,1) + i_vsvmin - 1); % depth of the Vsv minimum
Hkm = startmod(:,1); Hkm(Hkm == 0) = NaN; % layer thickness (kernels per km; not the half-space)

figure(10); clf; set(gcf,'position',[50 50 1500 650],'color','w');
fk = {'dcdvsv','dcdvsh'}; % Rayleigh: dc/dVsv; Love: dc/dVsh
xk = {'R: dc/dV_{SV} per km (norm.)','L: dc/dV_{SH} per km (norm.)'};
clrT = jet(length(T_k));
for iw = 1:2
    subplot(1,4,iw); box on; hold on;
    hl = [];
    for ip = 1:length(T_k)
        [~,jp] = min(abs(vec_T - T_k(ip)));
        for im = [1 3]
            kk = K{im,iw}.(fk{iw})(:,jp)./Hkm; % per km
            kk = kk/max(kk(imantle)); % normalized by the mantle maximum
            h = plot(kk, zmid, kst{im}, 'color', clrT(ip,:), 'linewidth', 2);
            if im == 3; hl(ip) = h; end
        end
    end
    xl0 = [0 1.3]; plot(xl0,[1 1]*z_vsvmin,'-','color',[1 1 1]*0.7,'linewidth',2); xlim(xl0);
    xlabel(xk{iw}); ylabel('Depth (km)'); ylim([0 350]);
    if iw == 1
        title('dotted: start, solid: true');
        legend(hl, arrayfun(@(t) sprintf('%.0f s',vec_T(find(abs(vec_T-t)==min(abs(vec_T-t)),1))), T_k, 'uniformoutput', false), 'location', 'southeast');
    end
    set(gca,'fontsize',13,'linewidth',1.5,'ydir','reverse');
end
% kernel-weighted mean depth (centroid) of the mantle part of the kernel vs period
subplot(1,4,3); box on; hold on;
clrw = [0.1 0.3 0.8; 0.85 0.2 0.1];
hl = []; lg = {};
for iw = 1:2
    for im = 1:3
        kk = K{im,iw}.(fk{iw})./Hkm;
        kk(~imantle | isnan(Hkm),:) = 0;
        zc = (zmid(:)'*kk)./sum(kk,1);
        hl(end+1) = plot(vec_T, zc, [kst{im},'o'], 'color', clrw(iw,:), 'linewidth', 2, 'markerfacecolor', clrw(iw,:));
        lg{end+1} = [data_types(iw),', ',kname{im}];
    end
end
plot([10 150],[1 1]*z_vsvmin,'-','color',[1 1 1]*0.7,'linewidth',3);
set(gca,'xscale','log','ydir','reverse','fontsize',13,'linewidth',1.5);
xlim([10 150]); ylim([0 250]);
xlabel('Period (s)'); ylabel('Kernel centroid depth, mantle (km)');
legend(hl, lg, 'location', 'southwest', 'fontsize', 9);
title('gray: V_{SV} minimum (true)');
% sensitivity to xi at fixed Voigt Vs: dc/dln(xi) per km.
%   Vsv^2 = 3 Vs^2/(2+xi), Vsh^2 = 3 xi Vs^2/(2+xi), at xi = xi_true:
%   dVsv/dxi = -Vsv/(2(2+xi)), dVsh/dxi = Vsh/(xi(2+xi)) - Vsh/(2(2+xi))
subplot(1,4,4); box on; hold on;
for iw = 1:2
    for im = [1 3]
        m = kmods{im}; xx = (m(:,6)./m(:,3)).^2; xx(~isfinite(xx)) = 1;
        dvsv = -m(:,3)./(2*(2+xx)); dvsh = m(:,6)./(xx.*(2+xx)) - m(:,6)./(2*(2+xx));
        % (Vph scaled with Vsh: dVph/dxi = Vph/Vsh dVsh/dxi)
        r_ph = m(:,5)./m(:,6); r_ph(~isfinite(r_ph)) = 0;
        dvph = r_ph.*dvsh;
        [~,jp] = min(abs(vec_T - 100));
        kx = xx.*(K{im,iw}.dcdvsv(:,jp).*dvsv + K{im,iw}.dcdvsh(:,jp).*dvsh + K{im,iw}.dcdvph(:,jp).*dvph)./Hkm;
        plot(kx, zmid, kst{im}, 'color', clrw(iw,:), 'linewidth', 2);
    end
end
xl0 = xlim; plot(xl0,[1 1]*z_vsvmin,'-','color',[1 1 1]*0.7,'linewidth',2); xlim(xl0);
xlabel('dc/dln\xi per km (100 s)'); ylim([0 350]);
set(gca,'xtick',[-0.006 -0.003 0 0.003 0.006]);
title('blue R, red L; dotted start, solid true');
set(gca,'fontsize',13,'linewidth',1.5,'ydir','reverse');
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'LVZ_kernels_start_vs_true.png']); end

%% Synthetic data and inversions
periods = []; wavetype = '';
for wt = data_types
    periods = [periods; vec_T(:)];
    wavetype = [wavetype; repmat(wt, length(vec_T), 1)];
end
I_R = wavetype=='R'; I_L = wavetype=='L';

cstart = zeros(size(periods));
cstart(I_R) = dispR_surf96_TI(periods(I_R),startmod,nmode);
cstart(I_L) = dispL_surf96_TI(periods(I_L),startmod,nmode);

cases = struct('name', {'LVZ, kernels updated', 'LVZ, kernels fixed (start)', 'no LVZ, kernels updated'}, ...
               'itrue', {1, 1, 2}, 'nit_recalc_c', {nit_recalc_c, nit+1, nit_recalc_c}, ...
               'nit_recalc_G', {nit_recalc_G, nit+1, nit_recalc_G});
if ~run_comparisons; cases = cases(1); end
res = struct([]);
for ic = 1:length(cases)
    tm = truemods{cases(ic).itrue};
    cobs = zeros(size(periods));
    cobs(I_R) = dispR_surf96_TI(periods(I_R),tm,nmode); % "observations"
    cobs(I_L) = dispL_surf96_TI(periods(I_L),tm,nmode);
    cstd = cobs * 0.01; % observation uncertainties
    fprintf('\n=== Inversion %d: %s ===\n', ic, cases(ic).name);
    tic;
    [fm, cpre, unc] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,wavetype,startmod,discs, ...
        eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,nit,cases(ic).nit_recalc_c,nmode,ic==1,vsv_thresh_xi,fref,cases(ic).nit_recalc_G,nmodes);
    toc
    if ic == 1 && ~isempty(figdir); figure(2); print(gcf,'-dpng',[figdir,'LVZ_iterations.png']); end
    % xi uncertainty from the model covariance: xi = Vsh^2/Vsv^2
    nl = size(fm,1);
    xf = (fm(:,6)./fm(:,3)).^2;
    xs = nan(nl,1);
    for k = 1:nl
        if fm(k,3) > 0
            Jx = [-2*xf(k)/fm(k,3), 2*xf(k)/fm(k,6)];
            ik = [k, nl+k];
            xs(k) = sqrt(Jx*unc.Cm(ik,ik)*Jx');
        end
    end
    res(ic).truemod = tm; res(ic).finalmod = fm; res(ic).cobs = cobs; res(ic).cstd = cstd;
    res(ic).cpre = cpre; res(ic).unc = unc; res(ic).xi = xf; res(ic).xi_std = xs;
    res(ic).chi2 = mean(((cobs-cpre)./cstd).^2);
    res(ic).chi2_R = mean(((cobs(I_R)-cpre(I_R))./cstd(I_R)).^2);
    res(ic).chi2_L = mean(((cobs(I_L)-cpre(I_L))./cstd(I_L)).^2);
end
finalmod = res(1).finalmod; cobs = res(1).cobs; cstd = res(1).cstd; cpre = res(1).cpre; unc = res(1).unc;

%% Summary: xi in depth bands
zb = [z_moho 60; 60 100; 100 140; 140 180; 180 220; 220 260];
xi_of = @(m) (m(:,6)./m(:,3)).^2;
fprintf('\nMean xi in depth bands (true / recovered +- 1 sigma)\n');
fprintf('%-28s', 'depth (km)'); fprintf('   %3.0f-%3.0f      ', zb'); fprintf('  chi2 (R, L)\n');
for ic = 1:length(res)
    xt = xi_of(res(ic).truemod);
    fprintf('%-28s', ['true (',cases(ic).name(1:strfind(cases(ic).name,',')-1),')']);
    for ib = 1:size(zb,1)
        ii = zmid >= zb(ib,1) & zmid < zb(ib,2);
        fprintf('   %6.3f         ', mean(xt(ii)));
    end
    fprintf('\n%-28s', cases(ic).name);
    for ib = 1:size(zb,1)
        ii = zmid >= zb(ib,1) & zmid < zb(ib,2);
        fprintf('   %6.3f+-%5.3f  ', mean(res(ic).xi(ii)), mean(res(ic).xi_std(ii)));
    end
    fprintf('  %.2f (%.2f, %.2f)\n', res(ic).chi2, res(ic).chi2_R, res(ic).chi2_L);
end

%% Plot final kernels (inversion 1)
lgd = {};
Npers = length(vec_T);
clr = jet(Npers);
flds = {'dcdvsv','dcdvsh','dcdvpv','dcdvph','dcdrho','dcdeta'};
xl = {'dc/dVsv','dc/dVsh','dc/dVpv','dc/dVph','dc/d\rho','dc/d\eta'};
for iw = 1:length(data_types)
    wt = data_types(iw);
    kernel = calc_kernel96_TI(finalmod, vec_T, wt, 0, 0, nmode);
    figure(100+iw); clf;
    set(gcf,'position',[100 100 1500 450],'color','w')
    for k = 1:6
        subplot(1,6,k); box on; hold on;
        for ip = 1:Npers
            plot(kernel.(flds{k})(:,ip),kernel.z,'-','color',clr(ip,:),'linewidth',2);
            lgd{ip} = [num2str(vec_T(ip),'%.1f'),' s'];
        end
        xlabel(xl{k});
        if k==1; ylabel('Depth (km)'); title([wt,' kernels, final model']); end
        ylim([0 400]);
        set(gca,'fontsize',13,'linewidth',1.5,'ydir','reverse');
    end
    pos = get(gca,'Position');
    legend(lgd,'location','eastoutside');
    set(gca,'position',pos);
    if ~isempty(figdir); print(gcf,'-dpng',[figdir,'LVZ_final_kernels_',wt,'.png']); end
end

%% Plot data fit (all inversions)
figure(103); clf; set(gcf,'position',[300 300 1200 800]);
clrc = [0.85 0.1 0.1; 0.95 0.6 0.1; 0.2 0.6 0.2];
for iw = 1:length(data_types)
    wt = data_types(iw);
    I_w = wavetype == wt;
    subplot(2,length(data_types),iw); box on; hold on;
    hs = plot(periods(I_w),cstart(I_w),'-ob','linewidth',2);
    hc = [];
    for ic = 1:length(res)
        hc(ic) = plot(periods(I_w),res(ic).cpre(I_w),'-o','color',clrc(ic,:),'linewidth',2);
    end
    he = errorbar(periods(I_w),cobs(I_w),2*cstd(I_w),'sk'); set(he,'markersize',6,'markerfacecolor','k','linewidth',1.5);
    if run_comparisons
        plot(periods(I_w),res(3).cobs(I_w),'sk','markersize',6,'linewidth',1.5);
    end
    if iw==1; legend([hs hc he],[{'c start'}, {cases(1:length(res)).name}, {'c obs (LVZ)'}],'Location','southeast','fontsize',10); end
    xlabel('Period (s)'); ylabel([wt,' phase velocity (km/s)']);
    set(gca,'FontSize',14,'linewidth',1.5,'xscale','log'); xlim([9 160]);
    subplot(2,length(data_types),length(data_types)+iw); box on; hold on;
    plot([9 160],[0 0],'-k');
    for ic = 1:length(res)
        plot(periods(I_w),(res(ic).cobs(I_w)-res(ic).cpre(I_w))./res(ic).cstd(I_w),'-o','color',clrc(ic,:),'linewidth',2);
    end
    xlabel('Period (s)'); ylabel([wt,' (c_{obs}-c_{pre})/\sigma']);
    set(gca,'FontSize',14,'linewidth',1.5,'xscale','log'); xlim([9 160]);
end
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'LVZ_data_fit.png']); end

%% Plot true and recovered models (inversion 1): Vsv, Vsh, Vpv, Vph, rho, eta, xi
vals = {@(m) m(:,3), @(m) m(:,6), @(m) m(:,2), @(m) m(:,5), @(m) m(:,4), @(m) m(:,7), xi_of};
stds = {unc.vsv_std, unc.vsh_std, unc.vpv_std, unc.vph_std, unc.rho_std, unc.eta_std, res(1).xi_std};
labels = {'V_{SV} (km/s)','V_{SH} (km/s)','V_{PV} (km/s)','V_{PH} (km/s)','\rho (g/cm^3)','\eta','\xi = (V_{SH}/V_{SV})^2'};
figure(100); clf;
set(gcf,'position',[50 100 1700 650]);
for k = 1:7
    subplot(1,7,k); box on; hold on;
    v_true = vals{k}(truemod); v_start = vals{k}(startmod); v_final = vals{k}(finalmod);
    h1 = plotlayermods(truemod(:,1), v_true, '-k'); set(h1,'linewidth',2);
    h2 = plotlayermods(startmod(:,1), v_start, '-b'); set(h2,'linewidth',2);
    h3 = plotlayermods(finalmod(:,1), v_final, '-r'); set(h3,'linewidth',2);
    h4 = plotlayermods(finalmod(:,1), v_final - stds{k}, '--r'); set(h4,'linewidth',1);
    h5 = plotlayermods(finalmod(:,1), v_final + stds{k}, '--r'); set(h5,'linewidth',1);
    xlabel(labels{k});
    if k==1; ylabel('Depth (km)'); title(['Data: ',data_types,', ',num2str(min(vec_T)),'-',num2str(max(vec_T)),' s']); end
    ylim([0 z_dampbot+50]);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
    if k==1; legend([h1 h2 h3],{'true','start','final'},'Location','southwest'); end
end
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'LVZ_recovered_model.png']); end

%% Compare Vsv, Vsh and xi for the three inversions
if run_comparisons
    figure(104); clf; set(gcf,'position',[50 100 1500 700],'color','w');
    pv = {@(m) m(:,3), @(m) m(:,6), xi_of};
    pl = {'V_{SV} (km/s)','V_{SH} (km/s)','\xi = (V_{SH}/V_{SV})^2'};
    for k = 1:3
        subplot(1,4,k); box on; hold on;
        hh = [];
        hh(1) = plotlayermods(startmod(:,1), pv{k}(startmod), '-b'); set(hh(1),'linewidth',2);
        hh(2) = plotlayermods(truemods{1}(:,1), pv{k}(truemods{1}), '-k'); set(hh(2),'linewidth',3);
        hh(3) = plotlayermods(truemods{2}(:,1), pv{k}(truemods{2}), '--k'); set(hh(3),'linewidth',2);
        for ic = 1:length(res)
            h = plotlayermods(res(ic).finalmod(:,1), pv{k}(res(ic).finalmod), '-'); set(h,'linewidth',2,'color',clrc(ic,:));
            hh(3+ic) = h;
            if k == 3
                h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi - res(ic).xi_std, ':'); set(h,'linewidth',1,'color',clrc(ic,:));
                h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi + res(ic).xi_std, ':'); set(h,'linewidth',1,'color',clrc(ic,:));
            end
        end
        xlabel(pl{k}); ylim([0 300]);
        if k == 1; ylabel('Depth (km)'); xlim([3.9 4.8]); end
        if k == 2; xlim([3.9 4.8]); end
        if k == 3
            xlim([0.95 1.2]);
            legend(hh, [{'start','true (LVZ)','true (no LVZ)'}, {cases.name}], 'location', 'southeast', 'fontsize', 10);
        end
        set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
    end
    % xi error
    subplot(1,4,4); box on; hold on;
    plot([0 0],[0 300],'-k');
    for ic = 1:length(res)
        h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi - xi_of(res(ic).truemod), '-'); set(h,'linewidth',2,'color',clrc(ic,:));
    end
    xlabel('\xi_{final} - \xi_{true}'); ylim([0 300]); xlim([-0.08 0.08]);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
    if ~isempty(figdir); print(gcf,'-dpng',[figdir,'LVZ_xi_comparison.png']); end
end
