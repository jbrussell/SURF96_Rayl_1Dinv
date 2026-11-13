% Two-stage inversion for PA5: can a Rayleigh-only isotropic Vsv model
% (which recovers the LVZ) serve as the starting model for the Rayleigh +
% Love inversion for Vsh / xi? (follow-up to
% test_kernel_inversion_vsvh_vpvh_rho_eta_PA5.m, same models, data and
% regularization)
%
% Three inversions of the same synthetic 10-150 s data (true model PA5,
% mode-0 Rayleigh and Love):
%   A. one stage: R + L from the smooth isotropic starting model (as
%      inversion 1 of test_kernel_inversion_vsvh_vpvh_rho_eta_PA5.m)
%   B. two stages: (1) Rayleigh only, xi held at 1 (eps_xi = 100), from the
%      smooth starting model -> isotropic model (Vsh = Vsv); (2) R + L from
%      the stage-1 model
%   C. best case for B: R + L starting from the true PA5 Vsv as an
%      isotropic model (Vsh = Vsv), with the smoothing broken at the PA5
%      68.5 and 166 km discontinuities, i.e. a perfect stage 1
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

nmode = 0; % mode branch number (0=fund)

nit = 8; % total number of iterations
nit_recalc_c = 1; % number of iterations after which to recalculate the phase velocities
nit_recalc_G = nit_recalc_c; % number of iterations after which to recalculate the kernels
nmodes = nmode+1; % number of modes computed by tdisp96 (more can help when spurious roots occur; see calc_der96_TI.m)
fref = 1; % [Hz] reference frequency of the Q model (only used if the models have Qp, Qs)

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


%% Data (mode 0, as in the PA5 test)
periods = [vec_T(:); vec_T(:)];
wavetype = [repmat('R',length(vec_T),1); repmat('L',length(vec_T),1)];
I_R = wavetype=='R'; I_L = wavetype=='L';
cobs = zeros(size(periods));
cobs(I_R) = dispR_surf96_TI(periods(I_R),truemod,nmode);
cobs(I_L) = dispL_surf96_TI(periods(I_L),truemod,nmode);
cstd = cobs * 0.01;

%% Starting models
% C: true Vsv as an isotropic model (PA5 Vpv/Vsv, Vph/Vsh, rho, eta)
start_C = pa5;
start_C(:,6) = pa5(:,3);                         % Vsh = Vsv
start_C(:,5) = pa5(:,3).*pa5(:,5)./max(pa5(:,6),eps); % Vph from PA5 Vph/Vsh
start_C(pa5(:,3)==0,5) = pa5(pa5(:,3)==0,5);     % water
discs_C = sort([discs(:); 68.54; 166.06]);

%% Inversions
names = {'A: one stage (R+L)','B: two stages (R, then R+L)','C: start = true V_{SV}, isotropic'};
res = struct([]);
for ic = 1:3
    fprintf('\n=== %s ===\n', names{ic}); tic;
    switch ic
        case 1
            m0 = startmod; d0 = discs;
        case 2
            % stage 1: Rayleigh only, xi held at the starting model (1)
            [m1, cpre1] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs(I_R),cstd(I_R),periods(I_R),'R',startmod,discs, ...
                eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,100,eps_eta,z_dampbot,nit,nit_recalc_c,nmode,0,vsv_thresh_xi,fref,nit_recalc_G,nmodes);
            fprintf('stage 1 (R only): chi2 = %.3f, max |xi-1| = %.4f\n', mean(((cobs(I_R)-cpre1)./cstd(I_R)).^2), ...
                max(abs((m1(m1(:,3)>0,6)./m1(m1(:,3)>0,3)).^2-1)));
            stage1 = m1;
            m0 = m1; d0 = discs;
        case 3
            m0 = start_C; d0 = discs_C;
    end
    [fm, cpre, unc] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,wavetype,m0,d0, ...
        eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,nit,nit_recalc_c,nmode,0,vsv_thresh_xi,fref,nit_recalc_G,nmodes);
    toc
    nl = size(fm,1); xf = (fm(:,6)./fm(:,3)).^2; xs = nan(nl,1);
    for k = 1:nl
        if fm(k,3) > 0
            Jx = [-2*xf(k)/fm(k,3), 2*xf(k)/fm(k,6)]; ik = [k, nl+k];
            xs(k) = sqrt(Jx*unc.Cm(ik,ik)*Jx');
        end
    end
    res(ic).start = m0; res(ic).finalmod = fm; res(ic).cpre = cpre; res(ic).xi = xf; res(ic).xi_std = xs;
    res(ic).chi2 = mean(((cobs-cpre)./cstd).^2);
end

%% Summary
xi_of = @(m) (m(:,6)./m(:,3)).^2;
zb = [z_moho 68.5; 68.5 100; 100 135; 135 166; 166 220; 220 260];
xt = xi_of(truemod);
fprintf('\nMean xi in depth bands\n%-36s', 'depth (km)'); fprintf('   %3.0f-%3.0f      ', zb'); fprintf('  chi2\n');
fprintf('%-36s', 'true (PA5)');
for ib = 1:size(zb,1); ii = zmid >= zb(ib,1) & zmid < zb(ib,2); fprintf('   %6.3f         ', mean(xt(ii))); end
fprintf('\n');
for ic = 1:3
    fprintf('%-36s', names{ic});
    for ib = 1:size(zb,1)
        ii = zmid >= zb(ib,1) & zmid < zb(ib,2);
        fprintf('   %6.3f+-%5.3f  ', mean(res(ic).xi(ii)), mean(res(ic).xi_std(ii)));
    end
    fprintf('  %.3f\n', res(ic).chi2);
end

%% Plot
figure(110); clf; set(gcf,'position',[50 100 1500 700],'color','w');
clrc = [0.85 0.1 0.1; 0.2 0.6 0.2; 0.5 0.2 0.8];
pv = {@(m) m(:,3), @(m) m(:,6), xi_of};
pl = {'V_{SV} (km/s)','V_{SH} (km/s)','\xi = (V_{SH}/V_{SV})^2'};
for k = 1:3
    subplot(1,4,k); box on; hold on; hh = [];
    hh(1) = plotlayermods(startmod(:,1), pv{k}(startmod), '-b'); set(hh(1),'linewidth',2);
    hh(2) = plotlayermods(truemod(:,1), pv{k}(truemod), '-k'); set(hh(2),'linewidth',3);
    hh(3) = plotlayermods(stage1(:,1), pv{k}(stage1), '--'); set(hh(3),'linewidth',2,'color',clrc(2,:));
    for ic = 1:3
        h = plotlayermods(res(ic).finalmod(:,1), pv{k}(res(ic).finalmod), '-'); set(h,'linewidth',2,'color',clrc(ic,:));
        hh(3+ic) = h;
        if k == 3
            h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi - res(ic).xi_std, ':'); set(h,'color',clrc(ic,:));
            h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi + res(ic).xi_std, ':'); set(h,'color',clrc(ic,:));
        end
    end
    xlabel(pl{k}); ylim([0 300]);
    if k == 1; ylabel('Depth (km)'); title('PA5: one vs two stages'); end
    if k < 3; xlim([3.9 5.0]); else; xlim([0.95 1.15]);
        legend(hh, [{'smooth start','true (PA5)','stage 1 (R only)'}, names], 'location', 'southeast', 'fontsize', 10); end
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
end
subplot(1,4,4); box on; hold on; plot([0 0],[0 300],'-k');
for ic = 1:3
    h = plotlayermods(res(ic).finalmod(:,1), res(ic).xi - xt, '-'); set(h,'linewidth',2,'color',clrc(ic,:));
end
xlabel('\xi_{final} - \xi_{true}'); ylim([0 300]); xlim([-0.1 0.1]);
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'PA5_two_stage.png']); end

%% Plot data fit (all inversions)
cstart = zeros(size(periods)); % smooth starting model (A, B stage 1)
cstart(I_R) = dispR_surf96_TI(periods(I_R),startmod,nmode);
cstart(I_L) = dispL_surf96_TI(periods(I_L),startmod,nmode);
cstage1 = zeros(size(periods)); % stage-1 model (start of B stage 2)
cstage1(I_R) = cpre1;
cstage1(I_L) = dispL_surf96_TI(periods(I_L),stage1,nmode);
figure(111); clf; set(gcf,'position',[300 300 1200 800],'color','w');
for iw = 1:length(data_types)
    wt = data_types(iw);
    I_w = wavetype == wt;
    subplot(2,length(data_types),iw); box on; hold on;
    hs = plot(periods(I_w),cstart(I_w),'-ob','linewidth',2);
    h1 = plot(periods(I_w),cstage1(I_w),'--o','color',clrc(2,:),'linewidth',2);
    hc = [];
    for ic = 1:length(res)
        hc(ic) = plot(periods(I_w),res(ic).cpre(I_w),'-o','color',clrc(ic,:),'linewidth',2);
    end
    he = errorbar(periods(I_w),cobs(I_w),2*cstd(I_w),'sk'); set(he,'markersize',6,'markerfacecolor','k','linewidth',1.5);
    if iw==1; legend([hs h1 hc he],[{'c start (smooth)','c stage 1 (R only)'}, names, {'c obs (PA5)'}],'Location','southeast','fontsize',10); end
    xlabel('Period (s)'); ylabel([wt,' phase velocity (km/s)']);
    set(gca,'FontSize',14,'linewidth',1.5,'xscale','log'); xlim([9 160]);
    subplot(2,length(data_types),length(data_types)+iw); box on; hold on;
    plot([9 160],[0 0],'-k');
    plot(periods(I_w),(cobs(I_w)-cstage1(I_w))./cstd(I_w),'--o','color',clrc(2,:),'linewidth',2);
    for ic = 1:length(res)
        plot(periods(I_w),(cobs(I_w)-res(ic).cpre(I_w))./cstd(I_w),'-o','color',clrc(ic,:),'linewidth',2);
    end
    xlabel('Period (s)'); ylabel([wt,' (c_{obs}-c_{pre})/\sigma']);
    set(gca,'FontSize',14,'linewidth',1.5,'xscale','log'); xlim([9 160]);
end
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'PA5_two_stage_data_fit.png']); end

%% Plot Rayleigh and Love kernels of the true model (PA5) at all periods
Npers = length(vec_T);
clr = jet(Npers);
lgd = arrayfun(@(T) sprintf('%.1f s',T), vec_T, 'uniformoutput', false);
for iw = 1:length(data_types)
    wt = data_types(iw);
    if wt == 'R'
        flds = {'dcdvsv','dcdvsh','dcdvpv','dcdvph','dcdrho','dcdeta'};
        xl = {'dc/dV_{SV}','dc/dV_{SH}','dc/dV_{PV}','dc/dV_{PH}','dc/d\rho','dc/d\eta'};
    else % Love waves are not sensitive to Vpv, Vph, eta
        flds = {'dcdvsv','dcdvsh','dcdrho'};
        xl = {'dc/dV_{SV}','dc/dV_{SH}','dc/d\rho'};
    end
    kernel = calc_kernel96_TI(truemod, vec_T, wt, 1, 0, nmode, fref, struct('phase_only',true,'nmodes',nmodes));
    figure(111+iw); clf;
    set(gcf,'position',[100 100 1500 450],'color','w')
    nf = length(flds);
    for k = 1:nf
        subplot(1,nf,k); box on; hold on;
        for ip = 1:Npers
            plot(kernel.(flds{k})(:,ip),kernel.z,'-','color',clr(ip,:),'linewidth',2);
        end
        xlabel(xl{k});
        if k==1; ylabel('Depth (km)'); title([wt,' kernels, true model (PA5)']); end
        ylim([0 400]);
        set(gca,'fontsize',13,'linewidth',1.5,'ydir','reverse');
    end
    pos = get(gca,'Position');
    legend(lgd,'location','eastoutside');
    set(gca,'position',pos);
    if ~isempty(figdir); print(gcf,'-dpng',[figdir,'PA5_kernels_true_',wt,'.png']); end
end
