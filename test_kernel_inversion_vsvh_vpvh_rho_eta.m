% Test inversion for a radially anisotropic (transversely isotropic, TI)
% model using the CPS TI programs (TI analog of
% test_kernel_inversion_vs_vp_rho.m)
%
% The inversion (run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta.m) solves for
% Vsv, Vsh, Vpv, Vph, Rho, and Eta in each layer, enforcing Vpv/Vsv,
% Vph/Vsh, and Rho/Vsv of the starting model and damping xi = (Vsh/Vsv)^2
% and eta toward the starting model. The kernels include all 6 parameters,
% so each iteration accounts for perturbations in all of them.
%
% Rayleigh + Love, Rayleigh-only, or Love-only phase velocities can be
% inverted (data_types below). Rayleigh waves are mostly sensitive to Vsv
% and Love waves to Vsh, so xi is only constrained by the data when both
% are used; with a single wave type xi is held close to the starting model
% (eps_xi large).
%
% Synthetic test: the true model is the starting model with Vsv, Vsh, Vpv,
% Vph, and Rho increased by 5% above 80 km and radial anisotropy added
% between 30 and 150 km by raising Vsh and Vph by sqrt(1.08)
% (xi = (Vsh/Vsv)^2 = 1.08, phi = (Vpv/Vph)^2 = 1/1.08).
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

%% Data to invert
data_types = 'RL'; % 'RL' = Rayleigh + Love, 'R' = Rayleigh only, 'L' = Love only

%% Inversion regularization parameters
eps_data = 1; % data fit
eps_H = 0.005; % norm damping
% Note: the norm damping pulls Vsv and Vsh toward their absolute starting
% values. Love waves need a larger change in Vsh (velocity change + xi)
% than Rayleigh waves need in Vsv, so a damping as strong as in
% test_kernel_inversion_vs_vp_rho.m (0.05) leaves the long-period Love
% data systematically underpredicted and xi underestimated.
eps_J = 0.001; % first derivative smoothing
eps_F = 2; % second derivative smoothing
eps_vpvs = 10; % enforce Vpv/Vsv and Vph/Vsh ratios
eps_rhovs = 10; % enforce rho/Vsv ratio
eps_eta = 10; % damp eta to starting model
if length(data_types) == 2
    eps_xi = 0.01; % damp xi=(Vsh/Vsv)^2 to starting model (weakly: R + L constrain xi)
else
    eps_xi = 100; % with a single wave type, xi is not constrained by the data
end

vsv_thresh_xi = 4; % [km/s] hold xi at the starting model where Vsv < this ("crust"); 0 = off

z_dampbot = 150; % [km] damp to starting model below this depth
zmax = 400; % [km] maximum depth of model space

nmode = 0; % mode branch number (0=fund, 1=1st overtone, 2=2nd overtone, etc.)

% Other inversion parameters
nit = 8; % total number of iterations
nit_recalc_c = 1; % number of iterations after which to recalculate phase velocity and kernels

%% Generate the synthetic dataset for this test
% Read in MINEOS model and convert to TI SURF96 layered model format
cardname = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_Line1_dist100.00km';

% Get MINEOS model [H Vpv Vsv rho Vph Vsh eta]
CARD = ['./MINEOS_CARD/',cardname,'.card'];
[truemod, discs] = card2mod_TI(CARD,zmax); % true model
% eta ~= 1 in the half-space can give spurious Rayleigh roots in tdisp96
% (see check_halfspace_eta_TI.m)
truemod(end,7) = 1;

% perturb true model
z = [0; cumsum(truemod(:,1))];
zmid = z(1:end-1) + truemod(:,1)/2;
Ibot = find(z<=80); Ibot = Ibot(end);
truemod(2:Ibot,[2 3 4 5 6]) = truemod(2:Ibot,[2 3 4 5 6])*1.05;
% add radial anisotropy between 30 and 150 km: Vsh and Vph are both raised
% by the factor sqrt(xi_true), so that xi = (Vsh/Vsv)^2 = xi_true and
% phi = (Vpv/Vph)^2 = 1/xi_true. (The inversion ties Vph to Vsh, so the P
% anisotropy of the true model is consistent with that parameterization.)
xi_true = 1.08;
I_ani = find(zmid>=30 & zmid<=150);
truemod(I_ani,6) = truemod(I_ani,6)*sqrt(xi_true); % Vsh
truemod(I_ani,5) = truemod(I_ani,5)*sqrt(xi_true); % Vph

% GENERATE SYNTHETIC DATASET
% Calculate dispersion for the true model (which is like the data)
vec_T = logspace(log10(10),log10(40),15);
periods = []; wavetype = '';
for wt = data_types
    periods = [periods; vec_T(:)];
    wavetype = [wavetype; repmat(wt, length(vec_T), 1)];
end
cobs = zeros(size(periods));
I_R = wavetype=='R'; I_L = wavetype=='L';
if any(I_R); cobs(I_R) = dispR_surf96_TI(periods(I_R),truemod,nmode); end % "observations"
if any(I_L); cobs(I_L) = dispL_surf96_TI(periods(I_L),truemod,nmode); end
cstd = cobs * 0.01; % observation uncertainties

%% Starting model
% load starting model
[startmod, discs] = card2mod_TI(CARD,zmax);
startmod(end,7) = 1; % eta = 1 in the half-space (as for the true model)
% Break the smoothing at the discontinuities of the card (water bottom,
% sediments, Moho). Smoothing across the Moho prevents the inversion from
% fitting the short-period (crust-sensitive) data, and the misfit then
% leaks into the mantle anisotropy. Set discs = [] to smooth across them.

figure(1); clf; set(gcf,'position',[50 50 1200 600]);
plot_model_TI({startmod, truemod}, {'start','true'}, zmax);
subplot(1,4,1); title('Starting and true models');

%% Do linearized inversion using the CPS TI programs
cstart = zeros(size(periods));
if any(I_R); cstart(I_R) = dispR_surf96_TI(periods(I_R),startmod,nmode); end % "predictions"
if any(I_L); cstart(I_L) = dispL_surf96_TI(periods(I_L),startmod,nmode); end

[finalmod, cpre, unc] = run_surf96_inv_RaylLove_Vsvh_Vpvh_Rho_Eta(cobs,cstd,periods,wavetype,startmod,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,eps_xi,eps_eta,z_dampbot,nit,nit_recalc_c,nmode,1,vsv_thresh_xi);

%% Plot final kernels
lgd = {};
Npers = length(vec_T);
clr = jet(Npers);
flds = {'dcdvsv','dcdvsh','dcdvpv','dcdvph','dcdrho','dcdeta'};
xl = {'dc/dVsv','dc/dVsh','dc/dVpv','dc/dVph','dc/d\rho','dc/d\eta'};
for iw = 1:length(data_types)
    wt = data_types(iw);
    kernel = calc_kernel96_TI(finalmod, vec_T, wt, 1, 0, nmode);
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
        ylim([0 zmax]);
        set(gca,'fontsize',13,'linewidth',1.5,'ydir','reverse');
    end
    pos = get(gca,'Position');
    legend(lgd,'location','eastoutside');
    set(gca,'position',pos);
end

%% Plot data fit
figure(103); clf; set(gcf,'position',[300 300 1000 420]);
for iw = 1:length(data_types)
    wt = data_types(iw);
    I_w = wavetype == wt;
    subplot(1,length(data_types),iw); box on; hold on;
    plot(periods(I_w),cstart(I_w),'-ob','linewidth',2);
    plot(periods(I_w),cpre(I_w),'-or','linewidth',2);
    he = errorbar(periods(I_w),cobs(I_w),2*cstd(I_w),'sk'); set(he,'markersize',8,'markerfacecolor','k','linewidth',2);
    legend({'c start','c final','c obs'},'Location','southeast')
    xlabel('Period (s)');
    ylabel([wt,' phase velocity (km/s)']);
    set(gca,'FontSize',15,'linewidth',1.5);
end

%% Plot true and recovered models: Vsv, Vsh, Vpv, Vph, rho, eta, xi
nlayer = size(finalmod,1);
xi = @(m) (m(:,6)./m(:,3)).^2;
% xi uncertainty from the model covariance: xi = Vsh^2/Vsv^2,
%   dxi/dVsv = -2 xi/Vsv, dxi/dVsh = 2 xi/Vsh
xi_final = xi(finalmod);
xi_std = nan(nlayer,1);
for k = 1:nlayer
    if finalmod(k,3) > 0
        Jx = [-2*xi_final(k)/finalmod(k,3), 2*xi_final(k)/finalmod(k,6)];
        ik = [k, nlayer+k];
        xi_std(k) = sqrt(Jx*unc.Cm(ik,ik)*Jx');
    end
end

vals = {@(m) m(:,3), @(m) m(:,6), @(m) m(:,2), @(m) m(:,5), @(m) m(:,4), @(m) m(:,7), xi};
stds = {unc.vsv_std, unc.vsh_std, unc.vpv_std, unc.vph_std, unc.rho_std, unc.eta_std, xi_std};
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
    if k==1; ylabel('Depth (km)'); title(['Data: ',data_types]); end
    ylim([0 z_dampbot+50]);
    set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');
    if k==1; legend([h1 h2 h3],{'true','start','final'},'Location','southwest'); end
end
