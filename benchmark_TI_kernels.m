% Benchmark: isotropic vs transversely isotropic (TI) sensitivity kernels
%
% Part 1: compares calc_kernel96 (isotropic; surf96 + srfker96) with
% calc_kernel96_TI (TI; tregn96/tlegn96 -DER) for a TI model converted from
% an isotropic model. For this model the TI isotropic-equivalent kernels
%     dc/dVs = dc/dVsv + dc/dVsh,   dc/dVp = dc/dVpv + dc/dVph
% (and the group-velocity equivalents) should match the isotropic kernels.
% The density kernels (calc_kernel96 outputs 6-7) and the srfker96
% thickness kernels dc/dh, dU/dh (outputs 8-9) are also compared.
%
%   Expected differences: ~1e-4 to 1e-3, set by the 4 significant figures
%   of the srfker96 text output and, for Rayleigh waves, by the slightly
%   different earth flattening of the CPS TI programs (see
%   benchmark_TI_dispersion.m). Before November 2026 calc_kernel96 did not
%   apply the spherical-earth mapping to the srfker96 phase-velocity
%   kernels, which then differed from the TI kernels by up to ~2.4% at
%   150 km depth, and its outputs 6-7 were srfker96's dc/dh and dU/dh
%   rather than density kernels.
%
% Part 2: validates every TI kernel (Vpv, Vph, Vsv, Vsh, eta, rho) against
% centered finite differences of dispR_surf96_TI / dispL_surf96_TI for a
% model WITH radial anisotropy (xi = 1.1, phi = 0.95, eta = 0.92 below
% 10 km; eta = 1 in the half-space, see check_halfspace_eta_TI.m). This also validates the group-velocity kernels, which are
% computed from the phase-velocity kernels (Rodi et al., 1975).
%
% jbrussell 11/2026
%
clear; close all;
benchmark_TI_setup;
fprintf('\n===== benchmark_TI_kernels =====\n');
fid = fopen(logfile,'a'); fprintf(fid,'\n===== benchmark_TI_kernels =====\n'); fclose(fid);

model = benchmark_model_TI();      % isotropic [H Vp Vs rho]
model_TI = convert_iso2TI(model);  % TI [H Vpv Vsv rho Vph Vsh eta]
vec_T = [5 10 20 50];
tol_iso = 1e-2;  % peak-normalized tolerance, iso vs TI
tol_vp  = 2e-2;  % same, for the Rayleigh Vp kernels (see below)
% Rayleigh-wave Vp-type kernels are compared with a looser tolerance
% (tol_vp): tdisp96/tregn96 do not apply the earth-flattening transform to
% the modulus F (see README_TI.md), which is what a Vp perturbation
% changes, so for SPHERICAL EARTH models the TI dc/dVp differs from the
% isotropic one by up to ~1.5% at long periods. (The corrected isotropic
% dc/dVp agrees with the sregn96 -DER kernel to < 1e-3.)
tol_fd  = 1e-2;  % kernel prediction vs finite difference, relative to the largest change over periods

%% Part 1: isotropic vs TI kernels
wavetypes = {'R','L'};
for iw = 1:2
    wt = wavetypes{iw};
    [dcdvs, dcdvp, dudvs, dudvp, z_mid, dcdrho_iso, dudrho_iso, dcdh_iso, dudh_iso] = calc_kernel96(model, vec_T, wt, 0, 0);
    [kernel, dispersion] = calc_kernel96_TI(model_TI, vec_T, wt, 0, 0);

    bench_compare_TI([wt,' dc/dVs  (dcdvsv+dcdvsh)'], kernel.dcdvs, dcdvs, tol_iso, 'peak', logfile);
    bench_compare_TI([wt,' dU/dVs  (dudvsv+dudvsh)'], kernel.dudvs, dudvs, tol_iso, 'peak', logfile);
    if strcmp(wt,'R')
        bench_compare_TI([wt,' dc/dVp  (dcdvpv+dcdvph)'], kernel.dcdvp, dcdvp, tol_vp, 'peak', logfile);
        bench_compare_TI([wt,' dU/dVp  (dudvpv+dudvph)'], kernel.dudvp, dudvp, tol_vp, 'peak', logfile);
    end
    bench_compare_TI([wt,' dc/drho'], kernel.dcdrho, dcdrho_iso, tol_iso, 'peak', logfile);
    bench_compare_TI([wt,' dU/drho'], kernel.dudrho, dudrho_iso, tol_iso, 'peak', logfile);
    bench_compare_TI([wt,' dc/dh'], kernel.dcdh, dcdh_iso, tol_iso, 'peak', logfile);
    bench_compare_TI([wt,' dU/dh'], kernel.dudh, dudh_iso, tol_iso, 'peak', logfile);

    figure(iw); clf; set(gcf,'position',[50 50 1200 700]);
    if strcmp(wt,'R')
        pl = {'dcdvs','dcdvp','dudvs','dudvp'}; iso = {dcdvs,dcdvp,dudvs,dudvp};
        lab = {'dc/dVs','dc/dVp','dU/dVs','dU/dVp'};
    else
        pl = {'dcdvs','dudvs','dcdh','dudh'}; iso = {dcdvs,dudvs,dcdh_iso,dudh_iso};
        lab = {'dc/dVs','dU/dVs','dc/dh','dU/dh'};
    end
    for jj = 1:4
        subplot(1,4,jj); hold on; box on;
        plot(iso{jj}, z_mid, '-k', 'linewidth', 2);
        plot(kernel.(pl{jj}), kernel.z, '--r', 'linewidth', 1.5);
        axis ij; ylim([0 150]); xlabel(lab{jj}); ylabel('Depth (km)');
        title([wt,': black iso, red TI']);
    end
    print(gcf,'-dpng',[figdir,'kernels_iso_vs_TI_',wt,'.png']);
end

%% Part 2: finite-difference validation of TI kernels (anisotropic model)
z = [0; cumsum(model(1:end-1,1))];
nlay = size(model,1);
xi = ones(nlay,1); phi = ones(nlay,1); eta = ones(nlay,1);
xi(z>=10) = 1.1; phi(z>=10) = 0.95; eta(z>=10) = 0.92;
eta(end) = 1; % eta ~= 1 in the half-space can give a spurious root in tdisp96 (see check_halfspace_eta_TI.m)
model_ani = convert_iso2TI(model, xi, phi, eta);

% smooth perturbation shape (fraction of the parameter), excluding the water layer
zm = z + model(:,1)/2; zm(end) = z(end);
shape = 0.01*exp(-((zm-40)/40).^2);
shape(1) = 0;

params = {'vpv','vph','vsv','vsh','eta','rho'};
cols   = [  2     5     3     6     7     4 ];
vec_T_fd = [8 15 25]; % (away from the Love cutoff near the half-space SH velocity, where 1% perturbations are strongly nonlinear)
for iw = 1:2
    wt = wavetypes{iw};
    [K] = calc_kernel96_TI(model_ani, vec_T_fd, wt, 0, 0);
    for it = 1:2
        types = {'C','U'};
        typ = types{it};
        for ip = 1:length(params)
            if strcmp(wt,'L') && any(strcmp(params{ip},{'vpv','vph','eta'}))
                continue % no Love sensitivity
            end
            if strcmp(wt,'R') && strcmp(params{ip},'vsh')
                continue % no Rayleigh sensitivity
            end
            dm = model_ani(:,cols(ip)) .* shape;
            mp = model_ani; mp(:,cols(ip)) = mp(:,cols(ip)) + dm;
            mm = model_ani; mm(:,cols(ip)) = mm(:,cols(ip)) - dm;
            if strcmp(wt,'R')
                dv_fd = (dispR_surf96_TI(vec_T_fd,mp,0,typ) - dispR_surf96_TI(vec_T_fd,mm,0,typ))/2;
            else
                dv_fd = (dispL_surf96_TI(vec_T_fd,mp,0,typ) - dispL_surf96_TI(vec_T_fd,mm,0,typ))/2;
            end
            if strcmp(typ,'C'); kname = ['dcd',params{ip}]; else; kname = ['dud',params{ip}]; end
            dv_pred = (K.(kname)' * dm);
            bench_compare_TI(sprintf('FD %s %s: d%s/d%s', wt, typ, typ, params{ip}), dv_pred, dv_fd, tol_fd, 'peak', logfile);
        end
    end
end

% plot the anisotropic-model kernels
figure(3); clf; set(gcf,'position',[50 50 1200 700]);
[KR] = calc_kernel96_TI(model_ani, 25, 'R', 0, 0);
[KL] = calc_kernel96_TI(model_ani, 25, 'L', 0, 0);
pl = {'dcdvsv','dcdvsh','dcdvpv','dcdvph','dcdeta','dcdrho'};
for jj = 1:length(pl)
    subplot(1,6,jj); hold on; box on;
    plot(KR.(pl{jj}), KR.z, '-b', 'linewidth', 1.5);
    plot(KL.(pl{jj}), KL.z, '-r', 'linewidth', 1.5);
    axis ij; ylim([0 150]); xlabel(pl{jj}); ylabel('Depth (km)');
    if jj==1; title('25 s; blue R, red L'); end
end
print(gcf,'-dpng',[figdir,'kernels_TI_anisotropic_model.png']);
