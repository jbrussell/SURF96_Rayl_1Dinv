% Benchmark: isotropic vs transversely isotropic (TI) anelastic kernels
%
% Part 1: compares calc_anelastic_kernel96 (isotropic; surf96 + srfker96)
% with calc_anelastic_kernel96_TI for a TI model converted from an
% isotropic model with Qp and Qs (Rayleigh waves only; the isotropic
% version has not been implemented for Love waves). As for the elastic
% kernels (benchmark_TI_kernels.m), the velocity and Q kernels differ by up
% to ~2% for a SPHERICAL EARTH model because srfker96 does not apply the
% earth-flattening chain rule; the tolerance allows for this.
%
% Part 2: self-consistency of calc_anelastic_kernel96_TI (Rayleigh and Love)
%   - its causal phase velocity and attenuation coefficient agree with the
%     values computed by tregn96/tlegn96 (dispR/dispL_surf96_TI and
%     calc_amplification_ellipticity_gamma_TI with Q)
%   - finite-difference checks of dc/dQs^-1, dc/dQp^-1, dgamma/dQs^-1 and
%     dgamma/dQp^-1
%
% jbrussell 11/2026
%
clear; close all;
benchmark_TI_setup;
fprintf('\n===== benchmark_TI_anelastic_kernels =====\n');
fid = fopen(logfile,'a'); fprintf(fid,'\n===== benchmark_TI_anelastic_kernels =====\n'); fclose(fid);

model_Q = benchmark_model_TI(1);         % [H Vp Vs rho Qp Qs]
model_TI_Q = convert_iso2TI(model_Q);    % [H Vpv Vsv rho Vph Vsh eta Qp Qs]
vec_T = [5 10 20 50];
tol_iso = 3e-2;
tol_self = 2e-3;
tol_fd = 2e-2;

%% Part 1: isotropic vs TI (Rayleigh)
[k_iso, d_iso] = calc_anelastic_kernel96(model_Q, vec_T, 'R', 0, 0);
[k_TI, d_TI]   = calc_anelastic_kernel96_TI(model_TI_Q, vec_T, 'R', 0, 0);
flds = {'dcdvs','dcdvp','dudvs','dudvp','dcdQs','dcdQp','dudQs','dudQp','dgdQs','dgdQp','dQdQmu','dQdQkap'};
for ii = 1:length(flds)
    bench_compare_TI(['R anelastic ',flds{ii}], k_TI.(flds{ii}), k_iso.(flds{ii}), tol_iso, 'peak', logfile);
end
% calc_anelastic_kernel96 stores srfker96's dc/dh and dU/dh as dcdrho and dudrho
bench_compare_TI('R anelastic dc/dh  (iso field dcdrho)', k_TI.dcdh, k_iso.dcdrho, tol_iso, 'peak', logfile);
bench_compare_TI('R anelastic dU/dh  (iso field dudrho)', k_TI.dudh, k_iso.dudrho, tol_iso, 'peak', logfile);
bench_compare_TI('R anelastic phase velocity', d_TI.phv, d_iso.phv, tol_self, 'rel', logfile);
bench_compare_TI('R anelastic group velocity', d_TI.grv, d_iso.grv, tol_self, 'rel', logfile);
bench_compare_TI('R attenuation coefficient gamma', d_TI.gamma, d_iso.gamma, tol_iso, 'rel', logfile);

figure; set(gcf,'position',[50 50 1200 600]);
pl = {'dcdvs','dcdvp','dQdQmu','dQdQkap'};
for jj = 1:4
    subplot(1,4,jj); hold on; box on;
    plot(k_iso.(pl{jj}), k_iso.z, '-k', 'linewidth', 2);
    plot(k_TI.(pl{jj}), k_TI.z, '--r', 'linewidth', 1.5);
    axis ij; ylim([0 150]); xlabel(pl{jj}); ylabel('Depth (km)'); title('R: black iso, red TI');
end
print(gcf,'-dpng',[figdir,'anelastic_kernels_iso_vs_TI_R.png']);

%% Part 2: self-consistency (Rayleigh and Love)
z = [0; cumsum(model_Q(1:end-1,1))];
zm = z + model_Q(:,1)/2; zm(end) = z(end);
shape = exp(-((zm-60)/40).^2); shape(1:2) = 0; % exclude water and sediment
for wt = {'R','L'}
    wt = wt{1};
    [k, d] = calc_anelastic_kernel96_TI(model_TI_Q, vec_T, wt, 0, 0);
    a = calc_amplification_ellipticity_gamma_TI(vec_T, model_TI_Q, wt, 0);
    bench_compare_TI([wt,' TI anelastic c vs t',lower(wt),'egn96'], d.phv, a.phv, tol_self, 'rel', logfile);
    bench_compare_TI([wt,' TI gamma vs t',lower(wt),'egn96'], d.gamma, a.gamma, tol_self, 'rel', logfile);

    % finite differences in Q^-1
    for iq = 1:2
        if iq == 1; col = 9; nm = 'Qs'; else; col = 8; nm = 'Qp'; end
        if strcmp(wt,'L') && iq == 2; continue; end
        qinv = 1./model_TI_Q(:,col); qinv(isinf(qinv)) = 0;
        dq = 0.05*qinv.*shape;
        mp = model_TI_Q; mm = model_TI_Q;
        mp(:,col) = 1./(qinv+dq); mm(:,col) = 1./(qinv-dq);
        mp(qinv==0,col) = 0; mm(qinv==0,col) = 0;
        ap = calc_amplification_ellipticity_gamma_TI(vec_T, mp, wt, 0);
        am = calc_amplification_ellipticity_gamma_TI(vec_T, mm, wt, 0);
        bench_compare_TI(sprintf('%s FD dc/d%s^-1', wt, nm), (k.(['dcd',nm])'*dq), (ap.phv(:)-am.phv(:))/2, tol_fd, 'peak', logfile);
        bench_compare_TI(sprintf('%s FD dgamma/d%s^-1', wt, nm), (k.(['dgd',nm])'*dq), (ap.gamma(:)-am.gamma(:))/2, tol_fd, 'peak', logfile);
    end
end
