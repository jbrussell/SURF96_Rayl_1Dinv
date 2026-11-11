% Benchmark: isotropic vs transversely isotropic (TI) anelastic kernels
%
% Part 1: compares calc_anelastic_kernel96 (isotropic; surf96 + srfker96)
% with calc_anelastic_kernel96_TI for a TI model converted from an
% isotropic model with Qp and Qs, for Rayleigh and Love waves.
%
% Part 2: self-consistency of calc_anelastic_kernel96_TI (Rayleigh and Love)
%   - its causal phase velocity and attenuation coefficient agree with the
%     values computed by tregn96/tlegn96 when run directly with the Q model
%     (the TI wrappers run CPS elastically and apply Q afterwards with the
%     tregn96 formulas, see calc_der96_TI.m)
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
tol_iso = 1e-2;
tol_vp  = 2e-2; % Rayleigh Vp / Qp kernels
% Rayleigh-wave Vp-type kernels are compared with a looser tolerance
% (tol_vp): tdisp96/tregn96 do not apply the earth-flattening transform to
% the modulus F (see README_TI.md), which is what a Vp perturbation
% changes, so for SPHERICAL EARTH models the TI dc/dVp differs from the
% isotropic one by up to ~1.5% at long periods. (The corrected isotropic
% dc/dVp agrees with the sregn96 -DER kernel to < 1e-3.)
tol_self = 2e-3;
tol_fd = 2e-2;

%% Part 1: isotropic vs TI (Rayleigh and Love)
for wt = {'R','L'}
    wt = wt{1};
    [k_iso, d_iso] = calc_anelastic_kernel96(model_Q, vec_T, wt, 0, 0);
    [k_TI, d_TI]   = calc_anelastic_kernel96_TI(model_TI_Q, vec_T, wt, 0, 0);
    if strcmp(wt,'R')
        flds = {'dcdvs','dcdvp','dudvs','dudvp','dcdrho','dudrho','dcdh','dudh', ...
            'dcdQs','dcdQp','dudQs','dudQp','dgdQs','dgdQp','dQdQmu','dQdQkap'};
    else
        flds = {'dcdvs','dudvs','dcdrho','dudrho','dcdh','dudh','dcdQs','dudQs','dgdQs','dQdQmu'};
    end
    for ii = 1:length(flds)
        tol = tol_iso;
        if any(strcmp(flds{ii},{'dcdvp','dudvp','dcdQp','dudQp','dgdQp','dQdQkap'})); tol = tol_vp; end
        bench_compare_TI([wt,' anelastic ',flds{ii}], k_TI.(flds{ii}), k_iso.(flds{ii}), tol, 'peak', logfile);
    end
    bench_compare_TI([wt,' anelastic phase velocity'], d_TI.phv, d_iso.phv, tol_self, 'rel', logfile);
    bench_compare_TI([wt,' anelastic group velocity'], d_TI.grv, d_iso.grv, tol_self, 'rel', logfile);
    bench_compare_TI([wt,' attenuation coefficient gamma'], d_TI.gamma, d_iso.gamma, tol_iso, 'rel', logfile);
    if strcmp(wt,'R'); k_iso_R = k_iso; k_TI_R = k_TI; end
end
k_iso = k_iso_R; k_TI = k_TI_R;

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
    % run tprep96/tdisp96/t?egn96 directly with the Q model, one period at
    % a time, and take the fundamental mode
    c_cps = nan(size(vec_T)); g_cps = nan(size(vec_T));
    writemod_surf96_TI(model_TI_Q, 'startQ.mod', 1);
    for ip = 1:length(vec_T)
        fid = fopen('periodsQ.txt','w'); fprintf(fid,'%.8f\n',vec_T(ip)); fclose(fid);
        system(['tprep96 -M startQ.mod -NMOD 1 -',wt,' -PARR periodsQ.txt > /dev/null; tdisp96 > /dev/null; t',lower(wt),'egn96 -DER > /dev/null']);
        e = read_der96_TI(['t',lower(wt),'egn96.der'], wt, 0);
        c_cps(ip) = e.phv(1); g_cps(ip) = e.gamma(1);
    end
    cleanup_cps96_TI(); delete('startQ.mod','periodsQ.txt');
    bench_compare_TI([wt,' TI anelastic c vs t',lower(wt),'egn96'], d.phv, c_cps, tol_self, 'rel', logfile);
    bench_compare_TI([wt,' TI gamma vs t',lower(wt),'egn96'], d.gamma, g_cps, tol_self, 'rel', logfile);
    a = calc_amplification_ellipticity_gamma_TI(vec_T, model_TI_Q, wt, 0);
    bench_compare_TI([wt,' TI amplification c vs t',lower(wt),'egn96'], a.phv, c_cps, tol_self, 'rel', logfile);

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
