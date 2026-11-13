% Benchmark: isotropic vs transversely isotropic (TI) dispersion
%
% Compares dispR_surf96 / dispL_surf96 (isotropic, surf96) with
% dispR_surf96_TI / dispL_surf96_TI (TI, tprep96/tdisp96/tregn96/tlegn96)
% for a TI model converted from an isotropic model (convert_iso2TI.m with
% xi = phi = eta = 1). The two should agree.
%
% Expected differences: with SPHERICAL EARTH models (as written by
% writemod_surf96 and writemod_surf96_TI) the TI programs use a slightly
% different Earth flattening than the isotropic programs (Earth radius 6371
% vs 6370 km, and for Rayleigh waves the modulus F is not flattened in
% tdisp96/tregn96). This gives relative differences of ~1e-4 (Rayleigh) and
% ~1e-5 (Love). For FLAT EARTH models the two agree to the output precision.
%
% Also shown (not a pass/fail test): Love and Rayleigh dispersion for a
% model with radial anisotropy xi = (Vsh/Vsv)^2 = 1.1 in the mantle.
%
% jbrussell 11/2026
%
clear; % (no 'close all', so figures from the other benchmark scripts are kept; this script uses figures 101-105)
benchmark_TI_setup;
fprintf('\n===== benchmark_TI_dispersion =====\n');
fid = fopen(logfile,'a'); fprintf(fid,'\n===== benchmark_TI_dispersion =====\n'); fclose(fid);

tol = 1e-3; % relative tolerance

model = benchmark_model_TI();      % isotropic [H Vp Vs rho]
model_TI = convert_iso2TI(model);  % TI [H Vpv Vsv rho Vph Vsh eta]

vec_T = [3:1:10 12:2:30 35:5:100];
modes = [0 1];
wavetypes = {'R','L'};
types = {'C','U'};

% plot the models
figure(101); clf; set(gcf,'position',[50 50 1200 600]);
plot_model_TI({model, benchmark_model_TI(1)}, {'benchmark model','with Q (Q tests)'}, 200);
subplot(1,5,1); title('benchmark_TI_dispersion: models','interpreter','none');
print(gcf,'-dpng',[figdir,'dispersion_models.png']);

for iw = 1:2
    figure(101+iw); clf;
    set(gcf,'position',[100 100 900 700]);
    for it = 1:2
        for im = 1:length(modes)
            nmode = modes(im);
            if strcmp(wavetypes{iw},'R')
                v_iso = dispR_surf96(vec_T,model,nmode,types{it});
                v_TI  = dispR_surf96_TI(vec_T,model_TI,nmode,types{it});
            else
                v_iso = dispL_surf96(vec_T,model,nmode,types{it});
                v_TI  = dispL_surf96_TI(vec_T,model_TI,nmode,types{it});
            end
            if isempty(v_iso); v_iso = nan(size(v_TI)); end
            bench_compare_TI(sprintf('elastic %s %s mode %d',wavetypes{iw},types{it},nmode), ...
                v_TI, v_iso, tol, 'rel', logfile);

            subplot(2,2,it); hold on; box on;
            plot(vec_T,v_iso,'-k','linewidth',2);
            plot(vec_T,v_TI,'--r','linewidth',1.5);
            xlabel('Period (s)'); ylabel([types{it},' (km/s)']);
            title(sprintf('%s %s (black: iso, red: TI)',wavetypes{iw},types{it}));
            subplot(2,2,it+2); hold on; box on;
            plot(vec_T,(v_TI-v_iso)./v_iso*100,'-o');
            xlabel('Period (s)'); ylabel('(TI - iso)/iso (%)');
            title('Relative difference, modes 0,1');
        end
    end
    print(gcf,'-dpng',[figdir,'dispersion_elastic_',wavetypes{iw},'.png']);
end

%% With Q (causal dispersion correction): Q-corrected phase velocity
model_Q = benchmark_model_TI(1);         % [H Vp Vs rho Qp Qs]
model_TI_Q = convert_iso2TI(model_Q);    % [H Vpv Vsv rho Vph Vsh eta Qp Qs]
for iw = 1:2
    if strcmp(wavetypes{iw},'R')
        c_iso = dispR_surf96(vec_T,model_Q,0,'C');
        c_TI  = dispR_surf96_TI(vec_T,model_TI_Q,0,'C');
    else
        % dispL_surf96 has no Q option; use the isotropic eigenfunction code
        d = calc_amplification_ellipticity_gamma(vec_T,model_Q,'L',0);
        c_iso = d.phv(:);
        c_TI  = dispL_surf96_TI(vec_T,model_TI_Q,0,'C');
    end
    bench_compare_TI(sprintf('anelastic %s C mode 0',wavetypes{iw}), c_TI, c_iso, tol, 'rel', logfile);
end

%% Demonstration: radial anisotropy in the mantle
xi = ones(size(model,1),1);
z = [0; cumsum(model(1:end-1,1))];
xi(z>=10) = 1.1;
model_TI_xi = convert_iso2TI(model,xi);
cR0 = dispR_surf96_TI(vec_T,model_TI); cR1 = dispR_surf96_TI(vec_T,model_TI_xi);
cL0 = dispL_surf96_TI(vec_T,model_TI); cL1 = dispL_surf96_TI(vec_T,model_TI_xi);
figure(105); clf; set(gcf,'position',[50 50 1200 600]);
plot_model_TI({model_TI, model_TI_xi}, {'isotropic','\xi = 1.1 below 10 km'}, 200);
subplot(1,4,1); title('Radial anisotropy demonstration: models');
print(gcf,'-dpng',[figdir,'dispersion_radial_anisotropy_models.png']);

figure(104); clf; hold on; box on;
plot(vec_T,cR0,'-b',vec_T,cR1,'--b',vec_T,cL0,'-r',vec_T,cL1,'--r','linewidth',1.5);
legend({'Rayleigh, isotropic','Rayleigh, \xi=1.1 mantle','Love, isotropic','Love, \xi=1.1 mantle'},'location','southeast');
xlabel('Period (s)'); ylabel('Phase velocity (km/s)');
title('Effect of radial anisotropy (Voigt Vs fixed)');
print(gcf,'-dpng',[figdir,'dispersion_radial_anisotropy.png']);
