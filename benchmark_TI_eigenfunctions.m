% Benchmark: isotropic vs transversely isotropic (TI) eigenfunctions,
% amplification/ellipticity/attenuation, and seafloor admittance
%
% Compares, for a TI model converted from an isotropic model:
%   calc_eigenfunctions96                vs calc_eigenfunctions96_TI
%   calc_amplification_ellipticity_gamma vs calc_amplification_ellipticity_gamma_TI
%   calc_admittance96_Ruan14             vs calc_admittance96_Ruan14_TI
%
% Notes: the overall amplitude of the eigenfunctions is arbitrary in CPS and
% is occasionally chosen differently by the isotropic and TI programs; the
% eigenfunctions are therefore compared after rescaling (the amplification
% factors A_R and A_L, which are compared separately, absorb the scaling).
% The isotropic eigenfunctions are read from the sdpder96 text file,
% which has 3 significant figures, and the isotropic dispersion values from
% the sdpegn96 table (4-5 digits); the tolerances reflect this. Small
% additional differences come from the different spherical-earth
% flattening used by the CPS TI programs (see benchmark_TI_dispersion.m).
%
% jbrussell 11/2026
%
clear; close all;
benchmark_TI_setup;
fprintf('\n===== benchmark_TI_eigenfunctions =====\n');
fid = fopen(logfile,'a'); fprintf(fid,'\n===== benchmark_TI_eigenfunctions =====\n'); fclose(fid);

model = benchmark_model_TI();      % isotropic [H Vp Vs rho]
model_TI = convert_iso2TI(model);  % TI [H Vpv Vsv rho Vph Vsh eta]
vec_T = [4 6 8 10 15 20 30 50];
tol_eig = 2e-2;  % peak-normalized
tol_disp = 2e-3; % relative

%% Eigenfunctions
for wt = {'R','L'}
    wt = wt{1};
    for nmode = 0:1
        eig_iso = calc_eigenfunctions96(vec_T, model, wt, nmode);
        eig_TI  = calc_eigenfunctions96_TI(vec_T, model_TI, wt, nmode);
        % the isotropic reader only returns periods at which the mode exists
        [~, I] = ismember(round(eig_iso.periods*1e3), round(eig_TI.periods*1e3));
        if strcmp(wt,'R'); flds = {'ur','uz','tr','tz'}; else; flds = {'ut','tt'}; end
        % The overall amplitude (and sign) of a CPS eigenfunction is arbitrary
        % and is occasionally chosen differently by slegn96/sregn96 and
        % tlegn96/tregn96 (the amplification factor A then differs by the
        % square of the same factor, so the excitation is unchanged). Rescale
        % each TI period by a single factor, from the vertical (R) or
        % transverse (L) displacement, before comparing shapes.
        if strcmp(wt,'R'); fref_eig = 'uz'; else; fref_eig = 'ut'; end
        nscaled = 0;
        for ip = 1:length(I)
            b = eig_iso.(fref_eig)(:,ip); a = eig_TI.(fref_eig)(:,I(ip));
            [~, imax] = max(abs(b));
            scl = b(imax)/a(imax);
            if abs(scl-1) > 1e-2; nscaled = nscaled + 1; end
            for ii = 1:length(flds)
                eig_TI.(flds{ii})(:,I(ip)) = eig_TI.(flds{ii})(:,I(ip))*scl;
            end
        end
        if nscaled > 0
            msg = sprintf('       (%s mode %d: TI eigenfunction amplitude rescaled at %d of %d periods)', wt, nmode, nscaled, length(I));
            fprintf('%s\n', msg); fid = fopen(logfile,'a'); fprintf(fid,'%s\n',msg); fclose(fid);
        end
        for ii = 1:length(flds)
            bench_compare_TI(sprintf('%s mode %d eigenfunction %s',wt,nmode,flds{ii}), ...
                eig_TI.(flds{ii})(:,I), eig_iso.(flds{ii}), tol_eig, 'peak', logfile);
        end
        if nmode == 0
            figure; set(gcf,'position',[50 50 1100 600]);
            for ii = 1:length(flds)
                subplot(1,length(flds),ii); hold on; box on;
                plot(eig_iso.(flds{ii}), eig_iso.z, '-k', 'linewidth', 2);
                plot(eig_TI.(flds{ii}), eig_TI.z, '--r', 'linewidth', 1.5);
                axis ij; ylim([0 150]); xlabel(flds{ii}); ylabel('Depth (km)');
                title([wt,' mode 0: black iso, red TI']);
            end
            print(gcf,'-dpng',[figdir,'eigenfunctions_',wt,'.png']);
        end
    end
end

%% Amplification, ellipticity, attenuation (with and without Q)
model_Q = benchmark_model_TI(1);
model_TI_Q = convert_iso2TI(model_Q);
for iq = 1:2
    if iq == 1
        m_iso = model; m_TI = model_TI; lab = 'default Q';
    else
        m_iso = model_Q; m_TI = model_TI_Q; lab = 'model Q';
    end
    for wt = {'R','L'}
        wt = wt{1};
        d_iso = calc_amplification_ellipticity_gamma(vec_T, m_iso, wt, 0);
        d_TI  = calc_amplification_ellipticity_gamma_TI(vec_T, m_TI, wt, 0);
        [~, I] = ismember(round(d_iso.period*1e3), round(d_TI.period*1e3));
        if strcmp(wt,'R'); flds = {'phv','grv','gamma','A_R','RZ'}; else; flds = {'phv','grv','gamma','A_L'}; end
        for ii = 1:length(flds)
            tol = tol_disp;
            if any(strcmp(flds{ii},{'gamma','A_R','A_L','RZ'}))
                tol = 1e-2; % 5 significant figures in the isotropic table
            end
            bench_compare_TI(sprintf('%s %s (%s)',wt,flds{ii},lab), d_TI.(flds{ii})(I), d_iso.(flds{ii}), tol, 'rel', logfile);
        end
    end
end

%% Seafloor admittance (Ruan et al., 2014)
adm_iso = calc_admittance96_Ruan14(vec_T, model, 0);
adm_TI  = calc_admittance96_Ruan14_TI(vec_T, model_TI, 0);
[~, I] = ismember(round(adm_iso.periods*1e3), round(adm_TI.periods*1e3));
bench_compare_TI('R admittance (Ruan et al., 2014)', adm_TI.admittance(I), adm_iso.admittance, tol_eig, 'rel', logfile);

figure; hold on; box on;
plot(adm_iso.periods, adm_iso.admittance, '-ok', 'linewidth', 2);
plot(adm_TI.periods, adm_TI.admittance, '--xr', 'linewidth', 1.5);
xlabel('Period (s)'); ylabel('Admittance (m/Pa)'); title('Black: iso, red: TI');
print(gcf,'-dpng',[figdir,'admittance.png']);
