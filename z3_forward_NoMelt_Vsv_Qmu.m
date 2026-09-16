% Load a Vsv(z) and Qmu(z) profile (e.g. from a NoMelt-style inversion),
% build a surf96 layered model, and use this repo's own CPS wrapper
% functions to forward calculate:
%   (1) Rayleigh wave fundamental-mode phase velocity dispersion, c(T)
%   (2) Rayleigh wave attenuation, 1/Q_R(T)
%
% This uses calc_anelastic_kernel96.m, which calls srfker96 (part of
% Herrmann's Computer Programs in Seismology, ./bin_v3.30/) on the actual
% model (including its real Qp/Qs values) and returns the physical
% dispersion (phv, grv) together with the anelastic attenuation
% coefficient gamma, converted to 1/Q_R via
%   dispersion.qinv = 2 .* dispersion.grv .* dispersion.gamma ./ omega;
% (see calc_anelastic_kernel96.m, bottom).
%
% Only Vs and Qmu are provided by the input files, so Vp, density, and
% Qkappa are derived using standard assumptions (see "Derived parameters"
% section below) -- adjust these if better constraints are available.
%
clear
path2BIN = './bin_v3.30/'; % path to surf96 binary
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
    setenv('PATH', [path2BIN,':',PATH]);
end
addpath('./functions/')
% Make binary files executable
!chmod ++x ./bin_v3.30/*

%% Input files: two columns, no header -> [value, depth_km]
vsv_file = './data/NoMelt_rev1/rev1_NoMelt_vsv.csv';
qmu_file = './data/NoMelt_rev1/rev1_NoMelt_Qmu.csv';

vsv_raw = readmatrix(vsv_file);
qmu_raw = readmatrix(qmu_file);

vs_val = vsv_raw(:,1); vs_z = vsv_raw(:,2);
qmu_val = qmu_raw(:,1); qmu_z = qmu_raw(:,2);

% sort by depth
[vs_z, I] = sort(vs_z); vs_val = vs_val(I);
[qmu_z, I] = sort(qmu_z); qmu_val = qmu_val(I);

% Both profiles repeat the same depth value at sharp step transitions
% (e.g. near the LAB) to draw a near-vertical jump. interp1 requires
% strictly increasing grid points, so nudge repeated depths apart by a
% negligible amount (1 mm) to preserve the step while keeping depths unique.
dz_tiny = 1e-6; % km
for i = 2:length(vs_z)
    if vs_z(i) <= vs_z(i-1)
        vs_z(i) = vs_z(i-1) + dz_tiny;
    end
end
for i = 2:length(qmu_z)
    if qmu_z(i) <= qmu_z(i-1)
        qmu_z(i) = qmu_z(i-1) + dz_tiny;
    end
end

%% Plot the input Vsv and Qmu profiles
figure(1); clf;
set(gcf,'position',[100 100 800 700]);

subplot(1,2,1); box on; hold on;
plot(vs_val,vs_z,'-b','linewidth',1.5);
xlabel('Vsv (km/s)');
ylabel('Depth (km)');
title('Vsv profile');
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');

subplot(1,2,2); box on; hold on;
plot(qmu_val,qmu_z,'-r','linewidth',1.5);
xlabel('Q_{\mu}');
title('Q_{\mu} profile');
set(gca,'FontSize',14,'linewidth',1.5,'ydir','reverse');

%% Build a surf96 layered model: [thickness, Vp, Vs, Rho, Qp, Qs]

% Merge the two depth grids and interpolate both profiles onto the union
z_common = unique([vs_z(:); qmu_z(:)]);
z_common = z_common([true; diff(z_common) > 1e-3]); % drop near-duplicate depths

vs = interp1(vs_z, vs_val, z_common, 'linear');
qs = interp1(qmu_z, qmu_val, z_common, 'linear'); % Qs == Qmu

dz = [diff(z_common); 0]; % last layer thickness = 0 -> half-space (CPS/surf96 convention, see e.g. z1_forward_models_vary_sed_thickness.m)

% --- Derived parameters (only Vs and Qmu were provided) ---
vpvs = 1.75; % constant Vp/Vs ratio (consistent with z1_/z2_ example scripts in this repo)
Qkappa = 9999; % ~elastic bulk modulus (standard assumption for mantle attenuation)

vp = vpvs * vs;

% Nafe-Drake density from Vp (Brocher, 2005, BSSA 95, eq. 1), Vp in km/s -> rho in g/cc
rho = 1.6612*vp - 0.4721*vp.^2 + 0.0671*vp.^3 - 0.0043*vp.^4 + 0.000106*vp.^5;

% Qp from Qmu (=Qs) and Qkappa via standard partitioning
% (e.g. Dahlen & Tromp, 1998, eq. 9.55): 1/Qp = L/Qmu + (1-L)/Qkappa,
% L = (4/3)*(Vs/Vp)^2
L = (4/3) * (1/vpvs)^2;
qp = 1 ./ ( L./qs + (1-L)/Qkappa );

startmod = [dz(:), vp(:), vs(:), rho(:), qp(:), qs(:)];
nlayer = size(startmod,1);
fprintf('Built %d-layer model (+ half-space), depth 0-%.1f km\n', nlayer, z_common(end));
fprintf('Assumptions: Vp = %.2f*Vs, density from Brocher (2005) Nafe-Drake fit to Vp,\n', vpvs);
fprintf('  Qkappa = %.0f (~elastic), Qp derived from Qmu & Qkappa (L=%.4f)\n', Qkappa, L);

%% Calculate Rayleigh wave dispersion + attenuation using CPS (srfker96)
vec_T = logspace(log10(3),log10(300),45); % periods [s]
Nmode = 0; % fundamental mode

ifnorm = 0; % kernels unnormalized by layer thickness
ifplot = 0;
[kernel, dispersion] = calc_anelastic_kernel96(startmod, vec_T, 'R', ifnorm, ifplot, Nmode);

phv_R = dispersion.phv;   % Rayleigh phase velocity [km/s], physically dispersed by Q
grv_R = dispersion.grv;   % Rayleigh group velocity [km/s]
qinv_R = dispersion.qinv; % 1/Q_R
Q_R = 1 ./ qinv_R;
T = dispersion.period;

%% Save results
out = [T(:), phv_R(:), grv_R(:), qinv_R(:), Q_R(:)];
writematrix(out, './data/NoMelt_rev1/rayleigh_dispersion_attenuation.csv');
fprintf('Saved ./data/NoMelt_rev1/rayleigh_dispersion_attenuation.csv\n');
fprintf('  columns: period_s, phase_vel_km_s, group_vel_km_s, Qinv_R, Q_R\n');

%% Plot dispersion and attenuation
figure(2); clf;
set(gcf,'position',[150 100 700 800]);

subplot(2,1,1); box on; hold on;
plot(T,phv_R,'-o','color',[0.5 0 0.5],'linewidth',1.5,'markersize',4);
ylabel('Phase velocity (km/s)');
title('Rayleigh wave fundamental-mode phase velocity dispersion');
set(gca,'FontSize',14,'linewidth',1.5,'xscale','log');
grid on;

subplot(2,1,2); box on; hold on;
plot(T,qinv_R,'-o','color',[0 0.5 0],'linewidth',1.5,'markersize',4);
xlabel('Period (s)');
ylabel('1/Q_R');
title('Rayleigh wave attenuation');
set(gca,'FontSize',14,'linewidth',1.5,'xscale','log');
grid on;
