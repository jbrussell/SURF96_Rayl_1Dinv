% Synthetic seismograms by mode summation for a TI oceanic model:
% all modes vs. the fundamental mode (mode 0) only vs. the first overtone
% (mode 1) only (uses calc_synth96_TI.m: CPS tprep96, tdisp96,
% tlegn96/tregn96, tpulse96, fmech96, f96tosac)
%
% Default: PA5 (Gaherty et al., 1999) and an Mw 6 earthquake 10 km below
% the sea floor, recorded by ocean-bottom seismometers at the sea floor.
% In PA5 the thin crust over the fast lid and the LVZ make the mode-0 Love
% wave at ~10-55 s trapped in the LVZ (small at the sea floor), while the
% crust/lid Love wave is on mode 1 and higher modes (mode osculation; see
% README_TI.md). The comparison shows which modes make up the "fundamental
% mode" Love and Rayleigh waves recorded at the sea floor in each period
% band.
%
% REQUIRES the CPS programs tpulse96 (VOLIII), fmech96 and f96tosac (VOLV),
% which are not in bin_v3.30/: build them from src.zip (PROGRAMS.330)
% and put them in the same directory as the other CPS programs (or set
% CPS_BIN).
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

figdir = './test_synth_figs/'; % output figures ('' = do not save)
if ~isempty(figdir) && ~exist(figdir,'dir'); mkdir(figdir); end

%% Model
CARD = './MINEOS_CARD/pa5_5km.card';
zmax = 600; % [km] model depth (half-space below)
[model, discs] = card2mod_TI(CARD, zmax);
model(end,7) = 1; % eta = 1 in the half-space (see check_halfspace_eta_TI.m)
h_water = sum(model(model(:,3)==0,1)); % water depth (km)

%% Source and receivers
src.strike = 30; src.dip = 60; src.rake = 45; % deg (oblique: excites Love and Rayleigh)
src.Mw = 6.0;
src.az = 60; % source-to-station azimuth (deg)
dist = [1000 2000]; % epicentral distances (km)

par.hs = h_water + 10; % source depth below the sea surface (km): 10 km below the sea floor
par.hr = h_water;      % receiver at the sea floor (moved 0.01 km into the solid)
par.dt = 1;            % s
par.npts = 4096;       % s of record (longer than the slowest arrival to avoid wrap-around)
par.fmax = 0.2;        % Hz: compute periods >= 5 s only (saves time)
par.nmod = 20;         % modes computed (complete for periods >~ 5 s in PA5)
par.wavetype = 'RL';
par.pulse = 't';       % triangle of base 2*L*dt = 6 s (roughly Mw 6)
par.L = 3;
par.output = 'D';      % displacement (cm)
par.modesel = {'all', 0, 1}; % all modes, fundamental mode only, first overtone only
lab = {'all modes','mode 0 (fundamental)','mode 1 (first overtone)'};
clr = {[0 0 0], [0.1 0.3 0.8], [0.85 0.2 0.1]};
lw = [2.5 1.5 1.5];

%% Synthetics
fprintf('Computing synthetics (dispersion for %d modes up to %.2f Hz) ...\n', par.nmod, par.fmax);
tic;
syn = calc_synth96_TI(model, dist, par, src);
toc

%% Plot: broadband (5-150 s) Z, R, T at each distance
comps = {'Z','R','T'};
Tband_bb = [5 150];
vwin = [5.5 2.8]; % [km/s] group-velocity window plotted (and used for the peak amplitudes)
figure(1); clf; set(gcf,'position',[50 50 1400 300*length(dist)+100],'color','w');
for k = 1:length(dist)
    t = syn.t(:,k);
    for ic = 1:3
        subplot(length(dist),3,(k-1)*3+ic); hold on; box on;
        X = squeeze(syn.(comps{ic})(:,k,:));
        Y = bandpass_fft(X, par.dt, Tband_bb(1), Tband_bb(2));
        A = max(abs(Y(:,1)));
        h = [];
        for is = 1:size(Y,2)
            h(is) = plot(t, Y(:,is)/A, '-', 'color', clr{is}, 'linewidth', lw(is));
        end
        xlim(dist(k)./vwin); ylim([-1.5 1.5]);
        title(sprintf('%s, %g km (max %.2g cm)', comps{ic}, dist(k), A));
        if k == length(dist); xlabel('Time (s)'); end
        if ic == 1 && k == 1; legend(h, lab, 'location', 'southwest', 'fontsize', 9); end
        set(gca,'fontsize',12,'linewidth',1.2);
    end
end
if ~isempty(figdir); print(gcf,'-dpng',[figdir,'synth_broadband.png']); end

%% Plot: period bands, transverse (Love) and vertical (Rayleigh), farthest distance
bands = [10 15; 15 30; 30 60; 60 150];
k = length(dist);
t = syn.t(:,k);
for icomp = [3 1]
    figure(10+icomp); clf; set(gcf,'position',[50 50 1300 900],'color','w');
    X = squeeze(syn.(comps{icomp})(:,k,:));
    for ib = 1:size(bands,1)
        subplot(size(bands,1),1,ib); hold on; box on;
        Y = bandpass_fft(X, par.dt, bands(ib,1), bands(ib,2));
        A = max(abs(Y(:,1)));
        h = [];
        for is = 1:size(Y,2)
            h(is) = plot(t, Y(:,is)/A, '-', 'color', clr{is}, 'linewidth', lw(is));
        end
        xlim(dist(k)./vwin); ylim([-1.5 1.5]);
        text(dist(k)/5.5 + 5, 1.25, sprintf('%g-%g s', bands(ib,:)), 'fontsize', 13, 'fontweight', 'bold');
        if ib == 1
            title(sprintf('%s component, %g km, az %g deg (normalized to all modes in each band)', comps{icomp}, dist(k), src.az));
            legend(h, lab, 'location', 'northeast', 'orientation', 'horizontal', 'fontsize', 10);
        end
        if ib == size(bands,1); xlabel('Time (s)'); end
        set(gca,'fontsize',12,'linewidth',1.2);
    end
    if ~isempty(figdir); print(gcf,'-dpng',[figdir,'synth_bands_',comps{icomp},'.png']); end
end

%% Summary: peak amplitude of each mode selection relative to all modes
fprintf('\nPeak amplitude relative to the all-modes synthetic (%g km)\n', dist(k));
fprintf('(within the surface-wave window %.0f-%.0f s)\n', dist(k)./vwin);
iw = t >= dist(k)/vwin(1) & t <= dist(k)/vwin(2);
fprintf('%-10s %-10s %10s %10s\n', 'component', 'band (s)', 'mode 0', 'mode 1');
for icomp = [3 1]
    X = squeeze(syn.(comps{icomp})(:,k,:));
    for ib = 1:size(bands,1)
        Y = bandpass_fft(X, par.dt, bands(ib,1), bands(ib,2));
        a = max(abs(Y(iw,:)), [], 1);
        fprintf('%-10s %4g-%-5g %10.2f %10.2f\n', comps{icomp}, bands(ib,:), a(2)/a(1), a(3)/a(1));
    end
end
