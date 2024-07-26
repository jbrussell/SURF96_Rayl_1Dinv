clear; close all;

% Setup paths to surf96 codes
path2BIN = './bin_v3.30/'; % path to surf96 binary
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
%     setenv('PATH', [PATH,':',path2BIN]);
    setenv('PATH', [path2BIN,':',PATH]);
end
addpath('./functions/')
% Make binary files executable
!chmod ++x ./bin_v3.30/*

% Setup paths to HVSR codes
path2BIN = './bin_HV_DFA/'; % path to HV binary
PATH = getenv('PATH');
if isempty(strfind(PATH,path2BIN))
%     setenv('PATH', [PATH,':',path2BIN]);
    setenv('PATH', [path2BIN,':',PATH]);
end
!chmod ++x ./bin_HV_DFA/*

%% Parameters for HVf

nf = 200; % Number of frequencies
fmin = 0.1; % [Hz] minimum frequency
fmax = 10; % [Hz] maximum frequency
nmr = 20; % number of Rayleigh modes
nml = 20; % number of Love modes
nks = 1000; % Number of k values for numeric integrals

%% Starting model

zmax = 1; % km

ylims = [0 zmax];

% Make homogeneous starting model
dzz = 50/1000; % km
z = [0 : dzz : zmax];
dz = [diff(z) 0];
vs = linspace(0.5,0.5,length(z));
vs(9:end) = vs(9:end)*2; % insert discontinuity
vp = 1.75*vs;
rho = 3*vs;
startmod = [dz(:), vp(:), vs(:), rho(:)];
startmod(1:end-1,:);
% discs = [zmoho]; % [km] depth to sharp discontinuities
discs = []; % for this simple example, don't allow discontinuities

figure(1); clf;
set(gcf,'position',[ 330         339        1150         686],'color','w');

subplot(2,3,[1 4]);
box on; hold on;
h = plotlayermods(startmod(:,1),startmod(:,3),'-b');
h.LineWidth = 2;
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title('Starting Model');
set(gca,'FontSize',18,'linewidth',1.5);
ylim(ylims);

%% Calculate HVSR

[~,hvsr_SW] = calc_HVf_SW(startmod, nf,fmin,fmax,nmr,nml);
[freq,hvsr_SW_BW] = calc_HVf_SW_BW(startmod, nf,fmin,fmax,nmr,nml,nks);
[~,hvsr_BW] = calc_HVf_SW_BW(startmod, nf,fmin,fmax,0,0,nks);
[~,hvsr_Rayl] = calc_HVf_SW(startmod, nf,fmin,fmax,nmr,0);

subplot(2,3,[2 3]);
box on; hold on;
plot(freq,hvsr_SW,'-r','linewidth',2,'DisplayName','HVSR (SW only)');
plot(freq,hvsr_Rayl,'--m','linewidth',2,'DisplayName','HVSR (Rayl only)');
plot(freq,hvsr_BW,'-k','linewidth',2,'DisplayName','HVSR (BW only)');
plot(freq,hvsr_SW_BW,'-b','linewidth',2,'DisplayName','HVSR (SW + BW)');
xlabel('Freq. (Hz)');
ylabel('H/V');
% title('Starting Model');
set(gca,'FontSize',18,'linewidth',1.5,'xscale','log');
pos = get(gca,'Position');
legend('location','northeast');
set(gca,'position',pos);

% [~,log] = system(['HVf -h']);

% HVf -f model.txt -nf 100 -fmin 0.1 -fmax 10 -logsam -nmr 20 -nml 20 -ph -gr -hv  > HV.dat


% Calculate Ellipticity
vec_T = flip(1 ./ freq);
nmode = 0; % fund mode
disper = calc_amplification_ellipticity_gamma(vec_T,startmod,'R',nmode);

% Calculate HVSR for fundamental mode Rayleigh only
nmr_fund_only = 1;
nml_none = 0;
[~,hvsr_fundrayl] = calc_HVf_SW(startmod, nf,fmin,fmax,nmr_fund_only,nml_none);

subplot(2,3,[5 6]);
box on; hold on;
plot(1 ./ disper.period,disper.RZ,'-m','linewidth',2,'DisplayName','Rayl. Ellipticity');
plot(1 ./ disper.period,abs(disper.RZ),'-r','linewidth',2,'DisplayName','|Rayl. Ellipticity|');
plot(freq,hvsr_fundrayl,'--b','linewidth',2,'DisplayName','HVSR (Fund Rayl only)');
xlabel('Freq. (Hz)');
ylabel('H/V');
% title('Starting Model');
set(gca,'FontSize',18,'linewidth',1.5,'xscale','log');
pos = get(gca,'Position');
legend('location','northeast');
set(gca,'position',pos);
