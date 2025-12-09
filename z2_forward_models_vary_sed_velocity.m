% Test inversion using surf96 kernels
%
% This version solves for Vs, Vp, and Rho by enforcing Vp/Vs and Rho/Vs of
% the starting model. The inversion includes kernels for Vs, Vp, Rho such
% that each iteration accounts for perturbations in all three.
%
% jbrussell 6/5/2020
% modified 11/23/2021
% modified 11/20/2022
%
%
%  WARNING!!!! writemod_surf96 set to FLAT EARTH
%
clear
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

nf = 100; %200; % Number of frequencies
fmin = 0.2; %0.1; % [Hz] minimum frequency
fmax = 10; % [Hz] maximum frequency
nmr = 10; %20; % number of Rayleigh modes
nml = 10; %20; % number of Love modes
nks = 500; %1000; % Number of k values for numeric integrals

%% Parameters for phv, grv, H/V
vec_T = [5:0.05:25]; %[3:0.05:25]; %[1:0.05:15];
Nmode = 0; % mode to plot

%% Starting model

% Sediment properties
zsed = 0.5; %0.25; % depth of base of sediment (i.e., sediment thickness)
v_vs_sed = [0.5 1.0 1.5 2.0 2.5]; %
rho_sed = 1.8;

% Crust properties
zcrust = 35; % Depth of base of crust
vs_crust = 3.5;
rho_crust = 2.7;
vpvs = 1.75;

% Mantle properties
vs_mantle = 4.5;
rho_mantle = 3.4;

ylims = [0 5];
xlims = [0 4.0];

%%
N = 1;
clrs_mod = flipud(brewermap(length(v_vs_sed)+N,'Greys'));
clrs_mod = clrs_mod(1:end-N,:);
clrs_hvsr = flipud(brewermap(length(v_vs_sed)+N,'Greens'));
clrs_hvsr = clrs_hvsr(1:end-N,:);
clrs_ellip = flipud(brewermap(length(v_vs_sed)+N,'Blues'));
clrs_ellip = clrs_ellip(1:end-N,:);
clrs_phv = flipud(brewermap(length(v_vs_sed)+N,'Reds'));
clrs_phv = clrs_phv(1:end-N,:);
clrs_grv = flipud(brewermap(length(v_vs_sed)+N,'Purples'));
clrs_grv = clrs_grv(1:end-N,:);

lw = 4;
for ii = 1:length(v_vs_sed)
    
    vs_sed = v_vs_sed(ii);
    disp(['vs_sed = ',num2str(vs_sed),' km/s'])
    
    % Make starting model
    dzz = 10; % km
    zmax = 100 + zcrust;
    z = [0 zsed [zcrust : dzz : zmax]];
    dz = [diff(z) 0];
    vs = [vs_sed vs_crust vs_mantle*ones(size([zcrust : dzz : zmax]))];
    % vp = [3.0 5.8*ones(size([zsed : dzz : 100]))]; %1.75*vs;
    vp = vpvs * vs;
    rho = [rho_sed rho_crust rho_mantle*ones(size([zcrust : dzz : zmax])) ];
    startmod = [dz(:), vp(:), vs(:), rho(:)];
    startmod(1:end-1,:);
    % discs = [zmoho]; % [km] depth to sharp discontinuities
    discs = []; % for this simple example, don't allow discontinuities
    
    % Calculate phase velocity
    disp('Calculating phase velocity');
    c_R = dispR_surf96(vec_T,startmod,Nmode); % "predictions";
    
    % Calculate group velocity
    disp('Calculating group velocity');
    grv_R = dispR_surf96(vec_T,startmod,Nmode,'U'); % "predictions";
        
    % Calculate Rayleigh wave ellipticity
    disp('Calculating ellipticity');
%     disper = calc_amplification_ellipticity_gamma(vec_T,startmod,'R',Nmode);
    % Fund mode Rayleigh ellipticity is faster using HVf!
    nmr_fund_only = 1;
    nml_none = 0;
    [freq_ellip,hvsr_fundrayl] = calc_HVf_SW(startmod, nf,1./max(vec_T),1./min(vec_T),nmr_fund_only,nml_none);
    vec_T_ellip = 1./freq_ellip;
    
    % Calculate HVSR
    disp('Calculating HVSR');
%     fmin = 1./max(vec_T);
%     fmax = 1./min(vec_T);
    [freq_hvsr,hvsr_SW_BW] = calc_HVf_SW_BW(startmod, nf,fmin,fmax,nmr,nml,nks);
    vec_T_hvsr = 1./freq_hvsr;
    
    if isempty(grv_R)
        continue
    end

    %%
    if ii == 1
        figure(1); clf;
        set(gcf,'position',[70          81        1526         842]);
    end

    subplot(1,3,1);
    box on; hold on;
    h = plotlayermods(startmod(:,1),startmod(:,3),'-');
    h.LineWidth = lw;
    h.Color = clrs_mod(ii,:);
    xlabel('Vs (km/s)');
    ylabel('Depth (km)');
    title('Varying Sed. Velocity');
    set(gca,'FontSize',18,'linewidth',1.5);
    ylim(ylims);
    xlim(xlims);
    
    % HVSR
    subplot(2,3,2);
    box on; hold on;
    plot(vec_T_hvsr,hvsr_SW_BW,'-','color',clrs_hvsr(ii,:),'linewidth',lw);
%     plot(vec_T(1:length(disper.RZ)),abs(disper.RZ),'-','color',clrs_ellip(ii,:),'linewidth',lw);
    plot(vec_T_ellip,hvsr_fundrayl,'-','color',clrs_hvsr(ii,:),'linewidth',lw);
    xlabel('Period (s)');
    ylabel('HVSR');
    set(gca,'FontSize',18,'linewidth',1.5);
    xlim([min(vec_T_hvsr) max(vec_T_hvsr)]);
    title('Spectral Ratio (0.1-5 s)','Color',clrs_hvsr(1,:));
    
    % Rayleigh wave ellipticity H/V
    subplot(2,3,3);
    box on; hold on;
%     plot(vec_T(1:length(disper.RZ)),abs(disper.RZ),'-','color',clrs_ellip(ii,:),'linewidth',lw);
    plot(vec_T_ellip,hvsr_fundrayl,'-','color',clrs_ellip(ii,:),'linewidth',lw);
    plot(vec_T_hvsr,hvsr_SW_BW,'-','color',clrs_hvsr(ii,:),'linewidth',lw);
    xlabel('Period (s)');
    ylabel('H/V');
    set(gca,'FontSize',18,'linewidth',1.5);
    xlim([min(vec_T) max(vec_T)]);
    title('Rayleigh Ellipticity (5-25 s)','Color',clrs_ellip(1,:));
    
    % Phase velocity
    subplot(2,3,5);
    box on; hold on;
    plot(vec_T(1:length(c_R)),c_R,'-','color',clrs_phv(ii,:),'linewidth',lw);
%     plot(vec_T(1:length(grv_R)),grv_R,'--','color',clrs(ii,:),'linewidth',2);
    xlabel('Period (s)');
    ylabel('Velocity (km/s)');
    set(gca,'FontSize',18,'linewidth',1.5);
    xlim([min(vec_T) max(vec_T)]);
    title('Phase Velocity (5-25 s)','Color',clrs_phv(1,:));
    
    % Group velocity
    subplot(2,3,6);
    box on; hold on;
%     plot(vec_T(1:length(c_R)),c_R,'-','color',clrs(ii,:),'linewidth',2);
    plot(vec_T(1:length(grv_R)),grv_R,'-','color',clrs_grv(ii,:),'linewidth',lw);
    xlabel('Period (s)');
    ylabel('Velocity (km/s)');
    set(gca,'FontSize',18,'linewidth',1.5);
    xlim([min(vec_T) max(vec_T)]);
    title('Group Velocity (5-25 s)','Color',clrs_grv(1,:));
    
    drawnow;
end

% save2pdf('./figs/z2_vary_sed_velocity.pdf',1,300);


