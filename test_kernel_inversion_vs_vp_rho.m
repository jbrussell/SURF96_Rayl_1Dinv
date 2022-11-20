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

%% Inversion regularization parameters
eps_data = 1; % data fit
eps_H = 0.05; % norm damping
eps_J = 0.001; % first derivative smoothing
eps_F = 2; % second derivative smoothing
eps_vpvs = 10; % enforce vp/vs ratio
eps_rhovs = 10; % enforce rho/vs ratio

z_dampbot = 150; % [km] damp to starting model below this depth
zmax = 400; % [km] maximum depth of model space

nmode = 0; % mode branch number (0=fund, 1=1st overtone, 2=2nd overtone, etc.)

% Other inversion parameters
nit = 8; % total number of iterations
nit_recalc_c = 1; % number of iterations after which to recalculate phase velocity and kernels

%% Generate the synthetic dataset for this test
%Read in MINEOS model and convert to SURF96 layered model format
cardname = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_Line1_dist100.00km';

% Get MINEOS model
CARD = ['./MINEOS_CARD/',cardname,'.card'];
[truemod, discs] = card2mod(CARD,zmax); % true model
% % Keep track of discontinuities
% waterdepth = discs(1);
% seddepth = discs(2);
% mohodepth = discs(end);

% perturb true model
z = [0; cumsum(truemod(:,1))];
Ibot=find(z<=80); Ibot=Ibot(end);
truemod(2:Ibot,2) = truemod(2:Ibot,2)*1.1;
truemod(2:Ibot,3) = truemod(2:Ibot,3)*1.1;
truemod(2:Ibot,4) = truemod(2:Ibot,4)*1.1;

% GENERATE SYNTHETIC DATASET
% Calculate dispersion for perturbed starting model and true model (which is like the data)
vec_T = logspace(log10(10),log10(40),20);
cobs = dispR_surf96(vec_T,truemod,nmode); % "observations"
cstd = cobs * 0.01; % observation uncertainties

%% Starting model
% load starting model
% Get MINEOS model
CARD = ['./MINEOS_CARD/',cardname,'.card'];
[startmod, discs] = card2mod(CARD,zmax); % true model

vp_vs = startmod(:,2) ./ startmod(:,3); vp_vs(isinf(vp_vs))=0;
rho_vs = startmod(:,4) ./ startmod(:,3); rho_vs(isinf(rho_vs))=0;

% % Make homogeneous starting model
% zh2o = 1.618; % [km] water depth
% zmoho = 20; % [km] moho depth
% zmax = 200;
% z = [0 zh2o:2.5:zmoho zmoho+1:5:zmax];
% dz = [diff(z) 0];
% vs = 4.1*ones(size(dz)); vs(1)=0; % water
% vp = 1.75*vs; vp(1)=1.5;
% rho = vp / 2.5; rho(1)=1.03;
% startmod = [dz(:), vp(:), vs(:), rho(:)];
% % discs = [zmoho]; % [km] depth to sharp discontinuities
discs = []; % for this simple example, don't allow discontinuities

figure(1); clf;
box on; hold on;
h = plotlayermods(startmod(:,1),startmod(:,3),'-b');
h.LineWidth = 2;
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title('Starting Model');
set(gca,'FontSize',18,'linewidth',1.5);

%% Do linearized inversion using surf96
cstart = dispR_surf96(vec_T,startmod,nmode); % "predictions";

[finalmod, cpre, vs_std] = run_surf96_inv_Rayl_Vs_Vp_Rho(cobs,cstd,vec_T,startmod,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,z_dampbot,nit,nit_recalc_c,vp_vs,rho_vs,nmode);

%% Plot final kernels
figure(101); clf;
set(gcf,'position',[300         547        1281         422],'color','w')
Npers = length(vec_T);
clr = jet(Npers);
lgd = {};
[dcdvs, dcdvp, dudvs, dudvp, zkern, dcdrho, dudrho] = calc_kernel96(finalmod, vec_T, 'R', 1, 0,nmode);
subplot(1,3,1); box on; hold on;
for ip = 1:Npers
    plot(dcdvs(:,ip),zkern,'-','color',clr(ip,:),'linewidth',2); hold on;
    lgd{ip} = [num2str(vec_T(ip)),' s'];
end
xlabel('dc/dVs');
ylabel('Depth (km)');
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
% legend(lgd,'location','eastoutside');

subplot(1,3,2);
for ip = 1:Npers
    plot(dcdvp(:,ip),zkern,'-','color',clr(ip,:),'linewidth',2); hold on;
    lgd{ip} = [num2str(vec_T(ip)),' s'];
end
xlabel('dc/dVp');
ylabel('Depth (km)');
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
% legend(lgd,'location','eastoutside');

subplot(1,3,3);
for ip = 1:Npers
    plot(dcdrho(:,ip),zkern,'-','color',clr(ip,:),'linewidth',2); hold on;
    lgd{ip} = [num2str(vec_T(ip)),' s'];
end
xlabel('dc/d\rho');
ylabel('Depth (km)');
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
pos = get(gca,'Position');
legend(lgd,'location','eastoutside');
set(gca,'position',pos);

%%
% Plot
figure(100); clf; 
set(gcf,'position',[370   372   967   580]);

subplot(2,2,[1 3]); box on; hold on;
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
h = plotlayermods(startmod(:,1),startmod(:,3),'-b');
h.LineWidth = 2;
h = plotlayermods(finalmod(:,1),finalmod(:,3),'-r');
h.LineWidth = 2;
h = plotlayermods(finalmod(:,1),finalmod(:,3)-vs_std,'--r');
h.LineWidth = 2;
h = plotlayermods(finalmod(:,1),finalmod(:,3)+vs_std,'--r');
h.LineWidth = 2;
xlabel('Velocity');
ylabel('Depth');
set(gca,'FontSize',18,'linewidth',1.5);
legend({'true','start','final'},'Location','southwest')
% legend({'start','final'},'Location','southwest')

subplot(2,2,2); box on; hold on;
cpre = dispR_surf96(vec_T,finalmod,nmode);
plot(vec_T,cstart,'-ob','linewidth',2);
plot(vec_T,cpre,'-or','linewidth',2);
errorbar(vec_T,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
legend({'c start','c final','c obs'},'Location','southeast')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

