% Test inversion using surf96 kernels
%
% jbrussell 6/5/2020
% modified 11/23/2021
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

z_dampbot = 150; % [km] damp to starting model below this depth

% Other inversion parameters
nit = 8; % total number of iterations
nit_recalc_c = 1; % number of iterations after which to recalculate phase velocity and kernels

%% Generate the synthetic dataset for this test
%Read in MINEOS model and convert to SURF96 layered model format
cardname = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_Line1_dist100.00km';

% Get MINEOS model
CARD = ['./MINEOS_CARD/',cardname,'.card'];
[truemod, discs] = card2mod(CARD,200); % true model
% % Keep track of discontinuities
% waterdepth = discs(1);
% seddepth = discs(2);
% mohodepth = discs(end);

% GENERATE SYNTHETIC DATASET
% Calculate dispersion for perturbed starting model and true model (which is like the data)
vec_T = logspace(log10(10),log10(40),20);
cobs = dispR_surf96(vec_T,truemod); % "observations"
cstd = cobs * 0.01; % observation uncertainties

%% Starting model
% Perturb the model to use as the starting model
vs = truemod(:,3);
vs_pre = vs;
% vs_pre(5:11) = vs_pre(5:11)*1.05;
% vs_pre(12:25) = vs_pre(12:25)*1.05;
vs_pre(2:end) = vs_pre(2:end)*0.9;
startmod = truemod; startmod(:,3) = vs_pre;

% % Make homogeneous starting model
% zh2o = 1.618; % [km] water depth
% zmoho = 20; % [km] moho depth
% zmax = 200;
% z = [0 zh2o:2.5:zmoho zmoho+1:5:zmax];
% dz = [diff(z) 0];
% vs = 4.1*ones(size(dz)); vs(1)=0; % water
% vp = 1.75*vs; vp(1)=1.5;
% rho = 3.3*ones(size(dz)); rho(1)=1.03;
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

%% Calculate kernels for G matrix using SURF96
ifnorm = 0; % for plotting only
ifplot = 0;
[dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(startmod, vec_T, 'R', ifnorm, ifplot);    
G = dcdvs';

% Data weighting
min_pct = 0.005; % minimum error percentage of observed
I_error_too_small = find(cstd./cobs < min_pct);
cstd(I_error_too_small) = cobs(I_error_too_small)*min_pct;
W = diag(1./cstd);

% Set up smoothing and damping matrices
nlayer = length(startmod(:,1));
dz = gradient(zkern);
dz_mat = repmat(dz,1,length(dz));
% Damping matrix
H00 = eye(nlayer);
h0 = startmod(:,3);
% first derviative flatness
J00 = build_flatness(nlayer) ./ dz_mat;
j0 = zeros(size(J00,2),1);
% second derivative smoothing
F00 = build_smooth( nlayer ) ./ dz_mat.^2;
f0 = zeros(size(F00,2),1);

% Break constraints at discontinuities
% z_brks = [waterdepth, seddepth, mohodepth];
% z_brks = [ seddepth, mohodepth];
z_brks = discs;
z = cumsum(startmod(:,1));
J00 = break_constraint(J00, z, z_brks);
F00 = break_constraint(F00, z, z_brks);

% Rescale the kernels
NA=norm(W*G,1);
NR=norm(H00,1);
eps_H0 = eps_H*NA/NR;
NR=norm(J00,1);
eps_J0 = eps_J*NA/NR;
NR=norm(F00,1);
eps_F0 = eps_F*NA/NR;

% Damp towards starting model
ind_dampstart = find(z > z_dampbot);
H00(ind_dampstart,ind_dampstart) = H00(ind_dampstart,ind_dampstart).*linspace(1,1000,length(ind_dampstart));
h0(ind_dampstart) = h0(ind_dampstart).*linspace(1,1000,length(ind_dampstart))';
% Kill water layers
ind_h2o = find(startmod(:,3)==0);
H00(ind_h2o,ind_h2o) = 1e9;

% combine all constraints
H = [H00*eps_H0; J00*eps_J0; F00*eps_F0];
h = [h0*eps_H0; j0*eps_J0; f0*eps_F0];

% Data vector
cstart = dispR_surf96(vec_T,startmod); % "predictions";
cpre = cstart;
dc = cobs - cpre;

% Least squares inversion
premod = startmod;
vs_pre = premod(:,3);
clrs = jet(nit);
isfigure = 1;
clear vs
for ii = 1:nit
    
    % reformulate inverse problem so that constraints apply directly to model
    dc = cobs - cpre;
    d = dc + G*vs_pre;
    
    % least squares
    F = [W*G*eps_data; H];
    f = [W*d*eps_data; h];
    m = (F'*F)\F'*f;
    vs = m(1:nlayer);
    
    % update data vector
    dvs = vs-vs_pre;
    dc_update = G * dvs;
    cpre = cpre + dc_update;
    
    % update model
    vs_pre = vs_pre + dvs;
    premod(:,3) = vs_pre;
    
    % Model uncertainties
    vs_std = diag(inv(F'*F)).^(1/2);
    
    if mod(ii,nit_recalc_c)==0
        cpre = dispR_surf96(vec_T,premod);
        dc = cobs - cpre;
        
        ifnorm = 0; % for plotting only
        ifplot = 0;
        [dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(premod, vec_T, 'R', ifnorm, ifplot);    
        G = dcdvs';
    end
    
    if isfigure
        if ii == 1 
            figure(2); clf; set(gcf,'position',[370   372   967   580]);
            
            subplot(2,2,[1 3]); box on; hold on;
            h2 = plotlayermods(startmod(:,1),startmod(:,3),'-k');
            h2.LineWidth = 4;
            
            subplot(2,2,2); box on; hold on;
            plot(vec_T,cstart,'-ok','linewidth',4);
        end
        subplot(2,2,[1 3]);
        h2 = plotlayermods(premod(:,1),premod(:,3),'-');
        h2.LineWidth = 2;
        h2.Color = clrs(ii,:);
        xlabel('Vs (km/s)');
        ylabel('Depth (km)');
        title('Starting Model');
        set(gca,'FontSize',18,'linewidth',1.5);
        
        subplot(2,2,2);
        plot(vec_T,cpre,'-o','color',clrs(ii,:),'linewidth',2);
        errorbar(vec_T,cobs,2*cstd,'or','linewidth',4);
        xlabel('Period');
        ylabel('Phase Velocity');
        set(gca,'FontSize',18,'linewidth',1.5);
        % pause;
        drawnow
    end
end

finalmod = premod;

%% Plot final kernels
figure(101); clf;
Npers = length(vec_T);
clr = jet(Npers);
lgd = {};
[dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(finalmod, vec_T, 'R', 1, 0);
for ip = 1:Npers
    plot(dcdvs(:,ip),zkern,'-','color',clr(ip,:),'linewidth',2); hold on;
    lgd{ip} = [num2str(vec_T(ip)),' s'];
end
xlabel('dc/dVs');
ylabel('Depth (km)');
title('Rayleigh-wave sensitivity kernels');
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
legend(lgd,'location','eastoutside');

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
cpre = dispR_surf96(vec_T,finalmod);
plot(vec_T,cstart,'-ob','linewidth',2);
plot(vec_T,cpre,'-or','linewidth',2);
errorbar(vec_T,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
legend({'c start','c final','c obs'},'Location','southeast')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

