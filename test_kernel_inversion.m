% Test inversion using surf96 kernels
%
% jbrussell 6/5/2020
% modified 11/23/2021
%
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
eps_data = 1e0; % data fit
eps_H = 5e-2; % norm damping
eps_J = 1e-1; % first derivative smoothing
eps_F = 3e-1; % second derivative smoothing

% Other inversion parameters
nit = 50; % total number of iterations
nit_recalc_c = 10; % number of iterations after which to recalculate phase velocity and kernels

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

%% Starting model
% % Perturb the model to use as the starting model
% vs = truemod(:,3);
% vs_pre = vs;
% % vs_pre(5:11) = vs_pre(5:11)*1.05;
% vs_pre(12:25) = vs_pre(12:25)*1.05;
% startmod = truemod; startmod(:,3) = vs_pre;

% Make homogeneous starting model
zh2o = 1.618; % [km] water depth
zmoho = 20; % [km] moho depth
zmax = 150;
z = [0 zh2o:1:zmoho zmoho+1:5:zmax];
dz = [diff(z) 0];
vs = 4*ones(size(dz)); vs(1)=0; % water
vp = 1.75*vs; vp(1)=1.5;
rho = 3.3*ones(size(dz)); rho(1)=1.03;
startmod = [dz(:), vp(:), vs(:), rho(:)];
% discs = [zmoho]; % [km] depth to sharp discontinuities
discs = []; % for this simple example, don't allow discontinuities

%% Calculate kernels for G matrix using SURF96
ifnorm = 0; % for plotting only
ifplot = 0;
[dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(startmod, vec_T, 'R', ifnorm, ifplot);    
G = dcdvs';

% Set up smoothing and damping matrices
nlayer = length(startmod(:,1));
% Damping matrix
H00 = eye(nlayer);
h0 = startmod(:,3);
% Kill water layers
ind_h2o = find(startmod(:,3)==0);
H00(ind_h2o,ind_h2o) = 1e9;
% first derviative flatness
J00 = build_flatness(nlayer);
j0 = zeros(size(J00,2),1);
% second derivative smoothing
F00 = build_smooth( nlayer );
f0 = zeros(size(F00,2),1);

% Break constraints at discontinuities
% z_brks = [waterdepth, seddepth, mohodepth];
% z_brks = [ seddepth, mohodepth];
z_brks = discs;
z = cumsum(startmod(:,1));
J00 = break_constraint(J00, z, z_brks);
F00 = break_constraint(F00, z, z_brks);

% combine all constraints
H = [H00*eps_H; J00*eps_J; F00*eps_F];
h = [h0*eps_H; j0*eps_J; f0*eps_F];

% Data vector
cstart = dispR_surf96(vec_T,startmod); % "predictions";
cpre = cstart;
dc = cobs - cpre;

% Least squares inversion
premod = startmod;
vs_pre = premod(:,3);
for ii = 1:nit
    
    % reformulate inverse problem so that constraints apply directly to model
    d = dc + G*vs_pre;
    
    % last squares
    F = [G*eps_data; H];
    f = [d*eps_data; h];
    m = (F'*F)\F'*f;
    vs = m(1:nlayer);
    
    % update data vector
    dvs = vs-vs_pre;
    dcpre = G * dvs;
    cpre = cpre + dcpre;
    dc = cobs - cpre;
    
    % update model
    vs_pre = vs_pre + dvs;
    premod(:,3) = vs_pre;
    
    if mod(ii,nit_recalc_c)==0
        cpre = dispR_surf96(vec_T,premod);
        dc = cobs - cpre;
        
        ifnorm = 0; % for plotting only
        ifplot = 0;
        [dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(premod, vec_T, 'R', ifnorm, ifplot);    
        G = dcdvs';
    end
end

%% Plot final kernels
figure(101); clf;
Npers = length(vec_T);
clr = jet(Npers);
lgd = {};
[dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(premod, vec_T, 'R', 1, 0);
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

subplot(2,2,[1 3]); box on; hold on;
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
h = plotlayermods(startmod(:,1),startmod(:,3),'-b');
h.LineWidth = 2;
h = plotlayermods(premod(:,1),premod(:,3),'-r');
h.LineWidth = 2;
xlabel('Velocity');
ylabel('Depth');
set(gca,'FontSize',18,'linewidth',1.5);

subplot(2,2,2); box on; hold on;
plot(vec_T,cobs,'-ok');
plot(vec_T,cstart,'-ob');
plot(vec_T,cpre,'--or');
legend({'c obs','c start','c pre'},'Location','southeast')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

