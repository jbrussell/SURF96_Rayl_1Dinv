function [outmodel, c_pre]= invdispR_lsqr(vec_T,c_obs,cstd,grv,grvstd,startmod,discs,par,maxdepth)
% Invert for Vs from Rayleigh wave phase velocities using surf96 and
% srfker96 codes (Vp and Density remain fixed). 
%
% Although surf96 has a surface wave inversion scheme built
% in, I've found that it deals with density in an unstable way and is
% difficult to control the general character of the models. This function
% calls surf96 to calculate phase velocities and srfker96 to calculated
% perturbation kernels, and an interative least squares inversion is
% carried out manually. Options are included for norm damping, first derivative
% smoothing, second derivative smoothing, and conserving layer velocity
% gradients. Layer boundaries can also be specified in order to break the
% constraint equations at desired depths.
%
%:INPUT: 
%	vec_T: vector of periods
%	c: phase velocity data (km/s)
%	cstd : error of phase velocity
%	grv: group velocity data (currently unused)
%	grvstd : error of group velocity (currently unused)
%	startmodel : nlayer x 4 array of data:
%		model(:,1): H(KM) layer thicknesses
%		model(:,2): VP(KM/S)
%		model(:,3): VS(KM/S)
%		model(:,4): DENSITY(g/cm^3)
%	discs: depth to discontinuities (km), if any. Constraint equations will
%          be broken at these depths
%   par: damping, smoothing, flatness, etc. parameters
%
%OUTPUT: 
%	outmodel: inverted 1D model nlayer x 3 array of data:
%	c_pre: predicted phase velocity dispersion (km/s) of final inverted model.
%
% jbrussell 6/5/2020
%

% Load parameters
eps_data = par.eps_data;% = 1e0; % data fit
eps_H = par.eps_H;% = 5e-2; % norm damping
eps_H_sedfac = par.eps_H_sedfac; % extra damping for sediments
eps_H_crustfac = par.eps_H_crustfac; % extra damping for sediments
eps_J = par.eps_J;% = 1e-1; % first derivative
eps_F = par.eps_F;% = 1e-1; % second derivative
eps_F_mantlefac = par.eps_F_mantlefac; % extra smoothing for mantle;
eps_K = par.eps_K;% = 1e0; % preserve layer velocity gradients;
nit = par.nit;% = 50; % number of least squares iterations
nit_recalc_c = par.nit_recalc_c;% = 10; % number of iteration before recalculating phase velocities
nit_recalc_kern = par.nit_recalc_kern;% = 25; % number of iterations before recalculating sensitivity kernels

if isempty(cstd)
    cstd = ones(size(vec_T));
end

disp('===========1D Rayleigh Dispersion Inversion============');

% Calculate perturbation kernels for G matrix
ifnorm = 0; % 1 for plotting only
ifplot = 0;
% [dcdvs, dcdvp, dudvs, dudvp, zkern] = calc_kernel96(startmod, vec_T, 'R', ifnorm, ifplot);    
[dcdvs, ~, ~, ~, ~] = calc_kernel96(startmod, vec_T, 'R', ifnorm, ifplot);    
G = dcdvs';

% Depth vector
z = cumsum(startmod(:,1));

% Set up smoothing and damping matrices
nlayer = length(startmod(:,1));

% Damping matrix
H00 = eye(nlayer);
h0 = startmod(:,3);
% constraint weighting matrix
wH = eye(size(H00))*eps_H;
% Pin water layer
ind_h2o = find(startmod(:,3)==0);
wH(ind_h2o,ind_h2o) = wH(ind_h2o,ind_h2o)*1e5;
% Add extra damping to seds
ind_discs = ismember(z,discs);
ind_noth2o = find(startmod(ind_discs,3)~=0);
ind_sed = ind_noth2o(1);
wH(ind_sed,ind_sed) = wH(ind_sed,ind_sed)*eps_H_sedfac;
% Add extra damping to crust (assumes last discontinuity is moho)
ind_crust = find(z<discs(end) & z>=discs(end-1));
wH(ind_crust,ind_crust) = wH(ind_crust,ind_crust)*eps_H_crustfac;
% Linear damp to ref model below maxdepth
dz = 50; % damp to background model below maxdepth + dz
ind_botmod_lin = find(z>=maxdepth & z<=maxdepth+dz);
wH(ind_botmod_lin,ind_botmod_lin) = wH(ind_botmod_lin,ind_botmod_lin) .* 0.1; %linspace(0.1,1,length(ind_botmod_lin));
ind_botmod_belowdz = find(z>maxdepth+dz);
wH(ind_botmod_belowdz,ind_botmod_belowdz) = wH(ind_botmod_belowdz,ind_botmod_belowdz)*100;

% first derviative flatness
J00 = build_flatness(nlayer);
j0 = zeros(size(J00,2),1);
% constraint weighting matrix
wJ = eye(size(J00))*eps_J;
wJ(ind_botmod_lin,ind_botmod_lin) = wJ(ind_botmod_lin,ind_botmod_lin)*0.1;

% second derivative smoothing
F00 = build_smooth( nlayer );
f0 = zeros(size(F00,2),1);
% constraint weighting matrix
wF = eye(size(F00))*eps_F;
% Add extra smoothing in mantle
ind_mantle = find(z<maxdepth & z>=discs(end));
wF(ind_mantle,ind_mantle) = wF(ind_mantle,ind_mantle)*eps_F_mantlefac;
%
wF(ind_botmod_lin,ind_botmod_lin) = wF(ind_botmod_lin,ind_botmod_lin)*10*10;

% preserve layer velocity gradients
K00 = build_flatness(nlayer);
k0 = [0; diff(startmod(:,3))];
% constraint weighting matrix
wK = eye(size(K00))*eps_K;
wK(ind_botmod_lin,ind_botmod_lin) = wK(ind_botmod_lin,ind_botmod_lin)*0;

% Break constraints at discontinuities
z_brks = discs;
J00 = break_constraint(J00, z, z_brks);
F00 = break_constraint(F00, z, z_brks);
K00 = break_constraint(K00, z, z_brks);

% combine all constraint equations
H = [wH*H00; wJ*J00; wF*F00; wK*K00];
h = [wH*h0;  wJ*j0;  wF*f0;  wK*k0];

% Initialize data vector
c_pre = dispR_surf96(vec_T,startmod);
dc = c_obs - c_pre;

% Data weighting
min_pct = 0.005; % minimum error percentage of observed
I_error_too_small = find(cstd./c_obs < min_pct);
cstd(I_error_too_small) = c_obs(I_error_too_small)*min_pct;
W = diag(1./cstd);

% Least squares inversion
premod = startmod;
vs_pre = startmod(:,3);
for ii = 1:nit
    
    % recalculate kernels?
    if mod(ii,nit_recalc_kern)==0
        [dcdvs, ~, ~, ~, ~] = calc_kernel96(premod, vec_T, 'R', 0, 0);    
        G = dcdvs';
    end
    % Recalculate phase velocities?
    if mod(ii,nit_recalc_c)==0
        c_pre = dispR_surf96(vec_T,premod);
        dc = c_obs - c_pre;
    end
    
    % reformulate inverse problem so that constraints apply directly to model
    d = dc + G*vs_pre;
    
    % least squares
    F = [W*G*eps_data; H];
    f = [W*d*eps_data; h];
    m = (F'*F)\F'*f;
    vs = m(1:nlayer);
    
    % Calculate model perturbation
    dvs = vs-vs_pre;
    
    % update data vector
    dc_pre = G * dvs;
    c_pre = c_pre + dc_pre;
    dc = c_obs - c_pre;
    
    % update model
    vs_pre = vs_pre + dvs;
    premod(:,3) = vs_pre;    
end
outmodel = premod;

disp('=================Inversion Finished========================');

end



