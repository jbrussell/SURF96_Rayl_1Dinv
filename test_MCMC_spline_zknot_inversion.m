% Test simple Markov chain Monte Carlo Bayesian inversion for 3-layer
% velocity model. The method mostly follows Shen et al. (2013) GJI doi:10.1093/gji/ggs050
% This version uses smooth splines rather than layers and allows the splines to
% shift up and down depending on requirements of the data.
%
% jbrussell 10/11/2022
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

clear;
warning('off','all')
fullMAINpath = mfilename('fullpath');
functionspath = [fullMAINpath(1:regexp(fullMAINpath,mfilename)-1),'functions'];
addpath(functionspath);

% % Compile the faster mex files for spline calculation
% % !!!!! This only needs to be compiled the first time !!!!!
% cd ./functions
% CompileMexFiles
% cd ..

is_save_mat = 1;
is_resume_progress = 1; % Resume from "_PROGRESS" file?

PROJ = 'test'; % Project name 
%% MCMC parameters
% Other inversion parameters
nit_mcmc = 10000; % total number of iterations
nit_restart = 1000; %250; %1e10; % number of iterations after which to restart with new random model (if never want to restart, set to giant number)
N_cooldown = 100; %50; % number of iterations over which temperature parameter (tau) decays
m_perturb_method = 'single'; %'all'; % 'single' (perturb one model parameter at a time) | 'all' (perturb all at once)
m_perturb_method_zknot = 'single'; %'all'; % 'single' (perturb one model parameter at a time) | 'all' (perturb all at once)
nit_plot = 250; %1e9; % number of iterations after which to plot
nit_save = 20; % Number of iterations after which to store model output
nit_saveprogress = 500; % Number of iterations after which to save progress to mat file

% Define bounds of allowed model space M relative to ref. model. For the spline
% inversion, this applies to the Vs spline coefficients, not the layers.
% Models occuring outside this space will not be allowed.
% (these values also act as the min and max of the uniform prior)
% If a water layer exists, it is held at fixed velocity/density
par.dv_M = [-0.5 +0.5]; %[-0.25 +0.25]; % pct of reference model
par.dv_M_bottomknot = [-0.05 0.05]; % Control bounds of bottom most knot

% Define vectors for plotting
par.min_vs = 0; % minimum velocity
par.max_vs = 7; % maximum velocity
par.dvs_vec = 0.025; %0.01; % spacing between velocities

% Define widths of gaussian perturbations made at each iteration
par.dv_std = 0.05; %0.05; % km/s
par.dv_std_bottomknot = 0.005; % km/s; % control perturbation of bottom most knot
par.dzknot_std = 10; % km; % pertrubation of knot depths

% Scale vp and density with vs
par.vp_vs = 1.75; % Vp/Vs
par.rho_vs = 0.74; % density/Vs

% Spline parameters
Nspline = 5; % Number of desired splines, evenly spaced from surface (or base of water layer) to zmax
dz_int = 5; % (km) interpolated layer thicknesses. If too small, surf96 will break...
zmax = 200; % Maximum depth of starting model
zbot = 500; % maximum depth of entire model

outname = [strrep(PROJ,' ',''),'_Vs_bayesian'];

%% Generate the synthetic dataset for this test

% Load what we will consider the "true model"

% % Get MINEOS model
% %Read in MINEOS model and convert to SURF96 layered model format
% cardname = 'Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_Line1_dist100.00km';
% CARD = ['./MINEOS_CARD/',cardname,'.card'];
% [truemod, discs] = card2mod(CARD,200); % true model
% % % Keep track of discontinuities
% % waterdepth = discs(1);
% % seddepth = discs(2);
% % mohodepth = discs(end);

% (3-layer model)
zh2o_true = 1.618; % [km] water depth
zsed_true = 7; 
zmoho_true = 20; % [km] moho depth
% zlab = 70;
z_true = [0 zh2o_true zsed_true zmoho_true zmoho_true+40 zmoho_true+100 zbot];
dz = [diff(z_true) 0];
vs = [0 2.5 3.5 4.7 4.2 4.4 4.4]; 
% vp = 1.75*vs; vp(1)=1.5;
vp = par.vp_vs * vs; vp(vs==0)=1.5;
% rho = vp / 2.5; rho(1)=1.03;
rho = par.rho_vs * vs; rho(vs==0)=1.03;
truemod = [dz(:), vp(:), vs(:), rho(:)];

% GENERATE SYNTHETIC DATASET
% Calculate dispersion for true model, which will be our "observations"
% periods = logspace(log10(10),log10(40),10);
periods = logspace(log10(10),log10(150),15);
cobs = dispR_surf96(periods,truemod); % "observations"
cstd = cobs * 0.01; % observation uncertainties

%% Starting model

% 3-layer model
zh2o = 1.618; % [km] water depth
zsed = 7; 
zmoho = 20; % [km] moho depth
% zlab = 70;
z = [0 zh2o zsed zmoho zmax zbot];
dz = [diff(z) 0];
vs = [0 2.5-0.2 3.5-0.4 4.4+0.4 4.4 4.4]; 
% vs = [0 3.5 3.5 3.5 4.4]; 
% vp = 1.75*vs; vp(1)=1.5;
% vp = truemod(:,2);
vp = par.vp_vs * vs; vp(vs==0)=1.5;
% rho = vp / 2.5; rho(1)=1.03;
% rho = truemod(:,4);
rho = par.rho_vs * vs; rho(vs==0)=1.03;
refmod = [dz(:), vp(:), vs(:), rho(:)];

%%
zmin = zh2o;

% Interpolate layered models
% zinterp = [0 zh2o linspace(zh2o,zmin,20) zmin:dz_int:zbot]';
zinterp = [0 zh2o zh2o:dz_int:zbot]';
% zsp = z(:);
% zsp = [0:10:zmoho zmoho+50:50:zmax]; 
[mod_true] = layerizemod_interp(truemod,zinterp);
[mod_ref] = layerizemod_interp(refmod,zinterp);
% [mod_true] = layermod2node_interp(truemod,zinterp);
% [mod_ref] = layermod2node_interp(refmod,zinterp);

% Do spline calculations for reference (and true) model
% zsp = [zh2o:50:zmax];
% zsp = [linspace(zh2o,zmoho,5-1) linspace(zmoho,zmax,6-1)]; % example of custom spline spacing with a discontinuity
zsp = linspace(zmin,zmax,Nspline-1);
Inoh2o = find(mod_true.z>=zmin & mod_true.z<=zmax); Inoh2o = Inoh2o(2:end);
Ih2o = find(mod_true.z<=zmin); Ih2o = Ih2o(1:end-1);
Ifixed = find(mod_true.z > zmax);
[spbasis_true,spcoeffs_true,spzz_true]=make_splines(zsp(:),[],mod_true.z(Inoh2o),mod_true.vs(Inoh2o));
vs_true_sp = spbasis_true * spcoeffs_true;
vs_true_sp = [mod_true.vs(Ih2o); vs_true_sp; mod_true.vs(Ifixed)];
Inoh2o = find(mod_ref.z>=zmin & mod_ref.z<=zmax); Inoh2o = Inoh2o(2:end);
Ih2o = find(mod_ref.z<=zmin); Ih2o = Ih2o(1:end-1);
Ifixed = find(mod_ref.z > zmax);
[spbasis,spcoeffs,spzz]=make_splines(zsp(:),[],mod_ref.z(Inoh2o),mod_ref.vs(Inoh2o));
vs_ref_sp = spbasis * spcoeffs;
vs_ref_sp = [mod_ref.vs(Ih2o); vs_ref_sp; mod_ref.vs(Ifixed)];

% Convert spline to mod structure
[refmod_sp] = spline2mod_zmin(mod_ref,vs_ref_sp,par.vp_vs,par.rho_vs,zmin);
refmod_sp_z = [0; cumsum(refmod_sp(1:end-1,1))];
[truemod_sp] = spline2mod_zmin(mod_true,vs_true_sp,par.vp_vs,par.rho_vs,zmin);

figure(1); clf;
set(gcf,'position',[370   372   967   580]);
subplot(2,2,[1 3]); box on; hold on;
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
plot(vs_true_sp,zinterp,'--k','linewidth',1.5);
h = plotlayermods(refmod(:,1),refmod(:,3),'-b');
h.LineWidth = 2;
plot(vs_ref_sp,zinterp,'--b','linewidth',1.5);
plot(spbasis'*2,mod_ref.z(Inoh2o));
plot(zeros(size(zsp)),zsp,'ok');
xlabel('Velocity');
ylabel('Depth');
set(gca,'FontSize',18,'linewidth',1.5);
legend({'true','true (spline)','start','start (spline)'},'Location','southwest')
% legend({'start','final'},'Location','southwest')

subplot(2,2,2); box on; hold on;
cref = dispR_surf96(periods,refmod); % predicted phase velocity
cobs_sp = dispR_surf96(periods,truemod_sp); % predicted phase velocity
cref_sp = dispR_surf96(periods,refmod_sp); % predicted phase velocity
errorbar(periods,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
plot(periods,cobs_sp,'--ok','linewidth',2);
plot(periods,cref,'-ob','linewidth',2);
plot(periods,cref_sp,'--b','linewidth',1.5);
legend({'c obs','c obs (spline)','c start','c start (spline)'},'Location','southeast')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

figure(10); clf; box on; hold on;
h = plotlayermods(truemod(:,1),truemod(:,3),'-b');
h.LineWidth = 2;
plot(vs_true_sp,zinterp,'--b','linewidth',1.5);
h = plotlayermods(truemod_sp(:,1),truemod_sp(:,3),'-c');
h.LineWidth = 2;
h = plotlayermods(refmod(:,1),refmod(:,3),'-r');
h.LineWidth = 2;
plot(vs_ref_sp,zinterp,'--r','linewidth',1.5);
h = plotlayermods(refmod_sp(:,1),refmod_sp(:,3),'-m');
h.LineWidth = 2;


%% Define priors for each layer
% Define edges of the model space M
Ncoeffs = length(spcoeffs);
model_bounds = nan(length(spcoeffs),2);
for ic = 1:Ncoeffs
    model_bounds(ic,1) = spcoeffs(ic)*(1+par.dv_M(1));
    model_bounds(ic,2) = spcoeffs(ic)*(1+par.dv_M(2));
end
model_bounds(Ncoeffs,1) = spcoeffs(ic)*(1+par.dv_M_bottomknot(1));
model_bounds(Ncoeffs,2) = spcoeffs(ic)*(1+par.dv_M_bottomknot(2));

% Uniform priors spanning M
priors.sample = @(N,ic) unifrnd(model_bounds(ic,1), model_bounds(ic,2) ,N,1);

% Function to perturb model
perturb_model = @(model,std_vec) normrnd(model(:)',std_vec)';
perturb_zknot = @(model,std_vec) normrnd(model(:)',std_vec)';

% Get pdf from distributions
vs_edges = [par.min_vs:par.dvs_vec:par.max_vs];
vs_vec = 0.5*(vs_edges(1:end-1)+vs_edges(2:end));
figure(1000);
for ic = 1:Ncoeffs
    priors.pdf_sp{ic} = 1./(model_bounds(ic,2)-model_bounds(ic,1)) * ones(size(vs_vec));
    priors.pdf_sp{ic}(vs_vec>model_bounds(ic,2)) = 0;
    priors.pdf_sp{ic}(vs_vec<model_bounds(ic,1)) = 0;
end
priors.vs_vec = vs_vec;

figure(999); clf;
for ic = 1:Ncoeffs
    plot(priors.vs_vec,priors.pdf_sp{ic}); hold on;
end
title('Priors on Coefficients');

% Project priors to layer space using spline basis
pdf_mat_sp = [];
for ic = 1:Ncoeffs
    pdf_mat_sp(ic,:) = priors.pdf_sp{ic};
end
pdf_mat = spbasis*pdf_mat_sp;
for ilay = 1:size(pdf_mat,1)
    priors.pdf{ilay} = pdf_mat(ilay,:);
end
figure(11); clf;
plot(priors.vs_vec,pdf_mat')
title('Priors on Layers');

%% Do MCMC
Nmodels = ceil(nit_mcmc / nit_save);
posterior = nan(size(refmod_sp,1),Nmodels);
posterior_sp = nan(Ncoeffs,Nmodels);
cpre = nan(length(cobs),Nmodels);
misfit = nan(1,Nmodels);
Likelihood = nan(1,Nmodels);
vs_models = nan(size(refmod_sp,1),Nmodels);
vs_models_sp = nan(Ncoeffs,Nmodels);
models = nan([size(refmod_sp),Nmodels]);
zsp_mat = nan(Nspline-1,Nmodels);

% Initiate
m_j = spcoeffs;
m_j(:,3) = sample_model(priors.sample,1,Ncoeffs);
m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
spbasis_j = spbasis;
zsp_j = zsp;
ii = 0;
ibad = 0;
ii_cooldown = 0;
ii_save = 0;
tic

outtemp = ['./bayesian_mcmc_Vs_spline_zknot/',outname,'_Nspline',num2str(Nspline),'_PROGRESS.mat'];
if exist(outtemp) && is_resume_progress==1
    load(outtemp);
end

while ii < nit_mcmc
    
    if ii>0 && mod(ii,nit_restart) == 0 % reinitialize mcmc, start over
        m_j(:,3) = sample_model(priors.sample,1,Ncoeffs);
        m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
        m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
        ii_cooldown = 0;
    end
    
    % Previous model
    %     Inoh2o = find(mod_ref.vs~=0);
    %     Ih2o = find(mod_ref.vs==0);
        Inoh2o = find(mod_ref.z>=zmin & mod_ref.z<=zmax); Inoh2o = Inoh2o(2:end);
        Ih2o = find(mod_ref.z<=zmin); Ih2o = Ih2o(1:end-1);
        Ifixed = find(mod_ref.z > zmax);
        vs_spline = spbasis_j * m_j(:,3);
        vs_spline = [mod_ref.vs(Ih2o); vs_spline; mod_ref.vs(Ifixed)];
        [splinemod_j] = spline2mod_zmin(mod_ref,vs_spline,par.vp_vs,par.rho_vs,zmin);
    c_j = dispR_surf96(periods,splinemod_j); % predicted phase velocity
    if length(c_j) ~= length(periods) % check if something is wrong...
        ibad = ibad+1;
        m_j(:,3) = sample_model(priors.sample,1,Ncoeffs);
        m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
        m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
        continue
    end
    S_j = sum((cobs(:)-c_j(:)).^2./cstd(:).^2); % misfit
    L_j = ((2 * pi)^(length(periods)) * prod(cstd(:).^2)).^(-0.5) .* exp(-0.5 * S_j); % likelihood
%     L_j = exp(-0.5 * S_j); % likelihood
    
    % Ensure that model is within model space M
    is_in_bounds = is_model_in_bounds(m_j,model_bounds);

    % If model is really bad, try a new one
    % if L_j < eps || isnan(L_j) || ~is_in_bounds
    if isinf(1./L_j) || isnan(L_j) || ~is_in_bounds
        ibad = ibad+1;
        m_j(:,3) = sample_model(priors.sample,1,Ncoeffs);
        m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
        m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
        dzknot = perturb_zknot(zsp,repmat(par.dzknot_std,1,length(zsp)));
        dzknot(dzknot<min(zsp)) = min(zsp);
        dzknot(dzknot>max(zsp)) = max(zsp);
        dzknot(1) = min(zsp);
        dzknot(end) = max(zsp);
        zsp_j = sort(dzknot);
        Inoh2o = find(mod_ref.z>=zmin & mod_ref.z<=zmax); Inoh2o = Inoh2o(2:end);
        [spbasis_j,~,~]=make_splines(zsp_j(:),[],mod_ref.z(Inoh2o),splinemod_j(Inoh2o,3));
        display(['Searching for stable starting model: ',num2str(ibad)]);
        continue
    end
    ii = ii + 1;
    
    if mod(ii,100) == 0
        display([num2str(ii),'/',num2str(nit_mcmc)]);
    end
    
    % Store output
    if ii>0 && mod(ii,nit_save) == 0
        ii_save = ii_save + 1;
        
        % Calculate posterior probability of model j (spline coefficients)
        for ic = 1:Ncoeffs
            [~,I] = min(abs(m_j(ic,3)-priors.vs_vec));
            posterior_sp(ic,ii_save) = L_j .* priors.pdf_sp{ic}(I);
        end
        % Calculate posterior for layered structure
        % Project priors to layer space using spline basis
        pdf_mat_sp = [];
        for ic = 1:Ncoeffs
            pdf_mat_sp(ic,:) = priors.pdf_sp{ic};
        end
        pdf_mat = spbasis_j*pdf_mat_sp;
        for ilay = 1:size(pdf_mat,1)
            priors.pdf{ilay} = pdf_mat(ilay,:);
        end
        ipdf = 0;
        for ilay = 1:size(posterior,1)
    %         if splinemod_j(ilay,3)==0 % water layer
            if refmod_sp_z(ilay)<zmin || refmod_sp_z(ilay)>zmax-dz_int
                posterior(ilay,ii_save) = L_j * 1;
                continue
            end
            ipdf = ipdf + 1;
            [~,I] = min(abs(splinemod_j(ilay,3)-priors.vs_vec));
            posterior(ilay,ii_save) = L_j .* priors.pdf{ipdf}(I);
        end
    %     posterior(:,ii) = L_j;
        
        % Save outputs
        misfit(ii_save) = S_j;
        Likelihood(ii_save) = L_j;
        cpre(:,ii_save) = c_j(:);
        vs_models(:,ii_save) = splinemod_j(:,3);
        vs_models_sp(:,ii_save) = m_j(:,3);
        models(:,:,ii_save) = splinemod_j;
        zsp_mat(:,ii_save) = zsp_j;
    end
    
    % Decaying thermal parameter (cool down parameter) from simulated 
    % annealing (Kirkpatrick et al. 1983). This allows larger changes 
    % between sequential models at early iterations. This premultiplies the
    % Gaussian distributions from which random model parameters are drawn
    % and also the likelihood of the trial model, so misfit increases are
    % more likely accepted early in the MCMC.
%     tau = 1 + 3 * erfc(ii/500); % denom = 500 means decays over ~1500 iterations (Eilon et al. 2018)
    tau = 1 + 3 * erfc(ii_cooldown/(N_cooldown/3)); % denom = 500 means decays over ~1500 iterations
    
    % Trial model
    is_in_bounds = 0;
    iiloop = 0;
    iirepeat = 0;
    is_restart = 0;
    % Get index for type of perturbation to perform
    I_perturbation_type = ceil(rand(1)*2);
    while is_in_bounds == 0
        m_i = m_j;
        zsp_i = zsp_j;
        spbasis_i = spbasis_j;
        
        if iiloop > 1e6
            % Get new index for type of perturbation to perform
            I_perturbation_type = ceil(rand(1)*2);
            iiloop = 0;
            iirepeat = iirepeat + 1;
            display(['Giving up on type ',num2str(I_perturbation_type),' ...'])
        end
        
        if iirepeat > 10
            % If get stuck, just restart from the beginning
            is_restart = 1;
            break
        end
        
        if I_perturbation_type == 1 % PERTURB VALUE OF COEFFICIENT
            std_pert = tau*repmat(par.dv_std,1,Ncoeffs);
            std_pert(Ncoeffs) = par.dv_std_bottomknot;
            dvs = perturb_model(m_i(:,3),std_pert); % perturb Vs
        %     dvs = sample_model(priors.sample,1,Ncoeffs);; % random Vs
            switch m_perturb_method
                case 'single'
                    I_pert = ceil(rand(1)*Ncoeffs); % randomly pick model parameter to perturb
                    m_i(I_pert,3) = dvs(I_pert);
                case 'all'
                    m_i(:,3) = dvs; % perturb all model parameters at once
                otherwise
                    error('m_perturb_method not a valid choice. must be ''single'' or ''all'' ');
            end
            m_i(:,2) = par.vp_vs*m_i(:,3); m_i(m_i(:,3)==0,2)=1.5;
            m_i(:,4) = par.rho_vs*m_i(:,3); m_i(m_i(:,3)==0,4)=1.03;
        elseif I_perturbation_type == 2 % PERTURB DEPTH OF KNOT
            dzknot = perturb_zknot(zsp_i,tau*repmat(par.dzknot_std,1,length(zsp_i)));
%             dzknot(dzknot<zh2o) = zh2o;
%             dzknot(dzknot>zmax_sp) = zmax_sp;
            switch m_perturb_method_zknot
                case 'single'
                    I_pert = randi([2,length(zsp_i)-1]); % randomly pick knot index (but avoid top and bottom knots)
                    if dzknot(I_pert)<min(zsp) || dzknot(I_pert)>max(zsp)
                        is_in_bounds = 0;
                        continue
                    end
                    zsp_i(I_pert) = dzknot(I_pert);
                case 'all'
                    if ~isempty(dzknot(dzknot<min(zsp))) || ~isempty(dzknot(dzknot>max(zsp)))
                        is_in_bounds = 0;
                        continue
                    end
                    zsp_i = dzknot; % perturb all knots at once
                otherwise
                    error('m_perturb_method not a valid choice. must be ''single'' or ''all'' ');
            end
            zsp_i = sort(zsp_i);
            Inoh2o = find(mod_ref.z>=zmin & mod_ref.z<=zmax); Inoh2o = Inoh2o(2:end);
            [spbasis_i,m_i(:,3),~]=make_splines(zsp_i(:),[],mod_ref.z(Inoh2o),splinemod_j(Inoh2o,3));
            m_i(:,2) = par.vp_vs*m_i(:,3); m_i(m_i(:,3)==0,2)=1.5;
            m_i(:,4) = par.rho_vs*m_i(:,3); m_i(m_i(:,3)==0,4)=1.03;
        end
        is_in_bounds = is_model_in_bounds(m_i,model_bounds);
        
        iiloop = iiloop + 1;
    end
    if is_restart
        ii = ii - 1;
        continue
    end
    Inoh2o = find(mod_ref.z>=zmin & mod_ref.z<=zmax); Inoh2o = Inoh2o(2:end);
    Ih2o = find(mod_ref.z<=zmin); Ih2o = Ih2o(1:end-1);
    Ifixed = find(mod_ref.z > zmax);
    vs_spline = spbasis_i * m_i(:,3);
    vs_spline = [mod_ref.vs(Ih2o); vs_spline; mod_ref.vs(Ifixed)];
    [splinemod_i] = spline2mod_zmin(mod_ref,vs_spline,par.vp_vs,par.rho_vs,zmin);
    c_i = dispR_surf96(periods,splinemod_i); % predicted phase velocity
    if length(c_i) ~= length(periods) % check if something is wrong...
        % skip
        continue
    end
    S_i = sum((cobs(:)-c_i(:)).^2./cstd(:).^2); % misfit
    L_i = ((2 * pi)^(length(periods)) * prod(cstd(:).^2)).^(-0.5) .* exp(-0.5 * S_i); % likelihood
%     L_i = exp(-0.5 * S_i); % likelihood
    L_i = tau * L_i;
    
    % Plot
    if mod(ii,nit_plot) == 0
        figure(2); clf;
        subplot(2,2,1); box on; hold on;
        yyaxis left
        plot(1:ii_save,misfit(1:ii_save) / length(periods),'o'); hold on;
        ylabel('Misfit');        
        yyaxis right
        plot(1:ii_save,log10(Likelihood(1:ii_save)),'o'); hold on;
        ylabel('log_{10}(Likelihood)');
        
        subplot(2,2,[2 4]); box on; hold on;
        for kk = 1:ii_save
            h = plotlayermods(models(:,1,kk),models(:,3,kk),'-r');
            h.LineWidth = 1;
        end
        h = plotlayermods(refmod(:,1),refmod(:,3),'-b');
        h.LineWidth = 2;
        xlabel('Vs (km/s)');
        ylabel('Depth (km)');
        set(gca,'FontSize',16,'linewidth',1.5);
%         legend({'start','ensemble'},'Location','southwest')

        subplot(2,2,3); box on; hold on;
        plot(periods,cpre(:,1:ii_save),'-or','linewidth',1);
        errorbar(periods,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
        plot(periods,cref,'-ob','linewidth',2);
        xlabel('Period');
        ylabel('Phase Velocity');
        set(gca,'FontSize',16,'linewidth',1.5);
        drawnow;
    end
    
    % Metropolis-Hastings acceptance criterion
    p_accept = min(L_i/L_j, 1);
    if rand <= p_accept % (rand always between [0 1])
        % Accept new model i
        m_j = m_i;
        spbasis_j = spbasis_i;
        zsp_j = zsp_i;
    else
        % Reject new model i
        % continue
    end
    
    % Save progress
    if mod(ii,nit_saveprogress) == 0
        if ~exist('./bayesian_mcmc_Vs_spline_zknot/')
            mkdir('./bayesian_mcmc_Vs_spline_zknot/');
        end
        outtemp = ['./bayesian_mcmc_Vs_spline_zknot/',outname,'_Nspline',num2str(Nspline),'_PROGRESS.mat'];
        save(outtemp,'posterior','posterior_sp','misfit','Likelihood','cpre','vs_models','vs_models_sp','models','zsp_mat','ii','ii_save');
    end
end
toc

%% Calculate marginal pdfs
vs_edges = [par.min_vs:par.dvs_vec:par.max_vs];
vs_vec = 0.5*(vs_edges(1:end-1)+vs_edges(2:end));
marginal_pdf = zeros(size(vs_models,1),length(vs_vec));
marginal_pdf_sp = zeros(Ncoeffs,length(vs_vec));

% Marginal for spline coefficients
for idim = 1:Ncoeffs
    ind_bin = discretize(vs_models_sp(idim,:),vs_edges);
    marginal = zeros(size(vs_vec));
    for ii = 1:length(ind_bin)
        if isnan(ind_bin(ii))
            continue
        end
        marginal(ind_bin(ii)) = marginal(ind_bin(ii)) + sum(posterior_sp(:,ii));
    end
    marginal_pdf_sp(idim,:) = marginal / sum(marginal); % normalize so sums to 1
end

% % Expand with basis function
% for ii = 1:size(marginal_pdf_sp,2)
%     marginal_pdf(2:end,ii) = spbasis * marginal_pdf_sp(:,ii);
% end

% Marginal for depth model
for idim = 1:size(vs_models,1)
    ind_bin = discretize(vs_models(idim,:),vs_edges);
    marginal = zeros(size(vs_vec));
    for ii = 1:length(ind_bin)
        if isnan(ind_bin(ii))
            continue
        end
        marginal(ind_bin(ii)) = marginal(ind_bin(ii)) + sum(posterior(:,ii));
    end
    marginal_pdf(idim,:) = marginal / sum(marginal); % normalize so sums to 1
end
z = [0; cumsum(refmod_sp(1:end-1,1))];
[vsgrid,zgrid] = meshgrid(vs_vec,z);

% 2-D MARGINAL PDFs
% Convert from layers defined by a center point and width to knots like mineos
z_lays = [];
z_lays(1,1) = z(1);
z_lays(2,1) = z(2);
icnt = 2;
for ii = 2:length(z)-1
    icnt = icnt + 1;
    z_lays(icnt,1) = z(ii);
    icnt = icnt + 1;
    z_lays(icnt,1) = z(ii+1);
end
Z_lays = repmat(z_lays,1,size(marginal_pdf,2));
VS = repmat(vs_vec,size(Z_lays,1),1);
marginal_pdf_lays = zeros(length(z_lays),size(marginal_pdf,2));
for ivs = 1:size(marginal_pdf,2)
    icnt = 0;
    for ii = 1:size(marginal_pdf,1)-1
        icnt = icnt + 1;
        marginal_pdf_lays(icnt,ivs) = marginal_pdf(ii,ivs);
        icnt = icnt + 1;
        marginal_pdf_lays(icnt,ivs) = marginal_pdf(ii+1,ivs);
    end
end

% Marginal for depth of spline knots
z_int_edges = [0:10:zbot];
z_int_vec = 0.5*(z_int_edges(1:end-1)+z_int_edges(2:end));
marginal_pdf_zsp = zeros(size(zsp_mat,1),length(z_int_vec));
for idim = 1:size(zsp_mat,1)
    ind_bin = discretize(zsp_mat(idim,:),z_int_edges);
    marginal = zeros(size(z_int_vec));
    for ii = 1:length(ind_bin)
        if isnan(ind_bin(ii))
            continue
        end
        marginal(ind_bin(ii)) = marginal(ind_bin(ii)) + sum(posterior_sp(:,ii));
    end
    marginal_pdf_zsp(idim,:) = marginal / sum(marginal); % normalize so sums to 1
end

% Marginal for phv predictions
par.dphv_vec = 0.01; % km/s 
% phv_edges = [3.8:par.dphv_vec:4.4];
phv_edges = [par.min_vs:par.dphv_vec:par.max_vs];
phv_vec = 0.5*(phv_edges(1:end-1)+phv_edges(2:end));
marginal_pre_pdf = zeros(size(cpre,1),length(phv_vec));
for idim = 1:size(cpre,1)
    ind_bin = discretize(cpre(idim,:),phv_edges);
    marginal = zeros(size(phv_vec));
    for ii = 1:length(ind_bin)
        if isnan(ind_bin(ii))
            continue
        end
        marginal(ind_bin(ii)) = marginal(ind_bin(ii)) + sum(posterior(:,ii));
    end
    marginal_pre_pdf(idim,:) = marginal / sum(marginal); % normalize so sums to 1
end
[phv_grid,periods_grid] = meshgrid(phv_vec,periods);

% Confidence Fields
vs_conf_mat = nan(size(marginal_pdf));
for ilay = 1:size(marginal_pdf,1)
    vs_conf_mat(ilay,:) = confidence_field(marginal_pdf(ilay,:));
end

phv_conf_mat = nan(size(marginal_pre_pdf));
for ip = 1:size(marginal_pre_pdf,1)
    phv_conf_mat(ip,:) = confidence_field(marginal_pre_pdf(ip,:));
end

%% Gather information

bayesian.PROJ = PROJ;
bayesian.params.zmax = zmax;
bayesian.params.zbot = zbot;
bayesian.params.par = par;
% bayesian.params.param = param;
bayesian.zh2o = zh2o;

bayesian.refmod = refmod;
bayesian.refmod_sp = refmod_sp;
bayesian.models_mat = models;
bayesian.vs_mat_sp = vs_models_sp;
bayesian.vs_mat = vs_models;
bayesian.zsp_mat = zsp_mat;

bayesian.priors = priors;
bayesian.zsp = zsp;

bayesian.post.Post_mat = posterior;
bayesian.post.Post_mat_sp = posterior_sp;
bayesian.post.vs_pdf_lays = marginal_pdf_lays;
bayesian.post.vs_pdf = marginal_pdf;
bayesian.post.vs_pdf_sp = marginal_pdf_sp;
bayesian.post.phv_pdf = marginal_pre_pdf;
bayesian.post.zsp_pdf = marginal_pdf_zsp;
bayesian.post.vs_conf_mat = vs_conf_mat;
bayesian.post.phv_conf_mat = phv_conf_mat;
% bayesian.post.disc_pdf = disc_pdf;
bayesian.vs_mat_lays = VS;
bayesian.z_mat_lays = Z_lays;
bayesian.vsgrid = vsgrid;
bayesian.zgrid = zgrid;
bayesian.z_int = z;
bayesian.vs_vec = vs_vec;
bayesian.vs_edges = vs_edges;
bayesian.z_int_vec = z_int_vec;
bayesian.phv_pre_mat = cpre;

bayesian.periods = periods;
bayesian.periods_grid = periods_grid;
bayesian.phv_grid = phv_grid;
bayesian.phv_edges = phv_edges;
bayesian.phv_vec = phv_vec;
bayesian.cobs = cobs;
bayesian.cstd = cstd;
% bayesian.data = mat.data;

bayesian.Likelihood = Likelihood;
bayesian.misfit = misfit;
% [chi2_mat_srt, isrt] = sort(bayesian.chi2_mat);
% igood = sort(isrt(1:bayesian.params.N_bestfitting));
% bayesian.qmu_inv_mat_good = qmu_inv_mat(:,igood);
% bayesian.qinv_pre_mat_good = bayesian.qinv_pre_mat(:,igood);
[~,imin] = min(bayesian.misfit);
bayesian.vs_mat_best = bayesian.vs_mat(:,imin);
bayesian.phv_pre_mat_best = bayesian.phv_pre_mat(:,imin);

% Posterior probabilities (histograms)
w = sum(posterior,1);
bayesian.post.vs_mean = sum(w.*bayesian.vs_mat,2)./sum(w);
% bayesian.post.qmu_inv_mean = sum(bayesian.post.qmu_inv_pdf.*qinvgrid,2);
% bayesian.post.phv_mean_pre = GG * [bayesian.post.sum(posterior,1)_mean(:); qkap_inv(:)];
bayesian.post.zsp_mean = sum(w.*bayesian.zsp_mat,2)./sum(w);
% Median
bayesian.post.vs_med = pdf_prctile(bayesian.post.vs_pdf,bayesian.vs_vec+par.dvs_vec/2,50);
bayesian.post.vs_med(isnan(bayesian.post.vs_med)) = bayesian.post.vs_mean(isnan(bayesian.post.vs_med));
bayesian.post.vs_l95 = pdf_prctile(bayesian.post.vs_pdf,bayesian.vs_vec+par.dvs_vec/2,2.5);
bayesian.post.vs_l95(isnan(bayesian.post.vs_l95)) = bayesian.post.vs_mean(isnan(bayesian.post.vs_l95));
bayesian.post.vs_u95 = pdf_prctile(bayesian.post.vs_pdf,bayesian.vs_vec+par.dvs_vec/2,97.5);
bayesian.post.vs_u95(isnan(bayesian.post.vs_u95)) = bayesian.post.vs_mean(isnan(bayesian.post.vs_u95));
bayesian.post.vs_l68 = pdf_prctile(bayesian.post.vs_pdf,bayesian.vs_vec+par.dvs_vec/2,16);
bayesian.post.vs_l68(isnan(bayesian.post.vs_l68)) = bayesian.post.vs_mean(isnan(bayesian.post.vs_l68));
bayesian.post.vs_u68 = pdf_prctile(bayesian.post.vs_pdf,bayesian.vs_vec+par.dvs_vec/2,84);
bayesian.post.vs_u68(isnan(bayesian.post.vs_u68)) = bayesian.post.vs_mean(isnan(bayesian.post.vs_u68));
% bayesian.post.qinv_med_pre = GG * [bayesian.post.qmu_inv_med(:); qkap_inv(:)];
% bayesian.post.qinv_l95_pre = GG * [bayesian.post.qmu_inv_l95(:); qkap_inv(:)];
% bayesian.post.qinv_u95_pre = GG * [bayesian.post.qmu_inv_u95(:); qkap_inv(:)];
% bayesian.post.qinv_l68_pre = GG * [bayesian.post.qmu_inv_l68(:); qkap_inv(:)];
% bayesian.post.qinv_u68_pre = GG * [bayesian.post.qmu_inv_u68(:); qkap_inv(:)];
bayesian.post.phv_med_pre = (pdf_prctile(bayesian.post.phv_pdf,bayesian.phv_vec+par.dphv_vec/2,50));
bayesian.post.phv_l95_pre = (pdf_prctile(bayesian.post.phv_pdf,bayesian.phv_vec+par.dphv_vec/2,2.5));
bayesian.post.phv_u95_pre = (pdf_prctile(bayesian.post.phv_pdf,bayesian.phv_vec+par.dphv_vec/2,97.5));
bayesian.post.phv_l68_pre = (pdf_prctile(bayesian.post.phv_pdf,bayesian.phv_vec+par.dphv_vec/2,16));
bayesian.post.phv_u68_pre = (pdf_prctile(bayesian.post.phv_pdf,bayesian.phv_vec+par.dphv_vec/2,84));

bayesian.nit_save = nit_save;

if is_save_mat
    if ~exist('./bayesian_mcmc_Vs_spline_zknot/')
        mkdir('./bayesian_mcmc_Vs_spline_zknot/');
    end
    outmat = ['./bayesian_mcmc_Vs_spline_zknot/',outname,'_Nspline',num2str(Nspline),'.mat'];
    save(outmat,'bayesian');
    
    delete(outtemp);
end

%% Plot misfit/likelihood evolution
figure(888); clf;
subplot(2,1,1);
plot([1:Nmodels] * bayesian.nit_save,bayesian.misfit / length(bayesian.periods),'o'); hold on;
xlabel('Model #');
ylabel('Misfit');

subplot(2,1,2);
plot([1:Nmodels] * bayesian.nit_save,log10(bayesian.Likelihood),'o'); hold on;
xlabel('Model #');
ylabel('log_{10}(Likelihood)');

%% Plot Vs models
figure(100); clf; 
set(gcf,'position',[370   372   967   580]);

subplot(2,2,[1 3]); box on; hold on;
for ii = 1:size(models,3)
%     if misfit(ii)/length(periods) > 2
%         continue
%     end
    h = plotlayermods(models(:,1,ii),models(:,3,ii),'-r');
    h.LineWidth = 1;
    if ii == 1
        h1(ii) = h;
    end
end
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
h1(2) = h;
h = plotlayermods(refmod(:,1),refmod(:,3),'-b');
h.LineWidth = 2;
h1(3) = h;
[~,imin] = min(misfit);
h = plotlayermods(models(:,1,imin),models(:,3,imin),'--g');
h.LineWidth = 2;
% w = sum(posterior,1);
% vs_med = sum(w.*vs_models,2)./sum(w);
% dz_med = sum(w.*squeeze(models(:,1,:)),2)./sum(w);
% h = plotlayermods(dz_med,vs_med,'--g');
% h.LineWidth = 2;
h1(4) = h;
xlabel('Velocity');
ylabel('Depth');
set(gca,'FontSize',18,'linewidth',1.5);
legend(h1,{'final','true','start','best'},'Location','southwest')
% legend({'start','final'},'Location','southwest')

subplot(2,2,2); box on; hold on;
cref = dispR_surf96(periods,refmod); % predicted phase velocity
h2(1) = errorbar(periods,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
h2(2) = plot(periods,cref,'-ob','linewidth',2);
h2(3) = plot(periods,cpre(:,imin),'--g','linewidth',2);
h = plot(periods,cpre,'-or','linewidth',1);
uistack(h,'bottom');
h2(4) = h(1);
legend(h2,{'c obs','c start','c best','c ensemble'},'Location','southeast')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

%% Histograms

figure(1001); clf;
for ic = 1:Ncoeffs
    
    subplot(ceil(sqrt(Ncoeffs)),ceil(sqrt(Ncoeffs)),ic);
    plot(bayesian.priors.vs_vec,bayesian.priors.pdf_sp{ic},'-k','linewidth',2); hold on;
    plot(vs_vec,bayesian.post.vs_pdf_sp(ic,:),'-r','linewidth',1.5); hold on;
    ylims = get(gca,'YLim');
%     plot(spcoeffs_true(ic)*[1 1],ylim,'--g','linewidth',1.5);
    title(['Coefficient ',num2str(ic)]);
    xlim([min(vs_edges) max(vs_edges)]);
    
end
    
%% Plot Vs Likelihood PDF

vs_pdf_pl = bayesian.post.vs_pdf;
% % qmu_inv_pdf_pl(log10(qmu_inv_pdf_pl)<-5) = nan;

figure(5); clf;
set(gcf,'Position',[189   507   689   518],'color','w')


ax1 = subplot(2,2,[1 3]); box on; hold on;
surface(bayesian.vsgrid,bayesian.zgrid,log10(vs_pdf_pl)); shading interp;
% surface(bayesian.vsgrid,bayesian.zgrid,zeros(size(vs_pdf_pl)),log10(vs_pdf_pl)); shading interp;
% surface(bayesian.vs_mat_lays-mean(diff(bayesian.vs_vec)),bayesian.z_mat_lays,zeros(size(bayesian.post.vs_pdf_lays)),log10(bayesian.post.vs_pdf_lays),'edgecolor','none');
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
h = plotlayermods(refmod(:,1),refmod(:,3),'-b');
h.LineWidth = 2;
clr = [0 0 1];
plot(bayesian.post.vs_med,bayesian.z_int,'-','color',clr,'linewidth',3); hold on;
plot(bayesian.vs_mat_best,bayesian.z_int,'-r','linewidth',3); hold on;
plot(bayesian.post.vs_mean,bayesian.z_int,'-c','linewidth',3);
plot(bayesian.post.vs_med,bayesian.z_int,'-','color',[1 0.7 0],'linewidth',3);
plot(bayesian.post.vs_l95,bayesian.z_int,'--','color',[1 0.7 0],'linewidth',3);
plot(bayesian.post.vs_u95,bayesian.z_int,'--','color',[1 0.7 0],'linewidth',3);
contour(bayesian.vsgrid,bayesian.zgrid,bayesian.post.vs_conf_mat,[0.95 0.95],'-m');
caxis([-5 0]);
pos = get(gca,'Position');
cb = colorbar;
ylabel(cb,'log_{10}(Probability)');
set(cb,'linewidth',1.5);
set(gca,'Position',pos);
colormap(viridis);
% colormap(flip(cptcmap('GMT_haxby')));
xlabel('Vs (km/s)');
ylabel('Depth (km)');
set(gca,'ydir','reverse','Position',[ax1.Position(1)-0.05 ax1.Position(2) ax1.Position(3) ax1.Position(4)]);
set(gca,'fontsize',15,'linewidth',1.5);
cbpos = get(cb,'position');
set(cb,'position',[cbpos(1)+0.09 cbpos(2) cbpos(3) 0.4]);
title(bayesian.PROJ,'fontweight','bold','fontsize',20);
xlim([3.7 5.5]);

pos = ax1.Position;
ax2 = axes('Position',[pos(1)+0.35 pos(2) pos(3)*0.25 pos(4)]);
% fill(bayesian.post.disc_pdf,bayesian.z_int,'-k','linewidth',2,'FaceColor',[0.9 0 0]);
ylim(ax1.YLim);
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');
yticklabels([]);


% PdF

phv_pdf_pl = bayesian.post.phv_pdf;
% % qinv_pdf_pl(log10(qinv_pdf_pl)<-5) = nan;
% qinv_meanpost = sum(qinv_pdf.*qinv_grid,1);

ax3 = subplot(2,2,2); box on; hold on;
surface(bayesian.periods_grid,bayesian.phv_grid,log10(phv_pdf_pl)); shading interp;
% scatter(bayesian.periods_grid(:),bayesian.qinv_grid(:),80,log10(qinv_pdf_pl(:)),'filled','markeredgecolor','k');
errorbar(bayesian.periods,bayesian.cobs,bayesian.cstd,'-ok','linewidth',2); hold on;
% plot(bayesian.periods,bayesian.phv_med_pre,'-','color',clr,'linewidth',3); hold on;
plot(bayesian.periods,bayesian.phv_pre_mat_best,'-r','linewidth',3);
% plot(periods_v,phv_meanpost,'-c','linewidth',3);
% plot(bayesian.periods,bayesian.post.phv_mean_pre,'-c','linewidth',3);
plot(bayesian.periods,bayesian.post.phv_med_pre,'-','color',[1 0.7 0],'linewidth',3);
plot(bayesian.periods,bayesian.post.phv_l95_pre,'--','color',[1 0.7 0],'linewidth',3);
plot(bayesian.periods,bayesian.post.phv_u95_pre,'--','color',[1 0.7 0],'linewidth',3);
contour(bayesian.periods_grid,bayesian.phv_grid,bayesian.post.phv_conf_mat,[0.95 0.95],'-m');
xlim([min(bayesian.periods) max(bayesian.periods)]);
ylim([2.5 4.4]);
caxis([-5 0]);
% pos = get(gca,'Position');
% cb = colorbar;
% ylabel(cb,'log_{10}(Probability)');
% set(cb,'linewidth',1.5);
% set(gca,'Position',pos);
% colormap(viridis);
ylabel('Phase Velocity (km/s)');
xlabel('Period (s)');
set(gca,'fontsize',15,'linewidth',1.5,'Position',[ax3.Position(1)+0.06 ax3.Position(2) ax3.Position(3) ax3.Position(4)]);

%% Plot PdF of knot locations

% PdF
figure(6); clf;
set(gcf,'color','w');
box on; hold on;

zsp_pdf_pl = bayesian.post.zsp_pdf ;

knot_ind = [1:size(bayesian.post.zsp_pdf,1)];
[z_int_mat,knot_ind_mat] = meshgrid(bayesian.z_int_vec,knot_ind);

% w = sum(posteriorsp,1);
% bayesian.post.qmu_inv_mean = sum(w.*bayesian.qmu_inv_mat,2)./sum(w);

% surface(knot_ind_mat-0.5,z_int_mat,log10(zsp_pdf_pl),'EdgeColor','none'); %shading interp;
imagesc(knot_ind,bayesian.z_int_vec,log10(zsp_pdf_pl')); %shading interp;
plot(knot_ind,zsp,'ok','MarkerFaceColor','w','markersize',15);
plot(knot_ind,bayesian.post.zsp_mean,'+','color',[0 0 0],'markersize',12,'linewidth',5);
plot(knot_ind,bayesian.post.zsp_mean,'+c','color',[1 0 0],'markersize',10,'linewidth',3);
xlim([0.5 max(knot_ind)+0.5]);
ylim([min(bayesian.z_int),max(bayesian.z_int)]);
caxis([-5 0]);
% pos = get(gca,'Position');
cb = colorbar;
ylabel(cb,'log_{10}(Probability)');
set(cb,'linewidth',1.5);
title('Knot Depth');
% set(gca,'Position',pos);
colormap(viridis);
ylabel('Depth (km)');
xlabel('Knot index');
set(gca,'fontsize',15,'linewidth',1.5,'ydir','reverse');


%% Plot 2-D marginal probabilities

figure(1002); clf; box on; hold on;
% marginal_pdf_lays(marginal_pdf_lays==0) = nan;
surface(bayesian.vs_mat_lays-mean(diff(bayesian.vs_vec)),bayesian.z_mat_lays,zeros(size(bayesian.post.vs_pdf_lays)),log10(bayesian.post.vs_pdf_lays),'edgecolor','none');
% surface(bayesian.vs_mat-mean(diff(bayesian.vs_vec)),bayesian.z_mat,zeros(size(bayesian.marginal_pdf_mat)),(bayesian.marginal_pdf_mat),'edgecolor','none');
cb = colorbar;
ylabel(cb,'log_{10}(Probability)')
set(cb,'linewidth',1.5,'fontsize',16);
set(gca,'fontsize',16,'ydir','reverse','linewidth',1.5);
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
h = plotlayermods(refmod(:,1),refmod(:,3),'-b');
h.LineWidth = 2;
w = sum(posterior,1);
vs_med = sum(w.*vs_models,2)./sum(w);
h = plotlayermods(refmod_sp(:,1),vs_med,'--r');
h.LineWidth = 2;
xlabel('Vs (km/s)');
ylabel('Depth (km)');
xlim([min(bayesian.vs_mat(:)) max(bayesian.vs_mat(:))]);
caxis([log10(1e-5) log10(1)]);

