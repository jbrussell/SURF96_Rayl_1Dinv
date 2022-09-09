% Test simple Markov chain Monte Carlo Bayesian inversion for 3-layer
% velocity model. The method mostly follows Shen et al. (2013) GJI doi:10.1093/gji/ggs050
% This version uses smooth splines rather than layers.
%
% jbrussell 9/7/2022
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

% % Compile the faster mex files for spline calculation
% % !!!!! This only needs to be compiled the first time !!!!!
% cd ./functions
% CompileMexFiles
% cd ..

%% MCMC parameters
% Other inversion parameters
nit_mcmc = 2000; % total number of iterations
nit_restart = 250; %1e10; % number of iterations after which to restart with new random model (if never want to restart, set to giant number)
N_cooldown = 100; %50; % number of iterations over which temperature parameter (tau) decays
m_perturb_method = 'all'; % 'single' (perturb one model parameter at a time) | 'all' (perturb all at once)
nit_plot = 250; % number of iterations after which to plot

% Define bounds of allowed model space M relative to ref. model. For the spline
% inversion, this applies to the Vs spline coefficients, not the layers.
% Models occuring outside this space will not be allowed.
% (these values also act as the min and max of the uniform prior)
% If a water layer exists, it is held at fixed velocity/density
par.dv_M = [-0.25 +0.25]; % pct of reference model

% Define widths of gaussian perturbations made at each iteration
par.dv_std = 0.05; % km/s

% Scale vp and density with vs
par.vp_vs = 1.75; % Vp/Vs
par.rho_vs = 0.74; % density/Vs

% Spline parameters
Nspline = 6; % Number of desired splines, evenly spaced from surface (or base of water layer) to zmax
dz_int = 5; % (km) interpolated layer thicknesses. If too small, surf96 will break...
zmax = 200; % Maximum depth of starting model

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
z_true = [0 zh2o_true zsed_true zmoho_true zmoho_true+40 zmoho_true+100 zmax];
dz = [diff(z_true) 0];
vs = [0 2.5 3.5 4.7 4.2 4.4 4.4]; 
% vp = 1.75*vs; vp(1)=1.5;
vp = par.vp_vs * vs; vp(vs==0)=1.5;
% rho = vp / 2.5; rho(1)=1.03;
rho = par.rho_vs * vs; rho(vs==0)=1.03;
truemod = [dz(:), vp(:), vs(:), rho(:)];

% GENERATE SYNTHETIC DATASET
% Calculate dispersion for true model, which will be our "observations"
periods = logspace(log10(10),log10(40),10);
cobs = dispR_surf96(periods,truemod); % "observations"
cstd = cobs * 0.01; % observation uncertainties

%% Starting model

% 3-layer model
zh2o = 1.618; % [km] water depth
zsed = 7; 
zmoho = 20; % [km] moho depth
% zlab = 70;
z = [0 zh2o zsed zmoho zmax];
dz = [diff(z) 0];
vs = [0 2.5-0.2 3.5-0.4 4.4+0.4 4.4]; 
% vs = [0 3.5 3.5 3.5 4.4]; 
% vp = 1.75*vs; vp(1)=1.5;
% vp = truemod(:,2);
vp = par.vp_vs * vs; vp(vs==0)=1.5;
% rho = vp / 2.5; rho(1)=1.03;
% rho = truemod(:,4);
rho = par.rho_vs * vs; rho(vs==0)=1.03;
refmod = [dz(:), vp(:), vs(:), rho(:)];

%%
% Interpolate layered models
zinterp = [0 zh2o zh2o:dz_int:zmax]';
% zsp = z(:);
% zsp = [0:10:zmoho zmoho+50:50:zmax]; 
[mod_true] = layerizemod_interp(truemod,zinterp);
[mod_ref] = layerizemod_interp(refmod,zinterp);

% Do spline calculations for reference (and true) model
% zsp = [zh2o:50:zmax];
% zsp = [linspace(zh2o,zmoho,5-1) linspace(zmoho,zmax,6-1)]; % example of custom spline spacing with a discontinuity
zsp = linspace(zh2o,zmax,Nspline-1);
Inoh2o = find(mod_true.vs~=0);
Ih2o = find(mod_true.vs==0);
[spbasis_true,spcoeffs_true,spzz_true]=make_splines(zsp(:),[],mod_true.z(Inoh2o),mod_true.vs(Inoh2o));
vs_true_sp = spbasis_true * spcoeffs_true;
vs_true_sp = [mod_true.vs(Ih2o); vs_true_sp];
Inoh2o = find(mod_ref.vs~=0);
Ih2o = find(mod_ref.vs==0);
[spbasis,spcoeffs,spzz]=make_splines(zsp(:),[],mod_ref.z(Inoh2o),mod_ref.vs(Inoh2o));
vs_ref_sp = spbasis * spcoeffs;
vs_ref_sp = [mod_ref.vs(Ih2o); vs_ref_sp];

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
xlabel('Velocity');
ylabel('Depth');
set(gca,'FontSize',18,'linewidth',1.5);
legend({'true','true (spline)','start','start (spline)'},'Location','southwest')
% legend({'start','final'},'Location','southwest')

subplot(2,2,2); box on; hold on;
cref = dispR_surf96(periods,refmod); % predicted phase velocity
[truemod_sp] = spline2mod(mod_true,vs_true_sp,par.vp_vs,par.rho_vs);
[refmod_sp] = spline2mod(mod_ref,vs_ref_sp,par.vp_vs,par.rho_vs);
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

% Uniform priors spanning M
priors.sample = @(N,ic) unifrnd(model_bounds(ic,1), model_bounds(ic,2) ,N,1);

% Function to perturb model
perturb_model = @(model,std_vec) normrnd(model(:)',std_vec)';

% Get pdf from distributions
vs_edges = [0:0.04:7];
vs_vec = 0.5*(vs_edges(1:end-1)+vs_edges(2:end));
figure(1000);
for ic = 1:Ncoeffs
    h = histogram(priors.sample(1000000,ic),vs_edges,'Normalization','probability');
    priors.pdf_sp{ic} = h.Values;
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
posterior = nan(size(refmod_sp,1),nit_mcmc);
posterior_sp = nan(Ncoeffs,nit_mcmc);
cpre = nan(length(cobs),nit_mcmc);
misfit = nan(1,nit_mcmc);
Likelihood = nan(1,nit_mcmc);
vs_models = nan(size(refmod_sp,1),nit_mcmc);
vs_models_sp = nan(Ncoeffs,nit_mcmc);
models = nan([size(refmod_sp),nit_mcmc]);

% Initiate
m_j = spcoeffs;
m_j(:,3) = sample_model(priors.sample,1,Ncoeffs);
m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
ii = 0;
ibad = 0;
tic
while ii < nit_mcmc
    
    if ii>0 && mod(ii,nit_restart) == 0 % reinitialize mcmc, start over
        m_j(:,3) = sample_model(priors.sample,1,Ncoeffs);
        m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
        m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
%         ibad = 0;
    end
    
    % Previous model
    Inoh2o = find(mod_ref.vs~=0);
    Ih2o = find(mod_ref.vs==0);
    vs_spline = spbasis * m_j(:,3);
    vs_spline = [mod_ref.vs(Ih2o); vs_spline];
    [splinemod_j] = spline2mod(mod_ref,vs_spline,par.vp_vs,par.rho_vs);
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
        display(['Searching for stable starting model: ',num2str(ibad)]);
        continue
    end
    ii = ii + 1;
    
    if mod(ii,100) == 0
        display([num2str(ii),'/',num2str(nit_mcmc)]);
    end
    
    % Calculate posterior probability of model j (spline coefficients)
    for ic = 1:Ncoeffs
        [~,I] = min(abs(m_j(ic,3)-priors.vs_vec));
        posterior_sp(ic,ii) = L_j .* priors.pdf_sp{ic}(I);
    end
    % Calculate posterior for layered structure
    ipdf = 0;
    for ilay = 1:size(posterior,1)
        if splinemod_j(ilay,3)==0 % water layer
            posterior(ilay,ii) = L_j * 1;
            continue
        end
        ipdf = ipdf + 1;
        [~,I] = min(abs(splinemod_j(ilay,3)-priors.vs_vec));
        posterior(ilay,ii) = L_j .* priors.pdf{ipdf}(I);
    end
%     posterior(:,ii) = L_j;
    
    % Save outputs
    misfit(ii) = S_j;
    Likelihood(ii) = L_j;
    cpre(:,ii) = c_j(:);
    vs_models(:,ii) = splinemod_j(:,3);
    vs_models_sp(:,ii) = m_j(:,3);
    models(:,:,ii) = splinemod_j;
    
    % Decaying thermal parameter (cool down parameter) from simulated 
    % annealing (Kirkpatrick et al. 1983). This allows larger changes 
    % between sequential models at early iterations. This premultiplies the
    % Gaussian distributions from which random model parameters are drawn
    % and also the likelihood of the trial model, so misfit increases are
    % more likely accepted early in the MCMC.
%     tau = 1 + 3 * erfc(ii/500); % denom = 500 means decays over ~1500 iterations (Eilon et al. 2018)
    tau = 1 + 3 * erfc(ii/(N_cooldown/3)); % denom = 500 means decays over ~1500 iterations
    
    % Trial model
    is_in_bounds = 0;
    while is_in_bounds == 0
        m_i = m_j;
        dvs = perturb_model(m_i(:,3),tau*repmat(par.dv_std,1,Ncoeffs)); % perturb Vs
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
        is_in_bounds = is_model_in_bounds(m_i,model_bounds);
    end
    Inoh2o = find(mod_ref.vs~=0);
    Ih2o = find(mod_ref.vs==0);
    vs_spline = spbasis * m_i(:,3);
    vs_spline = [mod_ref.vs(Ih2o); vs_spline];
    [splinemod_i] = spline2mod(mod_ref,vs_spline,par.vp_vs,par.rho_vs);
    c_i = dispR_surf96(periods,splinemod_i); % predicted phase velocity
    if length(c_i) ~= length(periods) % check if something is wrong...
        m_i = m_j; % revert back to previous model
        Inoh2o = find(mod_ref.vs~=0);
        Ih2o = find(mod_ref.vs==0);
        vs_spline = spbasis * m_i(:,3);
        vs_spline = [mod_ref.vs(Ih2o); vs_spline];
        [splinemod_i] = spline2mod(mod_ref,vs_spline,par.vp_vs,par.rho_vs);
        c_i = dispR_surf96(periods,splinemod_i); % predicted phase velocity
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
        plot(1:ii,misfit(1:ii) / length(periods),'o'); hold on;
        ylabel('Misfit');        
        yyaxis right
        plot(1:ii,log10(Likelihood(1:ii)),'o'); hold on;
        ylabel('log_{10}(Likelihood)');
        
        subplot(2,2,[2 4]); box on; hold on;
        for kk = 1:ii
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
        plot(periods,cpre(:,1:ii),'-or','linewidth',1);
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
    else
        % Reject new model i
        continue
    end
end
toc

%% Calculate marginal pdfs
vs_edges = [0:0.04:7];
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

% 2-D MARGINAL PDFs
% Convert from layers defined by a center point and width to knots like mineos
z = [0; cumsum(refmod_sp(1:end-1,1))];
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
bayesian.vs_mat = VS;
bayesian.z_mat = Z_lays;
bayesian.marginal_pdf_mat = marginal_pdf_lays;
bayesian.marginal_pdf_vec = marginal_pdf;
bayesian.marginal_pdf_vec_sp = marginal_pdf_sp;
bayesian.vs_vec = vs_vec;

%% Plot misfit/likelihood evolution
figure(888); clf;
subplot(2,1,1);
plot(1:nit_mcmc,misfit / length(periods),'o'); hold on;
xlabel('Model #');
ylabel('Misfit');

subplot(2,1,2);
plot(1:nit_mcmc,log10(Likelihood),'o'); hold on;
xlabel('Model #');
ylabel('log_{10}(Likelihood)');

%% Plot Vs models
figure(100); clf; 
set(gcf,'position',[370   372   967   580]);

subplot(2,2,[1 3]); box on; hold on;
for ii = 1:nit_mcmc
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
    
    subplot(3,3,ic);
    plot(priors.vs_vec,priors.pdf_sp{ic},'-k','linewidth',2); hold on;
    plot(vs_vec,bayesian.marginal_pdf_vec_sp(ic,:),'-r','linewidth',1.5); hold on;
    ylim = get(gca,'YLim');
    plot(spcoeffs_true(ic)*[1 1],ylim,'--g','linewidth',1.5);
    title(['Coefficient ',num2str(ic)]);
    xlim([min(vs_edges) max(vs_edges)]);
    
end
    
%% Plot 2-D marginal probabilities

figure(1002); clf; box on; hold on;
% marginal_pdf_lays(marginal_pdf_lays==0) = nan;
surface(bayesian.vs_mat-mean(diff(bayesian.vs_vec)),bayesian.z_mat,zeros(size(bayesian.marginal_pdf_mat)),log10(bayesian.marginal_pdf_mat),'edgecolor','none');
% surface(bayesian.vs_mat-mean(diff(bayesian.vs_vec)),bayesian.z_mat,zeros(size(bayesian.marginal_pdf_mat)),(bayesian.marginal_pdf_mat),'edgecolor','none');
cb = colorbar;
ylabel(cb,'log_{10}(Probability)')
set(cb,'linewidth',1.5,'fontsize',16);
set(gca,'fontsize',16,'ydir','reverse','linewidth',1.5);
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
w = sum(posterior,1);
vs_med = sum(w.*vs_models,2)./sum(w);
h = plotlayermods(refmod_sp(:,1),vs_med,'--r');
h.LineWidth = 2;
xlabel('Vs (km/s)');
ylabel('Depth (km)');
xlim([min(bayesian.vs_mat(:)) max(bayesian.vs_mat(:))]);
caxis([log10(1e-5) log10(1)]);

