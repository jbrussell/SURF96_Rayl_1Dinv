% Test simple Markov chain Monte Carlo Bayesian inversion for 3-layer 
% velocity model. The method mostly follows Shen et al. (2013) GJI doi:10.1093/gji/ggs050
%
% jbrussell 8/3/2022
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

%% MCMC parameters
% Other inversion parameters
nit_mcmc = 1000; % total number of iterations
nit_restart = 250; %1e10; % number of iterations after which to restart with new random model (if never want to restart, set to giant number)
N_cooldown = 100; %50; % number of iterations over which temperature parameter (tau) decays
m_perturb_method = 'all'; % 'single' (perturb one model parameter at a time) | 'all' (perturb all at once)
nit_plot = 250; % number of iterations after which to plot

% Define bounds of allowed model space M relative to ref. model.
% Models outside this space will not be allowed 
% (this also acts as the min and max of the uniform prior)
par.lay1.dv_M = [0 0]; % pct of reference model  (water) keep fixed
par.lay2.dv_M = [-0.3 +0.3]; % pct of reference model (sediment)
par.lay3.dv_M = [-0.3 +0.3]; % pct of reference model  (crust)
par.lay4.dv_M = [-0.3 +0.3]; % pct of reference model  (mantle)
par.lay5.dv_M = [0 0]; % pct of reference model  (halfspace) keep fixed

% Define widths of gaussian perturbations made at each iteration
par.lay1.dv_std = 0; % km/s (water) keep fixed
par.lay2.dv_std = 0.02; % km/s (sediment)
par.lay3.dv_std = 0.05; % km/s (crust)
par.lay4.dv_std = 0.05; % km/s (mantle)
par.lay5.dv_std = 0; % km/s (halfspace) keep fixed

% Scale vp and density with vs
par.vp_vs = 1.75; % Vp/Vs
par.rho_vs = 0.74; % density/Vs

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
zh2o = 1.618; % [km] water depth
zsed = 7; 
zmoho = 20; % [km] moho depth
% zlab = 70;
zmax = 200;
z = [0 zh2o zsed zmoho zmax];
dz = [diff(z) 0];
vs = [0 2.5 3.5 4.4 4.4]; 
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
zmax = 200;
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

figure(1); clf;
set(gcf,'position',[370   372   967   580]);
subplot(2,2,[1 3]); box on; hold on;
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
h = plotlayermods(refmod(:,1),refmod(:,3),'-b');
h.LineWidth = 2;
xlabel('Velocity');
ylabel('Depth');
set(gca,'FontSize',18,'linewidth',1.5);
legend({'true','start','final'},'Location','southwest')
% legend({'start','final'},'Location','southwest')

subplot(2,2,2); box on; hold on;
cref = dispR_surf96(periods,refmod); % predicted phase velocity
errorbar(periods,cobs,2*cstd,'sk','markersize',8,'markerfacecolor','k','linewidth',2);
plot(periods,cref,'-ob','linewidth',2);
legend({'c obs','c start'},'Location','southeast')
xlabel('Period');
ylabel('Phase Velocity');
set(gca,'FontSize',18,'linewidth',1.5);

%% Define priors for each layer
% Define edges of the model space M
model_bounds = [
                refmod(1,3)*(1+par.lay1.dv_M(1)) , refmod(1,3)*(1+par.lay1.dv_M(2));
                refmod(2,3)*(1+par.lay2.dv_M(1)) , refmod(2,3)*(1+par.lay2.dv_M(2));
                refmod(3,3)*(1+par.lay3.dv_M(1)) , refmod(3,3)*(1+par.lay3.dv_M(2));
                refmod(4,3)*(1+par.lay4.dv_M(1)) , refmod(4,3)*(1+par.lay4.dv_M(2));
                refmod(5,3)*(1+par.lay5.dv_M(1)) , refmod(5,3)*(1+par.lay5.dv_M(2));
               ];

% Uniform priors spanning M
priors.lay1 = @(N) unifrnd(model_bounds(1,1), model_bounds(1,2) ,N,1);
priors.lay2 = @(N) unifrnd(model_bounds(2,1), model_bounds(2,2) ,N,1);
priors.lay3 = @(N) unifrnd(model_bounds(3,1), model_bounds(3,2) ,N,1);
priors.lay4 = @(N) unifrnd(model_bounds(4,1), model_bounds(4,2) ,N,1);
priors.lay5 = @(N) unifrnd(model_bounds(5,1), model_bounds(5,2) ,N,1);

% Function to draw model from prior
sample_model = @(N) [priors.lay1(N) priors.lay2(N) priors.lay3(N) priors.lay4(N) priors.lay5(N)]';

% Function to perturb model
perturb_model = @(model,std_vec) normrnd(model(:)',std_vec)';

% Get pdf from distributions
vs_edges = [0:0.04:5];
vs_vec = 0.5*(vs_edges(1:end-1)+vs_edges(2:end));
figure(1000);
h = histogram(priors.lay1(1000000),vs_edges,'Normalization','probability');
priors.lay1_pdf = h.Values;
h = histogram(priors.lay2(1000000),vs_edges,'Normalization','probability');
priors.lay2_pdf = h.Values;
h = histogram(priors.lay3(1000000),vs_edges,'Normalization','probability');
priors.lay3_pdf = h.Values;
h = histogram(priors.lay4(1000000),vs_edges,'Normalization','probability');
priors.lay4_pdf = h.Values;
h = histogram(priors.lay5(1000000),vs_edges,'Normalization','probability');
priors.lay5_pdf = h.Values;
priors.vs_vec = vs_vec;

figure(999); clf;
plot(priors.vs_vec,priors.lay1_pdf); hold on;
plot(priors.vs_vec,priors.lay2_pdf);
plot(priors.vs_vec,priors.lay3_pdf);
plot(priors.vs_vec,priors.lay4_pdf);
plot(priors.vs_vec,priors.lay5_pdf);
title('Priors');

%% Do MCMC
posterior = nan(size(refmod,1),nit_mcmc);
cpre = nan(length(cobs),nit_mcmc);
misfit = nan(1,nit_mcmc);
Likelihood = nan(1,nit_mcmc);
vs_models = nan(size(refmod,1),nit_mcmc);
models = nan([size(refmod),nit_mcmc]);

% Initiate
m_j = refmod;
m_j(:,3) = sample_model(1);
m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
ii = 0;
ibad = 0;
ii_cooldown = 0;
tic
while ii < nit_mcmc
    
    if ii>0 && mod(ii,nit_restart) == 0 % reinitialize mcmc, start over
        m_j(:,3) = sample_model(1);
        m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
        m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
%         ibad = 0;
        ii_cooldown = 0;
    end
    
    % Previous model
    c_j = dispR_surf96(periods,m_j); % predicted phase velocity
    if length(c_j) ~= length(periods) % check if something is wrong...
        ibad = ibad+1;
        m_j(:,3) = sample_model(1);
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
        m_j(:,3) = sample_model(1);
        m_j(:,2) = par.vp_vs*m_j(:,3); m_j(m_j(:,3)==0,2)=1.5;
        m_j(:,4) = par.rho_vs*m_j(:,3); m_j(m_j(:,3)==0,4)=1.03;
        display(['Searching for stable starting model: ',num2str(ibad)]);
        continue
    end
    ii = ii + 1;
    
    if mod(ii,100) == 0
        display([num2str(ii),'/',num2str(nit_mcmc)]);
    end
    
    % Calculate posterior probability of model j
    [~,I] = min(abs(m_j(1,3)-priors.vs_vec));
    posterior(1,ii) = L_j .* priors.lay1_pdf(I);
    [~,I] = min(abs(m_j(2,3)-priors.vs_vec));
    posterior(2,ii) = L_j .* priors.lay2_pdf(I);
    [~,I] = min(abs(m_j(3,3)-priors.vs_vec));
    posterior(3,ii) = L_j .* priors.lay3_pdf(I);
    [~,I] = min(abs(m_j(4,3)-priors.vs_vec));
    posterior(4,ii) = L_j .* priors.lay4_pdf(I);
    [~,I] = min(abs(m_j(5,3)-priors.vs_vec));
    posterior(5,ii) = L_j .* priors.lay5_pdf(I);
    
    % Save outputs
    misfit(ii) = S_j;
    Likelihood(ii) = L_j;
    cpre(:,ii) = c_j(:);
    vs_models(:,ii) = m_j(:,3);
    models(:,:,ii) = m_j;
    
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
    while is_in_bounds == 0
        m_i = m_j;
        dvs = perturb_model(m_i(:,3),tau*[par.lay1.dv_std par.lay2.dv_std par.lay3.dv_std par.lay4.dv_std par.lay5.dv_std]); % perturb Vs
    %     dvs = sample_model(1); % random Vs
        switch m_perturb_method
            case 'single'
                I_pert = ceil(rand(1)*size(refmod,1)); % randomly pick model parameter to perturb
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
    c_i = dispR_surf96(periods,m_i); % predicted phase velocity
    if length(c_i) ~= length(periods) % check if something is wrong...
        m_i = m_j; % revert back to previous model
        c_i = dispR_surf96(periods,m_i); % predicted phase velocity
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
vs_edges = [0:0.04:5];
vs_vec = 0.5*(vs_edges(1:end-1)+vs_edges(2:end));
marginal_pdf = zeros(size(refmod,1),length(vs_vec));
for idim = 1:size(refmod,1)
    ind_bin = discretize(vs_models(idim,:),vs_edges);
    marginal = zeros(size(vs_vec));
    for ii = 1:length(ind_bin)
        marginal(ind_bin(ii)) = marginal(ind_bin(ii)) + sum(posterior(:,ii));
    end
    marginal_pdf(idim,:) = marginal / sum(marginal); % normalize so sums to 1
end

% 2-D MARGINAL PDFs
% Convert from layers defined by a center point and width to knots like mineos
z = [0; cumsum(refmod(1:end-1,1))];
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
subplot(2,3,1);
plot(priors.vs_vec,priors.lay1_pdf,'-k','linewidth',2); hold on;
plot(vs_vec,bayesian.marginal_pdf_vec(1,:),'-r','linewidth',1.5); hold on;
ylim = get(gca,'YLim');
plot(truemod(1,3)*[1 1],ylim,'--g','linewidth',1.5);
title('Layer 1 (water)');
xlim([0 1]);

subplot(2,3,2);
plot(priors.vs_vec,priors.lay2_pdf,'-k','linewidth',2); hold on;
plot(bayesian.vs_vec,bayesian.marginal_pdf_vec(2,:),'-r','linewidth',1.5); hold on;
% histogram('BinEdges',vs_edges,'BinCounts',marginal_pdf(2,:),'EdgeColor','r','FaceColor','r','FaceAlpha',1); hold on;
% histogram('BinEdges',vs_edges,'BinCounts',priors.lay2_pdf,'EdgeColor','k','FaceColor','none','linewidth',2); hold on;
ylim = get(gca,'YLim');
plot(truemod(2,3)*[1 1],ylim,'--g','linewidth',1.5);
title('Layer 2 (sediments)');
xlim([1 4]);

subplot(2,3,3);
plot(priors.vs_vec,priors.lay3_pdf,'-k','linewidth',2); hold on;
plot(bayesian.vs_vec,bayesian.marginal_pdf_vec(3,:),'-r','linewidth',1.5); hold on;
ylim = get(gca,'YLim');
plot(truemod(3,3)*[1 1],ylim,'--g','linewidth',1.5);
title('Layer 3 (crust)');
xlim([2 5]);

subplot(2,3,4);
plot(priors.vs_vec,priors.lay4_pdf,'-k','linewidth',2); hold on;
plot(bayesian.vs_vec,bayesian.marginal_pdf_vec(4,:),'-r','linewidth',1.5); hold on;
ylim = get(gca,'YLim');
plot(truemod(4,3)*[1 1],ylim,'--g','linewidth',1.5);
title('Layer 4 (mantle)');
xlim([4 5]);

subplot(2,3,5);
plot(priors.vs_vec,priors.lay5_pdf,'-k','linewidth',2); hold on;
plot(bayesian.vs_vec,bayesian.marginal_pdf_vec(5,:),'-r','linewidth',1.5); hold on;
ylim = get(gca,'YLim');
plot(truemod(5,3)*[1 1],ylim,'--g','linewidth',1.5);
title('Layer 5 (halfspace)');
xlim([4 5]);
    
%% Plot 2-D marginal probabilities

% temp = layerizemod(refmod);
% Z_lays = repmat(temp.z,1,size(marginal_pdf,2));
% VS = repmat(vs_vec,size(Z_lays,1),1);
% marginal_pdf_lays = [];
% for ii = 1:size(marginal_pdf,2)
%     temp = layerizemod([marginal_pdf(:,ii) marginal_pdf(:,ii) marginal_pdf(:,ii) marginal_pdf(:,ii)]);
%     marginal_pdf_lays(:,ii) = temp.z;
% end

figure(1002); clf; box on; hold on;
% marginal_pdf_lays(marginal_pdf_lays==0) = nan;
surface(bayesian.vs_mat-mean(diff(bayesian.vs_vec)),bayesian.z_mat,zeros(size(bayesian.marginal_pdf_mat)),log10(bayesian.marginal_pdf_mat),'edgecolor','none');
cb = colorbar;
ylabel(cb,'log_{10}(Probability)')
set(cb,'linewidth',1.5,'fontsize',16);
set(gca,'fontsize',16,'ydir','reverse','linewidth',1.5);
h = plotlayermods(truemod(:,1),truemod(:,3),'-k');
h.LineWidth = 2;
w = sum(posterior,1);
vs_med = sum(w.*vs_models,2)./sum(w);
h = plotlayermods(refmod(:,1),vs_med,'--r');
h.LineWidth = 2;
xlabel('Vs (km/s)');
ylabel('Depth (km)');
xlim([min(bayesian.vs_mat(:)) max(bayesian.vs_mat(:))]);
% caxis([log10(0.3) log10(1)]);

%% Plot tradeoffs
figure(1003); clf;
set(gcf,'position',[322          31        1022         994]);

% 1, 2
subplot(4,4,1); box on; axis square; hold on;
scatter(models(1,3,:),models(2,3,:),50,Likelihood,'.');
ylabel('Sed Vs (lay 2)')
xlabel('Water Vs (lay 1)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 1, 3
subplot(4,4,5); box on; axis square; hold on;
scatter(models(1,3,:),models(3,3,:),50,Likelihood,'.');
ylabel('Crust Vs (lay 3)')
xlabel('Water Vs (lay 1)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 2, 3
subplot(4,4,6); box on; axis square; hold on;
scatter(models(2,3,:),models(3,3,:),50,Likelihood,'.');
ylabel('Crust Vs (lay 3)');
xlabel('Sed Vs (lay 2)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 1,4
subplot(4,4,9); box on; axis square; hold on;
scatter(models(1,3,:),models(4,3,:),50,Likelihood,'.');
ylabel('Mantle Vs (lay 4)')
xlabel('Water Vs (lay 1)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 2,4
subplot(4,4,10); box on; axis square; hold on;
scatter(models(2,3,:),models(4,3,:),50,Likelihood,'.');
ylabel('Mantle Vs (lay 4)');
xlabel('Sed Vs (lay 2)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 3,4
subplot(4,4,11); box on; axis square; hold on;
scatter(models(3,3,:),models(4,3,:),50,Likelihood,'.');
ylabel('Mantle Vs (lay 4)');
xlabel('Crust Vs (lay 3)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 1,5
subplot(4,4,13); box on; axis square; hold on;
scatter(models(1,3,:),models(5,3,:),50,Likelihood,'.');
ylabel('Halfspace Vs (lay 5)')
xlabel('Water Vs (lay 1)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 2,5
subplot(4,4,14); box on; axis square; hold on;
scatter(models(2,3,:),models(5,3,:),50,Likelihood,'.');
ylabel('Halfspace Vs (lay 5)');
xlabel('Sed Vs (lay 2)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 3,5
subplot(4,4,15); box on; axis square; hold on;
scatter(models(3,3,:),models(5,3,:),50,Likelihood,'.');
ylabel('Halfspace Vs (lay 5)');
xlabel('Crust Vs (lay 3)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);

% 4,5
subplot(4,4,16); box on; axis square; hold on;
scatter(models(4,3,:),models(5,3,:),50,Likelihood,'.');
ylabel('Halfspace Vs (lay 5)');
xlabel('Mantle Vs (lay 4)');
set(gca,'fontsize',13,'linewidth',1.5);
colormap(gca,parula);
pos = get(gca,'Position');
cb = colorbar;
set(gca,'Position',pos);
ylabel(cb,'Likelihood');
set(cb,'Ticks',cb.Limits,'TickLabels',{'Low','High'},'linewidth',1.5,'fontsize',13);

% Covariance matrix
subplot(4,4,[3 4 7 8]);
cov_vs = cov(squeeze(models(:,3,:))');
corr_vs = corrcov(cov_vs);
N = size(cov_vs,1);
% imagesc([1:N],[1:N],cov_vs); hold on;
imagesc([1:N],[1:N],corr_vs); hold on;
set(gca,'xtick',[1:N+1],'ytick',[1:N+1])
for ii = 0:N
    plot((ii+0.5)*[1 1],[0.5 N+0.5],'-k','linewidth',1.5);
    plot([0.5 N+0.5],(ii+0.5)*[1 1],'-k','linewidth',1.5);
end
colormap(gca,redblue);
cb = colorbar;
ylabel(cb,'Correlation');
set(cb,'linewidth',1.5);
caxis([-1 1]);
set(gca,'fontsize',16,'linewidth',1.5);
% xlim([0.5 5.5]);
% ylim([0.5 5.5]);

