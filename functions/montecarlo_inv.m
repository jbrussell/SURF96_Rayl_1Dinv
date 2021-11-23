function inv_struct = montecarlo_inv(data_str)
% Do montecarlo inversion by perturbing a layerized starting model.
% First iteration inverts the reference model (i.e., no perturbations)
%
% jbrussell 5/28/2020
isplot = 1;
make_par_surf96('R');
setup_parameters;

% parameter setup
test_N = parameters.testN;
% crusth_var = parameters.crusth_var;
% sedh_var = parameters.sedh_var;
velocity_var = parameters.velocity_var;
laythick_var = parameters.laythick_var;
velocity_var_sed = parameters.velocity_var_sed;
laythick_var_sed = parameters.laythick_var_sed;
% r = parameters.r;
par = parameters.par;
error_tol = parameters.error_tol;
maxdepth = parameters.maxdepth;
depth_nodes = parameters.depth_nodes;
vec_nodes = parameters.vec_nodes;

% sedlayernum = parameters.sedlayernum;
% crustlayernum = parameters.crustlayernum;
% mantle_layer_thickness = parameters.mantle_layer_thickness;
% vpvs = parameters.vpvs;

% niteration = parameters.niteration;
% lambda = parameters.lambda;
% ddlay = parameters.ddlay;

% % cut reference model at max depth
% refmod = data_str.refmod;
% ind_cut = data_str.zlays <= maxdepth;
% refmod_cut = refmod(ind_cut,:);

% reference model
refmod = data_str.refmod;
h_ref = refmod(:,1);
vp_ref = refmod(:,2);
vs_ref = refmod(:,3);
rho_ref = refmod(:,4);
discs = data_str.discs;
zlays = data_str.zlays;

if vs_ref(discs(1)==zlays) == 0
    % Water layer
    waterdepth = discs(1);
    seddepth = discs(2);
    mohodepth = discs(end);
else
    % No water layer
    waterdepth = 0;
    seddepth = discs(1);
    mohodepth = discs(end);
end
ind_water = find(zlays<=waterdepth);
ind_sed = find(zlays<=seddepth & zlays>waterdepth);
ind_crust = find(zlays<=mohodepth & zlays>seddepth);
% ind_mantle = find(zlays>mohodepth);
ind_mantle = find(zlays>mohodepth & zlays<=maxdepth);
ind_belowbot = find(zlays>maxdepth);

% obtain the dispersion curve
velT = data_str.periods;
phv = data_str.phv;
phvstd = data_str.phvstd;
grv = data_str.grv;
grvstd = data_str.grvstd;
% topo = data_str.topo;

% run the tests
% for itest = 1:test_N
it = 0;
itest = 1;
while itest <= test_N && it < test_N*2
    it = it + 1;
	clear vec_h vec_vs vec_vp vec_rho initmodel mod_waterdepth mod_seddepth mod_mohodepth
	% reset the initial model
	vec_h = h_ref;
	vec_vp = vp_ref;
    vec_vs = vs_ref;
	vec_rho = rho_ref;
	
    % perturb model velocities within each layer (excluding water column)
    
    % Sediments
    % perturb layer thicknesses
    vec_h(ind_sed) = vec_h(ind_sed) .* (1 + normrnd(0,laythick_var_sed,1));
    % perturb layer velocities
    sed_perturb = (1 + normrnd(0,velocity_var_sed)) * ones(size(ind_sed));
%     sed_perturb = smooth(1 + normrnd(0,velocity_var,size(ind_sed)));
%     vec_vp(ind_sed) = vec_vp(ind_sed) .* sed_perturb;
    vec_vs(ind_sed) = vec_vs(ind_sed) .* sed_perturb;
%     vec_rho(ind_sed) = vec_rho(ind_sed) .* sed_perturb;
    
    % Crust
    % perturb layer thicknesses
    vec_h(ind_crust) = vec_h(ind_crust) .* (1 + normrnd(0,laythick_var,1));
    % perturb layer velocities
    crust_perturb = (1 + normrnd(0,velocity_var)) * ones(size(ind_crust));
%     crust_perturb = smooth(1 + normrnd(0,velocity_var,size(ind_crust)));
%     vec_vp(ind_crust) = vec_vp(ind_crust) .* crust_perturb;
    vec_vs(ind_crust) = vec_vs(ind_crust) .* crust_perturb;
%     vec_rho(ind_crust) = vec_rho(ind_crust) .* crust_perturb;
    
    % Mantle
    % perturb layer velocities
    mantle_perturb = (1 + normrnd(0,velocity_var)) * ones(size(ind_mantle));
%     mantle_perturb = smooth(1 + normrnd(0,velocity_var,size(ind_mantle)));
%     vec_vp(ind_mantle) = vec_vp(ind_mantle) .* mantle_perturb;
    vec_vs(ind_mantle) = vec_vs(ind_mantle) .* mantle_perturb;
%     vec_rho(ind_mantle) = vec_rho(ind_mantle) .* mantle_perturb;
    
    % Ensure that mantle velocities are faster than crust (otherwise surf96 sometimes hangs up)
    fac_v = vec_vs(ind_crust(end))/vec_vs(ind_mantle(1));
    if fac_v > 1
        fac_v = fac_v + 0.05;
%         vec_vp(ind_mantle) = vec_vp(ind_mantle) * fac_v;
        vec_vs(ind_mantle) = vec_vs(ind_mantle) * fac_v;
%         vec_rho(ind_mantle) = vec_rho(ind_mantle) * fac_v;
    end
    
    % Reset to reference model for 1st iteration
    if it == 1
        vec_h = h_ref;
        vec_vp = vp_ref;
        vec_vs = vs_ref;
        vec_rho = rho_ref;
    end
    
    % Get depths to boundaries
    zlays = cumsum(vec_h);
    if ~isempty(ind_water)
        mod_waterdepth = zlays(ind_water(end));
    else
        mod_waterdepth = 0;
    end
    mod_seddepth = zlays(ind_sed(end));
    mod_mohodepth = zlays(ind_crust(end));
    
    % Tie model into ref model below max depth
    maxdepit = zlays(ind_mantle(end));
    dz = 50;
    ind_botmod_lin = find(zlays>maxdepit & zlays<=maxdepit+dz);
    Nlays = length(ind_botmod_lin);
    vec_vp(ind_botmod_lin) = linspace(vec_vp(ind_botmod_lin(1)-1),vec_vp(ind_botmod_lin(end)),Nlays);
    vec_vs(ind_botmod_lin) = linspace(vec_vs(ind_botmod_lin(1)-1),vec_vs(ind_botmod_lin(end)),Nlays);
    vec_rho(ind_botmod_lin) = linspace(vec_rho(ind_botmod_lin(1)-1),vec_rho(ind_botmod_lin(end)),Nlays);
    
    discs_pert = [mod_waterdepth, mod_seddepth, mod_mohodepth];
    
	initmodel(:,1) = vec_h(:);
	initmodel(:,2) = vec_vp(:);
	initmodel(:,3) = vec_vs(:);
	initmodel(:,4) = vec_rho(:);
	err = disperr(velT,phv,phvstd,grv,grvstd,initmodel)
    if isempty(err)
        disp('Issue calculating error. Skipping...');
        continue
    end
	if err>500 %err > 0.3
% 		errors(itest) = 1000;
		continue
    end
    
    if isplot
        figure(1); clf;
        subplot(2,2,[1 3]); hold on;
        fsize=20;
        h = plotlayermods(h_ref,vs_ref);
        set(h,'linewidth',2,'color',[0 0 0]);
        h = plotlayermods(initmodel(:,1),initmodel(:,3));
        set(h,'linewidth',2,'color',[1 0 0]);
        ylim([-250 0])
        % xlim([2 5])
        legend('reference','starting','location','southwest');
        ylabel('Depth (km)','fontsize',fsize);
        xlabel('Shear Velocity (km/s)','fontsize',fsize);
        set(gca,'fontsize',fsize)
        drawnow;
    end
    
    % Run the inversion
	try
%     [outmod, phv_fwd] = invdispR(velT,phv,phvstd,grv,grvstd,initmodel,mod_mohodepth,mod_waterdepth,mod_seddepth,niteration,lambda,ddlay);
        [outmod, phv_fwd] = invdispR_lsqr(velT,phv,phvstd,grv,grvstd,initmodel,discs_pert,par,maxdepit);
	catch
		disp('Error in inversion');
		outmod = [];
		phv_fwd = [];
    end
	if ~isempty(outmod)
        misfit = (phv_fwd(:)-phv(:))./phvstd;
        rms=sqrt(sum(misfit.^2)/length(velT));
        chi2_red = sum(misfit.^2)/length(velT)
        final_mods(:,itest) = outmod(:,3);
        init_mods(:,itest) = initmodel(:,3);
        vec_hs(:,itest) = outmod(:,1);
        phv_fwds(:,itest) = phv_fwd(:);
    % 	errors(itest) = rms;
        errors(itest) = chi2_red;
        discs_mod(:,itest) = discs_pert;
	else
% 		errors(itest) = 999;
        continue
    end
    if isplot && ~isempty(outmod)
        figure(1);
        subplot(2,2,[1 3]);
        h = plotlayermods(outmod(:,1),outmod(:,3));
        set(h,'linewidth',2,'color',[0 0 1]);
        drawnow;
        
        phv_init = dispR_surf96(velT,initmodel);
        
        subplot(2,2,2); hold on;
        errorbar(velT,phv,phvstd,'ok');
        plot(velT,phv_init,'sr');
        plot(velT,phv_fwd,'sb');
        ylabel('Phase Velocity (km/s)','fontsize',fsize);
        xlabel('Period (s)','fontsize',fsize);
        set(gca,'fontsize',fsize)
        drawnow;
%         pause;
    end
    itest = itest + 1;
end
if size(init_mods,2) < length(errors)
	final_mods(:,test_N) = 0;
	init_mods(:,test_N) = 0;
	vec_hs(:,test_N) = 0;
	phv_fwds(:,test_N) = 0;
    discs_mod(:,test_N) = 0;
end

% % find bad inversions
% % error_tol = median(errors(find(errors<0.3)))*(1+r);
% % error_tol = median(errors(find(errors<1.5)))*(1+r);
% % error_tol = 1.5;
% bad_test_ind = find(errors > error_tol);
% errors(bad_test_ind) = [];
% init_mods(:,bad_test_ind) = [];
% final_mods(:,bad_test_ind) = [];
% vec_hs(:,bad_test_ind) = [];
% phv_fwds(:,bad_test_ind) = [];
% discs_mod(:,bad_test_ind) = [];
% test_N = test_N - length(bad_test_ind);

[depth_node init_avg init_std] = mod_statistic(vec_hs,init_mods,depth_nodes,vec_nodes);
[depth_node final_avg final_std] = mod_statistic(vec_hs,final_mods,depth_nodes,vec_nodes);

inv_struct.final_mods = final_mods;
inv_struct.init_mods = init_mods;
inv_struct.vec_hs = vec_hs;
inv_struct.phv_fwds = phv_fwds;
inv_struct.errors = errors;
inv_struct.init_avg = init_avg;
inv_struct.final_avg = final_avg;
inv_struct.init_std = init_std;
inv_struct.final_std = final_std;
inv_struct.depth_node = depth_node;
inv_struct.discs_mod = discs_mod;
