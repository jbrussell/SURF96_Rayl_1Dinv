function [finalmod,cpre,vs_std,vp_std,rho_std] = run_surf96_inv_Rayl_Vs_Vp_Rho(cobs,cstd,periods,startmod,discs,eps_data,eps_H,eps_J,eps_F,eps_vpvs,eps_rhovs,z_dampbot,nit,nit_recalc_c,vp_vs,rho_vs,nmode)
% Do linearized inversion of Rayleigh wave phase velocities for Vs, Vp, 
% and Density (rho) using surf96 to generate the kernels and to calculate 
% phase velocity.
%
% Both Vp/Vs and Rho/Vs are enforced through prior constraints on 
% Vp and Rho, respectively. Thus, although we technically do solve for
% all three, we solve primarily for Vs and simply scale the other two
% accordingly. Smoothing constraints on Vs will therefore apply also to 
% Vp and Rho. The inversion includes sensitivity kernels for all 3
% parameters, so that perturbations in Vs account for  perturbations in 
% Vp and Rho.
%
% INPUTS (N=number of data; M=number of layers in model)
% cobs - phase velocity observed km/s [N x 1]
% cstd - phase velocity uncertainty km/s [N x 1]
% periods - seconds [N x 1]
% startmod - starting model: dz (km), vp (km/s), vs (km/s), rho (kg/m^3) [M x 4]
% discs - depth of desired discontinuities in smoothing (km) [any length x 1]
% eps_data - weight for data fit as fraction of G norm [scalar]
% eps_H - weight for norm damping as fraction of G norm [scalar]
% eps_J - weight for 1st derivative smoothing as fraction of G norm [scalar]
% eps_F - weight for 2nd derivative smoothing as fraction of G norm [scalar]
% eps_vpvs - weight for enforcing vp/vs as fraction of G norm [scalar]
% eps_rhovs - weight for enforcing rho/vs as fraction of G norm [scalar]
% z_dampbot - depth below which to damp to starting model (km) [scalar]
% nit - number of iterations [scalar]
% nit_recalc_c - number of iterations after which to recalculate phase velocity and kernels [scalar]
% vp_vs - desired vp/vs scaling [M x 1 or scalar]
% rho_vs - desired rho/vs scaling [M x 1 or scalar]
% nmode - mode number [scalar]
%
% OUTPUTS
% finalmod - final surf96 model: dz (km), vp (km/s), vs (km/s), rho (kg/m^3) [M x 4]
% cpre - phase velocity predicted for final model [N x 1]
% vs_std - formal uncertainties on Vs [M x 1]
% vp_std - formal uncertainties on Vp [M x 1]
% rho_std - formal uncertainties on Rho [M x 1]
%
% jbrussell - 11/20/2022

eps_large = 1e9; % large weight to force constraint equation

% Calculate kernels for G matrix using SURF96
ifnorm = 0; % for plotting only
ifplot = 0;
[dcdvs, dcdvp, dudvs, dudvp, zkern, dcdrho, dudrho] = calc_kernel96(startmod, periods, 'R', ifnorm, ifplot,nmode);
G = [dcdvs' dcdvp' dcdrho'];

% Data weighting
% min_pct = 0.005; % minimum error percentage of observed
% I_error_too_small = find(cstd./cobs < min_pct);
% cstd(I_error_too_small) = cobs(I_error_too_small)*min_pct;
W = diag(1./cstd);

% Set up smoothing and damping matrices (applies directly to Vs only)
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
% z_brks = sort([discs; zdisc_Q(:)]);
z_brks = sort(discs);
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
H00(ind_h2o,ind_h2o) = H00(ind_h2o,ind_h2o)*eps_large;
h0(ind_h2o) = h0(ind_h2o)*eps_large;

% Add Vp and Rho dummy zeros
H00=[H00 zeros(size(H00)) zeros(size(H00))];
J00=[J00 zeros(size(J00)) zeros(size(J00))];
F00=[F00 zeros(size(F00)) zeros(size(F00))];

% Add Vp/Vs and Rho/Vs constraints
%                  Vs               Vp            Rho
VP_VS_mat = [-vp_vs.*eye(nlayer) eye(nlayer) zeros(nlayer)];
RHO_VS_mat = [-rho_vs.*eye(nlayer) zeros(nlayer) eye(nlayer)];
vp_vs_vec = zeros(nlayer,1);
rho_vs_vec = zeros(nlayer,1);
NR=norm(VP_VS_mat,1);
eps_vpvs0 = eps_vpvs*NA/NR;
NR=norm(RHO_VS_mat,1);
eps_rhovs0 = eps_rhovs*NA/NR;
% Set Vp of water layer
VP_h2o_mat = zeros(length(ind_h2o),3*nlayer);
for ii = 1:length(ind_h2o)
    VP_h2o_mat(ii,nlayer + ind_h2o) = 1;
end
vp_h2o_vec = 1.5 * ones(length(ind_h2o),1);
% Set Rho of water layer
RHO_h2o_mat = zeros(length(ind_h2o),3*nlayer);
for ii = 1:length(ind_h2o)
    RHO_h2o_mat(ii, 2*nlayer + ind_h2o) = 1;
end
rho_h2o_vec = 1.03 * ones(length(ind_h2o),1);

% combine all constraints
H = [H00*eps_H0; J00*eps_J0; F00*eps_F0; VP_VS_mat*eps_vpvs0; RHO_VS_mat*eps_rhovs0; VP_h2o_mat*eps_large; RHO_h2o_mat*eps_large];
h = [h0*eps_H0; j0*eps_J0; f0*eps_F0; vp_vs_vec*eps_vpvs0; rho_vs_vec*eps_rhovs0; vp_h2o_vec*eps_large; rho_h2o_vec*eps_large];

% Data vector
cstart = dispR_surf96(periods,startmod,nmode); % "predictions";
cpre = cstart;

% Least squares inversion
premod = startmod;
vs_pre = premod(:,3);
vp_pre = premod(:,2);
rho_pre = premod(:,4);
m_pre = [vs_pre; vp_pre; rho_pre];
clrs = jet(nit);
isfigure = 1;
clear vs
for ii = 1:nit
    
    % reformulate inverse problem so that constraints apply directly to model
    dc = cobs - cpre;
    d = dc + G*m_pre;
    
    % least squares
    F = [W*G*eps_data; H];
    f = [W*d*eps_data; h];
    m = (F'*F)\F'*f;
    vs = m(1:nlayer);
    vp = m(nlayer+1:2*nlayer);
    rho = m(2*nlayer+1:3*nlayer);
    vs(vs<eps) = 0;
    m = [vs; vp; rho];
    
    % update data vector
    dm = m-m_pre;
    dc_update = G * dm;
    cpre = cpre + dc_update;
    
    % update model
    m_pre = m_pre + dm;
    vs_pre = m_pre(1:nlayer);
    vp_pre = m_pre(nlayer+1:2*nlayer);
    rho_pre = m_pre(2*nlayer+1:3*nlayer);
    premod(:,3) = vs_pre;
    premod(:,2) = vp_pre;
    premod(:,4) = rho_pre;
    
    % Model uncertainties
    m_std = diag(inv(F'*F)).^(1/2);
    vs_std = m_std(1:nlayer);
    vp_std = m_std(nlayer+1:2*nlayer);
    rho_std = m_std(2*nlayer+1:3*nlayer);
    
    if mod(ii,nit_recalc_c)==0
        cpre = dispR_surf96(periods,premod,nmode);
        dc = cobs - cpre;
        
        ifnorm = 0; % for plotting only
        ifplot = 0;
        [dcdvs, dcdvp, dudvs, dudvp, zkern, dcdrho, dudrho] = calc_kernel96(startmod, periods, 'R', ifnorm, ifplot,nmode);
        G = [dcdvs' dcdvp' dcdrho'];
    end
    
    if isfigure
        if ii == 1 
            figure(2); clf; set(gcf,'position',[370   372   967   580]);
            
            subplot(2,2,[1 3]); box on; hold on;
            h2 = plotlayermods(startmod(:,1),startmod(:,3),'-k');
            h2.LineWidth = 4;
            
            subplot(2,2,2); box on; hold on;
            plot(periods,cstart,'-ok','linewidth',4);
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
        plot(periods,cpre,'-o','color',clrs(ii,:),'linewidth',2);
        errorbar(periods,cobs,2*cstd,'or','linewidth',4);
        xlabel('Period');
        ylabel('Phase Velocity');
        set(gca,'FontSize',18,'linewidth',1.5);
        % pause;
        drawnow
    end
end

finalmod = premod;
cpre = dispR_surf96(periods,finalmod,nmode);

end

