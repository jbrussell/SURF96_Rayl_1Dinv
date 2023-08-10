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

%% Starting model

Nmodes = 3; % number of modes to plot
zh2o = 5; % water depth
zsed = 0.5; %0.25; % sediment thickness

ylims = [0 100];

vec_T = [1:0.05:15];

% Make homogeneous starting model
dzz = 1; % km
z = [0 zh2o [zh2o+zsed : dzz : 100]];
dz = [diff(z) 0];
vs = [0 1.0 3.2*ones(size([zh2o+zsed : dzz : 100]))];
vp = [1.5 3.0 5.8*ones(size([zh2o+zsed : dzz : 100]))]; %1.75*vs;
rho = [1.0 1.8 3*ones(size([zh2o+zsed : dzz : 100])) ];
startmod = [dz(:), vp(:), vs(:), rho(:)];
startmod(1:end-1,:);
% discs = [zmoho]; % [km] depth to sharp discontinuities
discs = []; % for this simple example, don't allow discontinuities

% % Make homogeneous starting model
% z = [0 zh2o zh2o+zsed 10 20 30 40 50 60 70 80 90 100];
% dz = [diff(z) 0];
% vs = [0 1.0 3.2 3.2 3.2 3.2 3.2 3.2 3.2 3.2 3.2 3.2 4.5];
% vp = [1.5 3.0 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8 5.8]; %1.75*vs;
% rho = [1.0 1.8 3 3 3 3 3 3 3 3 3 3 3 ];
% startmod = [dz(:), vp(:), vs(:), rho(:)];
% startmod(1:end-1,:);
% % discs = [zmoho]; % [km] depth to sharp discontinuities
% discs = []; % for this simple example, don't allow discontinuities

% % Make homogeneous starting model
% z = [0 20 100];
% dz = [diff(z) 0];
% vs = [3 5 5];
% vp = [7 7 7]; %1.75*vs;
% rho = vp / 2.5;
% startmod = [dz(:), vp(:), vs(:), rho(:)];
% % discs = [zmoho]; % [km] depth to sharp discontinuities
% discs = []; % for this simple example, don't allow discontinuities

clrs = jet(Nmodes);
for ii = 1:Nmodes
    
    Nmode = ii-1;
    
    c_R = dispR_surf96(vec_T,startmod,Nmode); % "predictions";

    grv_R = dispR_surf96(vec_T,startmod,Nmode,'U'); % "predictions";
        
    disper = calc_amplification_ellipticity_gamma(vec_T,startmod,'R',Nmode);
    
    if isempty(grv_R)
        continue
    end

    %%
    if ii == 1
        figure(1); clf;
        set(gcf,'position',[70          81        1526         842]);
    end

    subplot(2,3,1);
    box on; hold on;
    h = plotlayermods(startmod(:,1),startmod(:,3),'-b');
    h.LineWidth = 2;
    xlabel('Vs (km/s)');
    ylabel('Depth (km)');
    title('Starting Model');
    set(gca,'FontSize',18,'linewidth',1.5);
    ylim(ylims);

    subplot(2,3,2);
    box on; hold on;
    plot(vec_T(1:length(c_R)),c_R,'-','color',clrs(ii,:),'linewidth',2);
    plot(vec_T(1:length(grv_R)),grv_R,'--','color',clrs(ii,:),'linewidth',2);
    xlabel('Period (s)');
    ylabel('Phase Velocity (km/s)');
    set(gca,'FontSize',18,'linewidth',1.5);
    lgd{ii} = num2str(Nmode);
    xlim([min(vec_T) max(vec_T)]);
    
    subplot(2,3,3);
    box on; hold on;
    plot(vec_T(1:length(disper.gamma)),disper.gamma,'-','color',clrs(ii,:),'linewidth',2);
    xlabel('Period (s)');
    ylabel('Gamma (1/km)');
    set(gca,'FontSize',18,'linewidth',1.5);
    lgd{ii} = num2str(Nmode);
    xlim([min(vec_T) max(vec_T)]);

    subplot(2,3,5);
    box on; hold on;
    plot(vec_T(1:length(disper.A_R)),disper.A_R,'-','color',clrs(ii,:),'linewidth',2);
    xlabel('Period (s)');
    ylabel('A_R');
    set(gca,'FontSize',18,'linewidth',1.5);
    xlim([min(vec_T) max(vec_T)]);
    
    subplot(2,3,6);
    box on; hold on;
    h4(ii) = plot(vec_T(1:length(disper.RZ)),disper.RZ,'-','color',clrs(ii,:),'linewidth',2);
    xlabel('Period (s)');
    ylabel('Ellipticity R/Z');
    set(gca,'FontSize',18,'linewidth',1.5);
    xlim([min(vec_T) max(vec_T)]);
end

pos = get(gca,'position');
legend(h4,lgd,'location','eastoutside');
set(gca,'position',pos);

%% Eigenfunctions

ipers = [1 50 100 200];

linetype = {'-','--',':'};
for ii = 1:Nmodes
    
    Nmode = ii-1;
    
    eig = calc_eigenfunctions96(vec_T,startmod,'R',Nmode);
        

    %%
    if ii == 1
        figure(2); clf;
        set(gcf,'position',[70          81        1526         842]);
    end
    
    for iper = 1:length(ipers)
    
        subplot(1,length(ipers),iper);
        box on; hold on;
        plot(eig.uz(:,ipers(iper)),eig.z,linetype{ii},'color',[1 0 0],'linewidth',2);
        plot(eig.ur(:,ipers(iper)),eig.z,linetype{ii},'color',[0 0 1],'linewidth',2);
        plot(eig.tz(:,ipers(iper)),eig.z,linetype{ii},'color',[1 0 1],'linewidth',2);
        plot(eig.tr(:,ipers(iper)),eig.z,linetype{ii},'color',[0 1 1],'linewidth',2);
        xlabel('Eigenfunction');
        ylabel('Depth (km)');
        title([num2str(eig.periods(ipers(iper))),' s']);
        set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
        ylim(ylims);
        legend({'Uz','Ur','Tz','Tr'},'location','southeast');
    end
    
end

%% Ellipticity

for ii = 1:Nmodes
    
    Nmode = ii-1;
    
    eig = calc_eigenfunctions96(vec_T,startmod,'R',Nmode);
        

    %%
    if ii == 1
        figure(3); clf;
        set(gcf,'position',[70          81        1526         842]);
    end
    
    for iper = 1:length(ipers)
    
        subplot(1,length(ipers),iper);
        box on; hold on;
        plot(eig.ur(:,ipers(iper)) ./ eig.uz(:,ipers(iper)),eig.z,linetype{ii},'color',[1 0 0],'linewidth',2);
        plot(eig.tr(:,ipers(iper)) ./ eig.tz(:,ipers(iper)),eig.z,linetype{ii},'color',[0 0 1],'linewidth',2);
        xlabel('Ellipticity');
        ylabel('Depth (km)');
        title([num2str(eig.periods(ipers(iper))),' s']);
        set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
        ylim(ylims);
        legend({'Ur/Uz','Tr/Tz'},'location','southeast');
    end
    
end


