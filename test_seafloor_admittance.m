% Plot seafloor admittance following Ruan et al. (2014)
%
% admittance = displacement / pressure
% adm = u_z(seafloor) / t_zz(seafloor)
%
% jbrussell - 7/2024
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

%% Recreate Figure 6

Nmodes = 1; %3; % number of modes to plot

ylims = [0 10];

% Periods for plotting eigenfunctions
vec_T = [6, 8, 10, 20]; % s

% Periods over which to calculate admittance
vec_T_adm = linspace(5,50,100);

% Sediment Vs to loop over
vs_sed_vec = [0.28 0.48 0.98]; % km/s

figure(99); clf;
set(gcf,'position',[70          81        1526         842]);
lgd = {}; h6=[]; h6_adm=[];
clrs_adm = lines(length(vs_sed_vec));
for ised = 1:length(vs_sed_vec)
    
    % Make homogeneous starting model
    dz_h2o = 2.5; % km
    dz_sed = 0.6; %0.6; % km
    vs_sed = vs_sed_vec(ised); %0.48; % km/s
    vp_sed = vs_sed * 2.5; % km/s
    z_h2o = dz_h2o;
    z_sed = dz_h2o + dz_sed;

    dz = [dz_h2o  dz_sed  2     5     20     0];
    vs = [0       vs_sed  2.63  3.89  3.27   3.27];
    vp = [1.5     vp_sed  5     6.8   7.913  7.913];
    rho = [1.03   2.0     2.45  3.05  4.326  4.326];
    startmod = [dz(:), vp(:), vs(:), rho(:)];


    % Interpolate to finer layer grid
    dzz = 0.1; % [km] desired layer thicknesses
    layers_to_interpolate = [1 2 3 4]; % layer indices to interpolate to dzz
    Nlayers_per_lay = round(dz ./ dzz);
    dz_int = [];
    vp_int = [];
    vs_int = [];
    rho_int = [];
    for ilay = 1:length(dz)
        if ismember(ilay,layers_to_interpolate)
            dz_int = [dz_int; repmat(dzz,Nlayers_per_lay(ilay),1)];
            vp_int = [vp_int; repmat(vp(ilay),Nlayers_per_lay(ilay),1)];
            vs_int = [vs_int; repmat(vs(ilay),Nlayers_per_lay(ilay),1)];
            rho_int = [rho_int; repmat(rho(ilay),Nlayers_per_lay(ilay),1)];
        else
            dz_int = [dz_int; dz(ilay)];
            vp_int = [vp_int; vp(ilay)];
            vs_int = [vs_int; vs(ilay)];
            rho_int = [rho_int; rho(ilay)];
        end
    end
    % Replace Vp in seds with Hamilton 1979
    D = 0:dzz:dz_sed;
    vp_hamilton = 1.511 + 1.304*D - 0.741*D.^2 + 0.257*D.^3;
    Ised = find(vp_int == vp_sed);
    vp_int(Ised) = vp_hamilton(1:length(Ised));
    
    startmod_int = [dz_int(:), vp_int(:), vs_int(:), rho_int(:)];


    % Eigenfunctions

    clrs = lines(length(vec_T));

    linetype = {'-','--',':'};

    for ii = 1:Nmodes

        Nmode = ii-1;

        % Calculate admittance
        adm = calc_admittance96_Ruan14(vec_T_adm,startmod_int,Nmode);

        % Calculate phase velocity
        c_R = dispR_surf96(vec_T_adm,startmod,Nmode); % "predictions";
    %     disper = calc_amplification_ellipticity_gamma(vec_T_adm,startmod,'R',Nmode);
    %     c_R = disper.phv;

        eig = calc_eigenfunctions96(vec_T,startmod_int,'R',Nmode);
        
        for iper = 1:length(vec_T)

    %         if ipers(iper) > size(eig.uz,2)
    %             continue
    %         end

            subplot(2,3,[1 4]);
            box on; hold on;
            h6(iper) = plot(eig.uz(:,iper),eig.z,linetype{ised},'color',clrs(iper,:),'linewidth',2,'displayname', [num2str(eig.periods(iper)),' s']);
    %         plot(eig.ur(:,iper),eig.z,linetype{ii},'color',clrs(iper,:),'linewidth',2,'displayname', [num2str(eig.periods(iper)),' s']);
            xlabel('u_z (km)');
            ylabel('Depth (km)');
            title('Vertical Displacement');
            set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
            ylim(ylims);
            xvals = get(gca,'XLim');
            if ii == 1 && iper == 1
                plot(xvals,z_h2o*[1 1],'--k','linewidth',2);
                plot(xvals,z_sed*[1 1],'--k','linewidth',2);
            end
            lgd{iper} = [num2str(eig.periods(iper)),' s'];
            legend(h6,lgd,'location','northwest');

            subplot(2,3,[2 5]);
            box on; hold on;
            plot(eig.tz(:,iper),eig.z,linetype{ised},'color',clrs(iper,:),'linewidth',2,'displayname', [num2str(eig.periods(iper)),' s']);
    %         plot(eig.tr(:,iper),eig.z,linetype{ii},'color',clrs(iper,:),'linewidth',2,'displayname', [num2str(eig.periods(iper)),' s']);
            xlabel('\tau_{zz} (GPa)');
            ylabel('Depth (km)');
            title('Normal Stress');
            set(gca,'FontSize',18,'linewidth',1.5,'ydir','reverse');
            ylim(ylims);
            xvals = get(gca,'XLim');
            if ii == 1 && iper == 1
                plot(xvals,z_h2o*[1 1],'--k','linewidth',2);
                plot(xvals,z_sed*[1 1],'--k','linewidth',2);
            end
        end

        % Phase velocity
        subplot(2,3,3);
        box on; hold on;
        vec_f_adm = 1./vec_T_adm;
        plot(vec_f_adm,c_R,linetype{ised},'color',clrs_adm(ised,:),'linewidth',2,'displayname', [num2str(vs_sed),' km/s']);
        xlabel('Frequency (Hz)');
        ylabel('Phase velocity (km/s)');
        set(gca,'FontSize',18,'linewidth',1.5);
    %     ylim(ylims);
        xlim([0 max(vec_f_adm)]);
%         legend();

        % Admittance
        subplot(2,3,6);
        box on; hold on;
        plot(1./adm.periods(:),abs(adm.admittance),linetype{ised},'color',clrs_adm(ised,:),'linewidth',2,'displayname', [num2str(vs_sed),' km/s']);
        xlabel('Frequency (Hz)');
        ylabel('|Admittance| (m/Pa)');
        set(gca,'FontSize',18,'linewidth',1.5,'YScale','log');
    %     ylim(ylims);
        xlim([0 max(vec_f_adm)]);
        legend();

    end
end

