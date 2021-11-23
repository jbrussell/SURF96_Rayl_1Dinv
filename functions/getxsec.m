function [ V,Z,DISTkm,lats,lons ] = getxsec( v_mat,z_mat,xi,yi,latxsec,lonxsec, z_max )
% Get model parameter v along cross-section defined by lat/lon endpoints.
%
% INPUT:
% v_mat - parameter Vs, Vp, rho [Nlat x Nlon x Nz]
% z_mat - depth node values     [Nlat x Nlon x Nz]
% xi - latitude grid mesh [Nlat x Nlon]
% yi - longitude grid mesh [Nlat x Nlon]
% latxsec - latitude of endpoints for cross-section [2 x 1]
% lonxsec - longitude of endpoints for cross-section [2 x 1]
%
% OUTPUT:
% V - output velocity profile
% Z - output depth nodes
% DIST - output distance along profile
% lats - longitudes along xsection
% lons - latitudes along xsections
% 
% jbrussell 6/29/2020
%
isplot = 0;

% Get profile coordinates
[~,ix(1)] = min(abs(xi(:,1)-latxsec(1)));
[~,ix(2)] = min(abs(xi(:,1)-latxsec(2)));
[~,iy(1)] = min(abs(yi(1,:)-lonxsec(1)));
[~,iy(2)] = min(abs(yi(1,:)-lonxsec(2)));
[cx,cy,~] = improfile(squeeze(v_mat(1,1,:)),ix,iy);
cx = round(cx);
cy = round(cy);

% define lats and lons along profile
lats = xi(cx,1);
lons = yi(1,cy)';
dist_km = deg2km(distance(lats(1),lons(1),lats,lons));

% define indeces of output model
I = find(squeeze(z_mat(cx(1),cy(1),:))<=z_max);

% Loop over all points along profile
V = zeros(length(I),length(cx));
Z = V;
DISTkm = Z;
for idis = 1:length(cx)
    v = squeeze(v_mat(cx(idis),cy(idis),:));
    z = squeeze(z_mat(cx(idis),cy(idis),:));
%     % Remove repeated nodes
%     Irep = [false; (diff(model(idis).z)==0 & diff(model(idis).vs)==0)];
%     dis = model(1).dis;
%     vs = [flip(model(idis).vs(~Irep))];
%     vp = [flip(model(idis).vp(~Irep))];
%     rho = [flip(model(idis).rho(~Irep))];
%     z = [flip(model(idis).z(~Irep))];

    % Add an extra discontinuity if no water column otherwise plotting functions get confused
    if v(1)~=0
        v = [0; 0; v];
%         vp = [1500; 1500; vp];
%         rho = [1030; 1030; rho];
        z = [-1; 0; z];
    end
    
    V(:,idis) = v(I);
    Z(:,idis) = z(I);
    DISTkm(:,idis) = ones(size(z(I)))*dist_km(idis);
end

if isplot
    figure(99); clf;
    set(gcf,'color','w')
    subplot(3,1,[2 3]); box on; hold on;
    plot(lons,lats,'-ok')
    plot(lonxsec,latxsec,'-r');
    xlim([min(yi(:)),max(yi(:))]);
    ylim([min(xi(:)),max(xi(:))]);
    
    subplot(3,1,1); hold on;
    V(V==0) = nan;
    contourf(DISTkm,Z,V/1000,100,'edgecolor','none'); hold on;
    [C,h] = contour(DISTkm,Z,V/1000,[3:0.2:4.4],'-w','linewidth',1);
    clabel(C,h,'color',[1 1 1]);
    xlabel('Distance (km)');
    ylabel('Depth (km)');
    set(gca,'YDir','reverse','linewidth',1.5,'Layer','top','fontsize',16,'TickDir','out');
    set(gca,'Color',[0.8 0.8 0.8]);
    ylim([0 z_max]);
    xlim([min(DISTkm(:)) max(DISTkm(:))]);
    cb = colorbar;
    colormap(tomo_cmap(100));
    caxis([1 max(V(:)/1000)]);
    set(cb,'linewidth',1.5);
    ylabel(cb,'Vs (km/s)');
end

end

