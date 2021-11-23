function ax = montecarl_plot(str)

setup_parameters;
lalim = parameters.lalim;
lolim = parameters.lolim;

final_mods=str.result.final_mods;
init_mods=str.result.init_mods;
vec_hs=str.result.vec_hs;
phv_fwds=str.result.phv_fwds;
errors=str.result.errors;
init_avg=str.result.init_avg;
final_avg=str.result.final_avg;
init_std=str.result.init_std;
final_std=str.result.final_std;
velT = str.datastr.periods;
phv = str.datastr.phv;
phvstd = str.datastr.phvstd;
test_N = length(str.result.errors);
depth_nodes = parameters.depth_nodes;


Nx = 2; Ny = 1;
sidegap = 0.10; topgap = 0.10; botgap = 0.10; vgap = 0.05; hgap = 0.05;
cbar_bot = 0.04;

width = (1 - vgap*(Nx-1)-2*sidegap)/Nx;
height = (1 - topgap - botgap - (Ny-1)*hgap)/Ny;

figure(57)
clf
set(gcf,'position',[ 118   290   800   600]);
set(gcf,'color','w');
ix = 1;
iy = 1;
left = sidegap + (iy-1)*(vgap+width);
bot = botgap + (ix-1)*(hgap+height);
ax(1) = subplot('position',[left,bot,width,height]);
hold on
for itest=1:test_N
	h = plotlayermods(vec_hs(:,itest),init_mods(:,itest));
	color = [0.6 0.6 0.6];
	set(h,'color',color);
end
avgvel = init_avg;
upvel = init_avg + init_std;
lowvel = init_avg - init_std;
plot(avgvel,-depth_nodes,'k','linewidth',2);
plot(upvel,-depth_nodes,'k--','linewidth',2);
plot(lowvel,-depth_nodes,'k--','linewidth',2);
title('Initial Models','fontsize',24)
ylabel('Depth (km)','fontsize',18)
xlabel('Shear Velocity (km/s)','fontsize',18)
xlim([1.5 5]);
ylim([-100 0]);
set(gca,'fontsize',15,'box','on','LineWidth',1.5);

ix = 1;
iy = 2;
left = sidegap + (iy-1)*(vgap+width);
bot = botgap + (ix-1)*(hgap+height);
ax(2) = subplot('position',[left,bot,width,height]);
hold on
for itest=1:test_N
	h = plotlayermods(vec_hs(:,itest),final_mods(:,itest));
	color = [0.6 0.6 0.6];
	set(h,'color',color);
end
avgvel = final_avg;
upvel = final_avg + final_std;
lowvel = final_avg - final_std;
plot(avgvel,-depth_nodes,'k','linewidth',2);
plot(upvel,-depth_nodes,'k--','linewidth',2);
plot(lowvel,-depth_nodes,'k--','linewidth',2);
title('Final Models','fontsize',22)
xlabel('Shear Velocity (km/s)','fontsize',18)
xlim([1.5 5]);
ylim([-100 0]);
set(gca,'fontsize',15,'box','on','LineWidth',1.5);

ix = 1;
iy = 2;
left = sidegap + (iy-1)*(vgap+width)+0.02;
bot = botgap + (ix-1)*(hgap+height);
ax(3) = axes('position',[left+width/7,bot+height/10,width/2,height/3]);
box on
hold on
for itest=1:test_N
	h = plot(velT(:), phv_fwds(:,itest));
	set(h,'color',[0.7 0.7 0.7],'linewidth',2);
end
% plot(velT, phv,'ko','markersize',10,'linewidth',2);
errorbar(velT,phv,phvstd,'ko','markerfacecolor',[0 0 0],'markersize',6,'linewidth',1.5);
% xlim([8 65])
xlim(ax(3),[min(velT)*0.95 max(velT)*1.05]);
xlabel('Period (s)');
ylabel('Phase Velocity (km/s)');
set(gca,'FontSize',13,'box','on','LineWidth',1.5);


if 1
    ix = 1;
    iy = 1;
    left = sidegap + (iy-1)*(vgap+width)+0.02;
    bot = botgap + (ix-1)*(hgap+height);
    ax(4) = axes('position',[left+width/9,bot+height/10,width/2,height/3]);
    worldmap(lalim,lolim);
    setm(gca,'FLineWidth',1.5,'FontSize',13)
    % drawpng
    mlat = str.lat; mlon = str.lon;
    plotm(mlat,mlon,'rx','markersize',15,'linewidth',2);
    plotm(mlat,mlon,'ro','markersize',15,'linewidth',2);
end
