
clear all

% load /Users/marq789/Downloads/N.xyz
% load /Users/marq789/Downloads/O.xyz
load /Volumes/passport/Argentina_geo/N.xyz
load /Volumes/passport/Argentina_geo/O.xyz



nlon = N(1:1:end,1);
nlat = N(1:1:end,2);
nele = N(1:1:end,3);
olon = O(1:1:end,1);
olat = O(1:1:end,2);
oele = O(1:1:end,3);

clear O N

save Terr_for_Adam.mat

plat = vertcat(nlat,olat);
plon = vertcat(nlon,olon);
pele = vertcat(nele,oele);

la = [-34:0.01:-30];
lo = [-67:0.01:-63];

la = [-38:0.01:-26];
lo = [-72:0.01:-58];

[lat lon] = meshgrid(la, lo);


terr = griddata(plat,plon,pele,lat,lon);

terr2 = terr;

zero = find(terr == 0);
terr2(zero) = -100;

save('/Users/marq789/Documents//PROJECTS/ICLASS/LASSO_stuff/Terr_for_JIMLASSO.mat','terr','terr2','plat','plon','lat','lon','la','lo','pele','-v7.3')

clear all

load('/Users/marq789/Documents//PROJECTS/ICLASS/LASSO_stuff/Terr_for_JIMLASSO.mat')


terr(terr<= 0) = -01;

figure
contourf(lon,lat,terr,[-100:100:3000],'LineColor','none'); axis equal;
%colormap(flipud(copper))
demcmap(terr)
hold on
%contourf(lon,lat,terr,[-10:0],'LineColor','none'); 
borders('argentina','k','LineWidth',2)
borders('chile','k','LineWidth',2)
borders('uruguay','k','LineWidth',2)
borders('paraguay','k','LineWidth',2)
borders('brazil','k','LineWidth',2)
hold on

% plot(-65,-31,'k.','MarkerSize',1)  %wrf
% plot(-65,-32,'k.','MarkerSize',1)  %wrf
% plot(-64,-31,'k.','MarkerSize',1)  %wrf
% plot(-64,-32,'k.','MarkerSize',1)  %wrf

plot(-70,-37,'+r','MarkerSize',8)      %D1
plot(-59,-37,'+r','MarkerSize',8)      %D1
plot(-69.4,-27.8,'+r','MarkerSize',8)  %D1
plot(-59.7,-27.8,'+r','MarkerSize',8)  %D1

plot(-68,-36,'+m','MarkerSize',8)    %D2
plot(-60.7,-36,'+m','MarkerSize',8)  %D2
plot(-67.7,-29,'+m','MarkerSize',8)  %D2
plot(-61,-29,'+m','MarkerSize',8)    %D2

plot(-66.3,-34,'+c','MarkerSize',8)    %D3
plot(-62.2,-34,'+c','MarkerSize',8)    %D3
plot(-66.2,-30.1,'+c','MarkerSize',8)  %D3
plot(-62.3,-30.1,'+c','MarkerSize',8)  %D3

plot(-65.3,-33.25,'+w','MarkerSize',8)    %D4
plot(-63.1,-33.25,'+w','MarkerSize',8)    %D4
plot(-65.3,-30.8,'+w','MarkerSize',8)  %D4
plot(-63.1,-30.8,'+w','MarkerSize',8)  %D4

% plot(-64.192,-31.441,'dk','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %COR rad
% plot(-64.213,-31.312,'ok','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %COR sonde
% plot(-64.170,-31.635,'dk','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %CSU rad
% plot(-64.223,-31.987,'dk','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %COW
% %plot(-64.500,-31.421,'ko','MarkerSize',5,'MarkerFaceColor','w','LineWidth',1)  %VCP
% plot(-65.150,-31.951,'ok','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %villa dol
% plot(-64.927,-32.052,'^k','MarkerSize',6,'MarkerFaceColor','w','LineWidth',1)  %met1
% plot(-64.881,-32.067,'^k','MarkerSize',6,'MarkerFaceColor','w','LineWidth',1)  %met2
% hold on


% 
 colorbar
 caxis([-100 2600])
% cbh=colorbar('h');
% set(cbh,'YTick',[250:500:2750])
% figout = horzcat('/Users/marq789/Documents/PROJECTS/RELAMPAGO/TrappEtAl_KumjianEtAl_VarbleetAl_items/CACTI_BAMS_colorbar.eps');
% EPSprint = horzcat('print -painters -depsc ',figout);
% eval([EPSprint]);

% %wheels
% th = 0:pi/30:2*pi;
% degth = 0:2:360;
% [latcirc100,loncirc100] = reckon(-32.126,-64.729,km2deg(100),degth) ;
% plot(loncirc100, latcirc100,'r','LineWidth',1.5); hold on
% [latcirc20,loncirc20] = reckon(-32.126,-64.729,km2deg(20),degth) ;
% plot(loncirc20, latcirc20,'g','LineWidth',1.5);

% 
% %spokes
% for degth = 0:30:360
%     leng = 0:1:100
%     [latspoke,lonspoke] = reckon(-32.126,-64.729,km2deg(leng),degth) ;
%     plot(lonspoke, latspoke,'r','LineWidth',1.5); hold on
%     leng = 0:1:20
%     [latspoke,lonspoke] = reckon(-32.126,-64.729,km2deg(leng),degth) ;
%     plot(lonspoke, latspoke,'g','LineWidth',1.5); hold on
% end



% plot(-64.729,-32.126,'kd','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %amf
% plot(-64.192,-31.441,'dk','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %COR rad
% plot(-64.213,-31.312,'ok','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %COR sonde
% plot(-64.170,-31.635,'dk','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %CSU rad
% plot(-64.223,-31.987,'dk','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %COW
% plot(-63.883,-31.677,'ks','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %dial
% plot(-65.150,-31.951,'ok','MarkerSize',8,'MarkerFaceColor','w','LineWidth',1)  %villa dol
% plot(-64.927,-32.052,'^k','MarkerSize',6,'MarkerFaceColor','w','LineWidth',1)  %met1
% plot(-64.881,-32.067,'^k','MarkerSize',6,'MarkerFaceColor','w','LineWidth',1)  %met2
% plot(-64.496,-31.421,'ko','MarkerSize',4,'MarkerFaceColor','k','LineWidth',1)  %vcp
% plot(-64.754,-32.106,'ko','MarkerSize',4,'MarkerFaceColor','k','LineWidth',1)  %yacanto
hold on

axis equal
axis([-72 -58 -38 -26])
% xticks([-65.5 -65 -64.5 -64])
% yticks([-33 -32.5 -32 -31.5])
grid on
% ax = gca;
% ax.LineWidth = 0.75;
% ax.GridColor = 'k';
% ax.GridAlpha = 0.5;

xlabel('Lon [deg]')
ylabel('Lat [deg]')

figout = horzcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/LASSO_domain.eps');
EPSprint = horzcat('print -painters -depsc ',figout);
eval([EPSprint]);





