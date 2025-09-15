


clear all


f = figure ; 


% [ha, pos] = tight_subplot(2,3,0.005,0.1) ;
% 
% axes(ha(1))

metfile = '/Users/marq789/Downloads/restest/trimxy_corlasso.methagl.2018112900gefs09d4.base.M1.m1.20181129.153000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy


hold on
KK = 3;
%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])

nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_29N1530_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);


%%%%%%%%%%%%%%%%%%%%%%%%%



f = figure ; 

metfile = '/Users/marq789/Downloads/restest/trimxy_corlasso.methagl.2018120400gefs19d4.base.M1.m1.20181204.153000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy


hold on
KK = 3;
%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])

colorbar 

nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_04D1530_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);

%%%%%%%%%%%%%%%%%%%%%%%%%



f = figure ; 

metfile = '/Users/marq789/Downloads/restest/trimxy_corlasso.methagl.2019012200gefs18d4.base.M1.m1.20190122.150000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy


hold on
KK = 3;

x = [-65.2  -64.5   -64.5  -65.2]; % x-coordinates
y = [-32.5  -31.6   -32.5  -31.6 ]; % y-coordinates

x = [ -65.2 -64.5   -64.5  -65.2  ]; % x-coordinates
y = [ -32.5 -32.5 -31.6  -31.6]; % y-coordinates

% Create a polyshape object
pgon = polyshape(x, y);


figure

% Plot the polyshape
% plot(pgon,'FaceColor',[0 0 1]);
hold on

%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])


nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_22J1500_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);



%%%%%%%%%%%%%%%%%%%%%%%%%



f = figure ; 

metfile = '/Users/marq789/Downloads/restest/trimxy_corlasso.methagl.2019012300gefs18d4.base.M1.m1.20190123.140000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy


hold on
KK = 3;
%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])



nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_23J1400_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);

%%%%%%%%%%%%%%%%%%%%%%%%%


f = figure ; 

metfile = '/Users/marq789/Downloads/restest/trimxy_corlasso.methagl.2019012500eda07d4.base.M1.m1.20190125.150000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy


hold on
KK = 3;
%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])



nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_25J1500_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);





%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;

metfile = '/Users/marq789/Downloads/restest/corlasso_methagl_2019012900gefs11d4_base_M1.m1.20190129.150000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy

hold on
KK = 3;
%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])



nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_29J1500_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);



%%%%%%%%%%%%%%%%%%%%%%%%

f = figure;

metfile = '/Users/marq789/Downloads/restest/trimxy_corlasso.methagl.2019020800eda08d4.base.M1.m1.20190208.153000.nc';
dx = 100.;

Vagl = ncread(metfile,'VA');
Uagl = ncread(metfile,'UA');
% Vagl10 = ncread(metfile,'V10');
% Uagl10 = ncread(metfile,'U10');
Terr = ncread(metfile,'HGT');
Lon = ncread(metfile,'XLONG');
Lat = ncread(metfile,'XLAT');
[dudy dudx dudz] = gradient(Uagl,dx,dx,1);
[dvdy dvdx dvdz] = gradient(Vagl,dx,dx,1);
clear dudy dudz dvdx dvdz
Convagl = -( dudx +  dvdy );
clear dudx dvdy

hold on
KK = 3;
%contourf(Lon,Lat, Convagl10(:,:),[0.001:0.0005:0.012],'LineColor','none')
contourf(Lon,Lat, Convagl(:,:,KK), [0.001:0.0005:0.012],'LineColor','none')
colormap(flipud(autumn(22)))
caxis([0.001 0.012])
th = 50;
scal = 0.5;
quiver(Lon(1:th:end,1:th:end),Lat(1:th:end,1:th:end),scal*Uagl(1:th:end,1:th:end,KK),scal*Vagl(1:th:end,1:th:end,KK),'k','linewidth',1)
contour(Lon,Lat,Terr,[1500:250:2500],'Color','b','LineWidth',0.75)  %[0.8 0.5 0]
axis equal
axis([-65.1 -64.6 -32.5 -31.8])



nfraclab = strcat('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/D4_08F1530_conv300m.eps');
EPSprint = horzcat('print -painters -depsc ',nfraclab);
eval([EPSprint]);

