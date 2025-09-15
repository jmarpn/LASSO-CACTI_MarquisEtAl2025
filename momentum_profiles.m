



clear


%{

metfile = '/Users/marq789/Downloads/29Nov100m/20181129/gefs09/base/les/subset_d4/f15min/trimxy_asl_corlasso.methasl.2018112900gefs09d4.base.M1.m1.20181129.144500.nc' ;
cloudfile = '/Users/marq789/Downloads/29Nov100m/20181129/gefs09/base/les/subset_d4/f15min/trimxy_asl_corlasso.cldhaasl.2018112900gefs09d4.base.M1.m1.20181129.144500.nc' ;

dt = '29Nov-1445' ;

metfile = '/Users/marq789/Downloads/29Nov100m/20181129/gefs09/base/les/subset_d4/f15min/trimxy_asl_corlasso.methasl.2018112900gefs09d4.base.M1.m1.20181129.120000.nc' ;
cloudfile = '/Users/marq789/Downloads/29Nov100m/20181129/gefs09/base/les/subset_d4/f15min/trimxy_asl_corlasso.cldhaasl.2018112900gefs09d4.base.M1.m1.20181129.120000.nc' ;

dt = '29Nov-1200' ;

%}



%{

metfile = '/Users/marq789/Downloads/29Jan100m/20190129/gefs11/base/les/subset_d4/f15min/trimxy_asl_corlasso.methasl.2019012900gefs11d4.base.M1.m1.20190129.151500.nc' ;
cloudfile = '/Users/marq789/Downloads/29Jan100m/20190129/gefs11/base/les/subset_d4/f15min/trimxy_asl_corlasso.cldhaasl.2019012900gefs11d4.base.M1.m1.20190129.151500.nc' ;

dt = '29Jan-1515' ;

metfile   = '/Users/marq789/Downloads/29Jan100m/20190129/gefs11/base/les/subset_d4/f15min/trimxy_asl_corlasso.methasl.2019012900gefs11d4.base.M1.m1.20190129.120000.nc' ;
cloudfile = '/Users/marq789/Downloads/29Jan100m/20190129/gefs11/base/les/subset_d4/f15min/trimxy_asl_corlasso.cldhaasl.2019012900gefs11d4.base.M1.m1.20190129.120000.nc' ;

dt = '29Jan-1200' ;

%}



%{

metfile   = '/Users/marq789/Downloads/04Dec100m/trimxy_asl_corlasso.methasl.2018120400gefs19d4.base.M1.m1.20181204.144500.nc' ;
cloudfile = '/Users/marq789/Downloads/04Dec100m/trimxy_asl_corlasso.cldhaasl.2018120400gefs19d4.base.M1.m1.20181204.144500.nc' ;

dt = '4Dec-1445' ;

metfile   = '/Users/marq789/Downloads/04Dec100m/trimxy_asl_corlasso.methasl.2018120400gefs19d4.base.M1.m1.20181204.120000.nc' ;
cloudfile = '/Users/marq789/Downloads/04Dec100m/trimxy_asl_corlasso.cldhaasl.2018120400gefs19d4.base.M1.m1.20181204.120000.nc' ;

dt = '4Dec-1200' ;
%}


%{

metfile   = '/Users/marq789/Downloads/22Jan100m/trimxy_asl_corlasso.methasl.2019012200gefs18d4.base.M1.m1.20190122.155500.nc' ;
cloudfile = '/Users/marq789/Downloads/22Jan100m/trimxy_asl_corlasso.cldhaasl.2019012200gefs18d4.base.M1.m1.20190122.155500.nc' ;

dt = '22Jan-1555' ;

metfile   = '/Users/marq789/Downloads/22Jan100m/trimxy_asl_corlasso.methasl.2019012200gefs18d4.base.M1.m1.20190122.120000.nc' ;
cloudfile = '/Users/marq789/Downloads/22Jan100m/trimxy_asl_corlasso.cldhaasl.2019012200gefs18d4.base.M1.m1.20190122.120000.nc' ;

dt = '22Jan-1200' ;

%}



%{

metfile   = '/Users/marq789/Downloads/08Feb100m/trimxy_asl_corlasso.methasl.2019020800eda08d4.base.M1.m1.20190208.153000.nc' ;
cloudfile = '/Users/marq789/Downloads/08Feb100m/trimxy_asl_corlasso.cldhaasl.2019020800eda08d4.base.M1.m1.20190208.153000.nc' ;

dt = '08Feb-1530' ;

metfile   = '/Users/marq789/Downloads/08Feb100m/trimxy_asl_corlasso.methasl.2019020800eda08d4.base.M1.m1.20190208.120000.nc' ;
cloudfile = '/Users/marq789/Downloads/08Feb100m/trimxy_asl_corlasso.cldhaasl.2019020800eda08d4.base.M1.m1.20190208.120000.nc' ;

dt = '08Feb-1200' ;

%}



%{

metfile   = '/Users/marq789/Downloads/25Jan100m/trimxy_asl_corlasso.methasl.2019012500eda07d4.base.M1.m1.20190125.153000.nc' ;
cloudfile = '/Users/marq789/Downloads/25Jan100m/trimxy_asl_corlasso.cldhaasl.2019012500eda07d4.base.M1.m1.20190125.153000.nc' ;

dt = '25Jan-1530' ;

metfile   = '/Users/marq789/Downloads/25Jan100m/trimxy_asl_corlasso.methasl.2019012500eda07d4.base.M1.m1.20190125.120000.nc' ;
cloudfile = '/Users/marq789/Downloads/25Jan100m/trimxy_asl_corlasso.cldhaasl.2019012500eda07d4.base.M1.m1.20190125.120000.nc' ;

dt = '25Jan-1200' ;

%}




imout = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/';

U = ncread(metfile,'UA');
W = ncread(metfile,'WA');
QV = ncread(metfile,'QVAPOR');
Terr = ncread(metfile,'HGT');
Zasl = ncread(metfile,'HAMSL');
QC = ncread(cloudfile,'QCLOUD');

figure; contourf(Terr',20)

Uxz = permute(U,[1 3 2]);
Wxz = permute(W,[1 3 2]);
Cxz = permute(QC,[1 3 2]);
Qxz = permute(QV,[1 3 2]);
U0xz = Uxz;  U0xz(:) = 0;


dualpol_colmap
ff = figure('position',[1161,370,435,440]); hold on
%pgo = polyshape([0 0 60 60],[0 60 60 0]);  plot(pgo,'FaceColor',[0 0 0 ]);
contourf(([1:951]*100)/1000, Zasl/1000,  mean( Uxz(:,:,1000:1200),3,'omitnan'  )',20,'LineColor','none') ;
colormap(flipud(pepsi2))
caxis([-10 10])
colorbar
hold on
contour(([1:951]*100)/1000, Zasl/1000,  mean( Cxz(:,:,1000:1200),3,'omitnan'  )',[0.00001 0.00001],'k','LineWidth',1.5 )
thin = 10;
thin2 = 3;
contour(([1:951]*100)/1000, Zasl/1000,  mean( Wxz(:,:,1000:1200),3,'omitnan'  )',[0.5 0.5],'m','LineWidth',1 )
% MU = mean(U0xz(:,:,1000:1200),3,'omitnan');
% MW = mean(Wxz(:,:,1000:1200),3,'omitnan') ;
% quiver( ([1:thin:951]*100)/1000 ,   Zasl(1:thin2:49)/1000,   MU(1:thin:951,1:thin2:49)'   ,     MW(1:thin:951,1:thin2:49)'   )
axis([10 50 0 6])
title(dt)
saveas(ff,[imout,dt,'_U.png'])


dualpol_colmap
fq = figure('position',[1161,370,435,440]); 
%pgo = polyshape([0 0 60 60],[0 60 60 0]);  plot(pgo,'FaceColor',[0 0 0]);
contourf(([1:951]*100)/1000, Zasl/1000,  mean( Qxz(:,:,1000:1200),3,'omitnan'  )',20,'LineColor','none') ;
colormap((melyel2))
caxis([.000 .012])
colorbar
hold on
contour(([1:951]*100)/1000, Zasl/1000,  mean( Cxz(:,:,1000:1200),3,'omitnan'  )',[0.00001 0.00001],'k','LineWidth',1.5 )
thin = 10;
thin2 = 3;
contour(([1:951]*100)/1000, Zasl/1000,  mean( Wxz(:,:,1000:1200),3,'omitnan'  )',[0.5 0.5],'m','LineWidth',1 )
% MU = mean(U0xz(:,:,1000:1200),3,'omitnan');
% MW = mean(Wxz(:,:,1000:1200),3,'omitnan') ;
% quiver( ([1:thin:951]*100)/1000 ,   Zasl(1:thin2:49)/1000,   MU(1:thin:951,1:thin2:49)'   ,     MW(1:thin:951,1:thin2:49)'   )
axis([10 50 0 6])
title(dt)
saveas(fq,[imout,dt,'_QV.png'])









clear

%{
metlist = ls(['/Users/marq789/Downloads/04Dec100m/trimxy_asl_corlasso.methasl.2018120400gefs19d4.base.M1.m1.20181204.*nc']);
metlist = split(metlist) ;  metlist(1,:) = [];   metlist(end,:) = [];
[ta tb] = size(metlist) ; clear tb ;

cldlist = ls(['/Users/marq789/Downloads/04Dec100m/trimxy_asl_corlasso.cldhaasl.2018120400gefs19d4.base.M1.m1.20181204.*nc']);
cldlist = split(cldlist) ;  cldlist(1,:) = [];   cldlist(end,:) = [];
[ta tb] = size(cldlist) ; clear tb ;
%}



metlist = ls(['/Users/marq789/Downloads/29Nov100m/20181129/gefs09/base/les/subset_d4/f15min/more/trimxy_asl_corlasso.methasl.2018112900gefs09d4.base.M1.m1.20181129.*nc']);
metlist = split(metlist) ;  metlist(1,:) = [];   metlist(end,:) = [];
[ta tb] = size(metlist) ; clear tb ;

cldlist = ls(['/Users/marq789/Downloads/29Nov100m/20181129/gefs09/base/les/subset_d4/f15min/more/trimxy_asl_corlasso.cldhaasl.2018112900gefs09d4.base.M1.m1.20181129.*nc']);
cldlist = split(cldlist) ;  cldlist(1,:) = [];   cldlist(end,:) = [];
[ta tb] = size(cldlist) ; clear tb ;



imout = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/images/';

dualpol_colmap

for  t = 1:1:ta

    metfile = char(metlist(t)) ;   cloudfile = char(cldlist(t) );

    dt = [ metfile(end-17:end-10),'-',metfile(end-8:end-3)]

    U = ncread(metfile,'UA');
    W = ncread(metfile,'WA');
    QV = ncread(metfile,'QVAPOR');
    Terr = ncread(metfile,'HGT');
    Zasl = ncread(metfile,'HAMSL');
    QC = ncread(cloudfile,'QCLOUD');

    % figure; contourf(Terr',20)

    Uxz = permute(U,[1 3 2]);
    Wxz = permute(W,[1 3 2]);
    Cxz = permute(QC,[1 3 2]);
    Qxz = permute(QV,[1 3 2]);
    U0xz = Uxz;  U0xz(:) = 0;


    Nn = 1250;
    Ss = 800;

    Nn = 1200;
    Ss = 1000;

    ff = figure('position',[883,483,713,327]);
    %pgo = polyshape([0 0 60 60],[0 60 60 0]);  plot(pgo,'FaceColor',[0 0 0]);
    contourf(([1:951]*100)/1000, Zasl/1000,  mean( Qxz(:,:,Ss:Nn),3,'omitnan'  )',20,'LineColor','none') ;
    colormap((melyel2))
    caxis([.001 .011])
    colorbar
    hold on
    contour(([1:951]*100)/1000, Zasl/1000,  mean( Qxz(:,:,Ss:Nn),3,'omitnan'  )',[0.006 0.0065 0.007],'-.','LineColor',[0 0.3 0],'LineWidth',1.5) ;
    contour(([1:951]*100)/1000, Zasl/1000,  mean( Cxz(:,:,Ss:Nn),3,'omitnan'  )',[0.00001 0.00001],'k','LineWidth',1.5 )
    thin = 10;
    thin2 = 3;
    contour(([1:951]*100)/1000, Zasl/1000,  mean( Wxz(:,:,Ss:Nn),3,'omitnan'  )',[0.5 0.5],'m','LineWidth',1 )
    
    %MU = mean(U0xz(:,:,1000:1200),3,'omitnan');
    MU = mean(Uxz(:,:,Ss:Nn),3,'omitnan');
    MW = mean(Wxz(:,:,Ss:Nn),3,'omitnan') ;
    fact = 0.2
    quiver( ([1:thin:951]*100)/1000 ,   Zasl(1:thin2:49)/1000,   MU(1:thin:951,1:thin2:49)' * fact  ,     MW(1:thin:951,1:thin2:49)' * fact,'off' ,'k')
    axis([10 50 0.5 5])
    
    title(dt)
    %saveas(ff,[imout,dt,'_QV.png'])

    outlab = horzcat(imout,'/Nov29_moisture_xz_',dt,'.eps') ;
    EPSprint = horzcat('print -painters -depsc ',outlab);
    eval([EPSprint]);




end























