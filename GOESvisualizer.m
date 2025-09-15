

clear

%%%%% terrain file:
load('Terr_for_me.mat','terr','lon','lat')
lon_terr = lon;
lat_terr = lat;
clear lat lon





R = 4  ;




%%%%%%  load LASSO swaths:
rootdir =  '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/';
RUNS = ['20181129_gefs09_base_d4';
        '20190122_gefs18_base_d4';
        '20190125_eda07_base__d4';
        '20190129_gefs11_base_d4'];

%swathmatfilein = [rootdir,'/', RUNS(R,:) ,'/SwathIDthermoOutput_c1p11h1p13tsmftopf_wraw_wmag025_v12_d4_',RUNS(R,:),'.mat'] ;
swathmatfilein = [rootdir,'/', RUNS(R,:) ,'/SwathIDthermoOutput_c1p11h1p13tsmftopf_wraw_wmag025_v12_d4_',RUNS(R,:),'.mat'] ;
load(swathmatfilein,'W_object','Conv_object','HHMMSS_all','Lat','Lon')
swathHH = str2num(HHMMSS_all(:,1:2)) ;
swathMM = str2num(HHMMSS_all(:,3:4)) ;
swath_hr = swathHH + swathMM/60;

if(R==1)

elseif(R==2)
    W_object(:,:,31:end,1400:end) = 0;
    W_object(:,:,49:end,1:200) = 0;

elseif(R==3)
    W_object(:,:,:,1200:end) = 0;
    W_object(:,:,1:29,1:100) = 0;

elseif(R==4)
    W_object(:,:,1:13,1:175) = 0;
    W_object(:,:,1:25,1:150) = 0;
    W_object(:,:,:,1200:end) = 0;

end



% ignore the 05,10,20,25,35,40,50,55 min analyses:
kill = find( swathMM == 5  |  swathMM == 10  |  swathMM == 20  |  swathMM == 25  | ...
             swathMM == 35  |  swathMM == 40  |  swathMM == 50  |  swathMM == 55 ) ;
swathHH(kill) = [];
swathMM(kill) = [];
swath_hr(kill) = [];
W_object(:,:,kill,:) = [];
Conv_object(:,:,kill,:) = [];
HHMMSS_all(kill,:) = [];




%%%% load LASSO cloud & w fields
lassocldfiles = ls(strcat(rootdir,'/', RUNS(R,:) , '/trimxy_asl_corlasso.cldhaasl.*nc'))
lassocld_list = split(lassocldfiles) ;
lassocld_list(end) = []  ;



%%%% load obs GOES satellite data:
satdir = '/Users/marq789/Downloads/GOES_cacti/';
if(R==1)
    satfiles = ls(strcat(satdir,'corvisstpx2drectg16v4minnisX1.parallaxcorrected.c1.20181129*'))
elseif(R==2)
    satfiles = ls(strcat(satdir,'corvisstpx2drectg16v4minnisX1.parallaxcorrected.c1.20190122*'))
elseif(R==3)
    satfiles = ls(strcat(satdir,'corvisstpx2drectg16v4minnisX1.parallaxcorrected.c1.20190125*'))
elseif(R==4)
    satfiles = ls(strcat(satdir,'corvisstpx2drectg16v4minnisX1.parallaxcorrected.c1.20190129*'))
end
sat_list = split(satfiles) ;
sat_list(end) = []  ;
%[sa sb] = size(sat_list); clear sb;
%%% trim out unwanteds from list:
KILL = [];
if(R==1)
    KILL = cat(1, find(contains(sat_list,'1000')), find(contains(sat_list,'1015')), find(contains(sat_list,'1030')), find(contains(sat_list,'1045')),...
                  find(contains(sat_list,'1100')), find(contains(sat_list,'1115')), find(contains(sat_list,'1130')), find(contains(sat_list,'1145')),...
                  find(contains(sat_list,'1200')), find(contains(sat_list,'1215')), find(contains(sat_list,'1230')), find(contains(sat_list,'1245')),...
                                                                                                                     ...
                                                                                                                     find(contains(sat_list,'1645')),...
                  find(contains(sat_list,'1700')), find(contains(sat_list,'1715')), find(contains(sat_list,'1730')), find(contains(sat_list,'1745')),...
                  find(contains(sat_list,'1800')), find(contains(sat_list,'1815')), find(contains(sat_list,'1830')), find(contains(sat_list,'1845')),...
                  find(contains(sat_list,'1900')), find(contains(sat_list,'1915')), find(contains(sat_list,'1930')), find(contains(sat_list,'1945')),...
                  find(contains(sat_list,'2000')), find(contains(sat_list,'2015')), find(contains(sat_list,'2030')), find(contains(sat_list,'2045')),...
                  find(contains(sat_list,'2100')), find(contains(sat_list,'2115')), find(contains(sat_list,'2130')), find(contains(sat_list,'2145')),...
                  find(contains(sat_list,'2200')), find(contains(sat_list,'2215')), find(contains(sat_list,'2230')), find(contains(sat_list,'2245')),...
                  find(contains(sat_list,'2300')), find(contains(sat_list,'2315')), find(contains(sat_list,'2330')), find(contains(sat_list,'2345'))    );
    sat_list(KILL) = [];

elseif(R==2)
    KILL = cat(1, find(contains(sat_list,'1000')), find(contains(sat_list,'1015')), find(contains(sat_list,'1030')), find(contains(sat_list,'1045')),...
                  find(contains(sat_list,'1100')), find(contains(sat_list,'1115')), find(contains(sat_list,'1130')), find(contains(sat_list,'1145')),...
                  find(contains(sat_list,'1200')), find(contains(sat_list,'1215')), find(contains(sat_list,'1230')), find(contains(sat_list,'1245')),...
                                                                                                                    ...
                                                                                                                     find(contains(sat_list,'1845')),...
                  find(contains(sat_list,'1900')), find(contains(sat_list,'1915')), find(contains(sat_list,'1930')), find(contains(sat_list,'1945')),...
                  find(contains(sat_list,'2000')), find(contains(sat_list,'2015')), find(contains(sat_list,'2030')), find(contains(sat_list,'2045')),...
                  find(contains(sat_list,'2100')), find(contains(sat_list,'2115')), find(contains(sat_list,'2130')), find(contains(sat_list,'2145')),...
                  find(contains(sat_list,'2200')), find(contains(sat_list,'2215')), find(contains(sat_list,'2230')), find(contains(sat_list,'2245')),...
                  find(contains(sat_list,'2300')), find(contains(sat_list,'2315')), find(contains(sat_list,'2330')), find(contains(sat_list,'2345'))    );
    sat_list(KILL) = [];

elseif(R==3)
    KILL = cat(1, find(contains(sat_list,'1000')), find(contains(sat_list,'1015')), find(contains(sat_list,'1030')), find(contains(sat_list,'1045')),...
                  find(contains(sat_list,'1100')), find(contains(sat_list,'1115')), find(contains(sat_list,'1130')), find(contains(sat_list,'1145')),...
                  find(contains(sat_list,'1200')), find(contains(sat_list,'1215')), find(contains(sat_list,'1230')), find(contains(sat_list,'1245')),...
                                                                                                                    ...
                                                                                                                     find(contains(sat_list,'1845')),...
                  find(contains(sat_list,'1900')), find(contains(sat_list,'1915')), find(contains(sat_list,'1930')), find(contains(sat_list,'1945')),...
                  find(contains(sat_list,'2000')), find(contains(sat_list,'2015')), find(contains(sat_list,'2030')), find(contains(sat_list,'2045')),...
                  find(contains(sat_list,'2100')), find(contains(sat_list,'2115')), find(contains(sat_list,'2130')), find(contains(sat_list,'2145')),...
                  find(contains(sat_list,'2200')), find(contains(sat_list,'2215')), find(contains(sat_list,'2230')), find(contains(sat_list,'2245')),...
                  find(contains(sat_list,'2300')), find(contains(sat_list,'2315')), find(contains(sat_list,'2330')), find(contains(sat_list,'2345'))    );
    sat_list(KILL) = [];

elseif(R==4)
    KILL = cat(1, find(contains(sat_list,'1000')), find(contains(sat_list,'1015')), find(contains(sat_list,'1030')), find(contains(sat_list,'1045')),...
                  find(contains(sat_list,'1100')), find(contains(sat_list,'1115')), find(contains(sat_list,'1130')), find(contains(sat_list,'1145')),...
                  find(contains(sat_list,'1200')), find(contains(sat_list,'1215')), find(contains(sat_list,'1230')), find(contains(sat_list,'1245')),...
                                                                                                                     ...
                                                   find(contains(sat_list,'1615')), find(contains(sat_list,'1630')), find(contains(sat_list,'1645')),...
                  find(contains(sat_list,'1700')), find(contains(sat_list,'1715')), find(contains(sat_list,'1730')), find(contains(sat_list,'1745')),...
                  find(contains(sat_list,'1800')), find(contains(sat_list,'1815')), find(contains(sat_list,'1830')), find(contains(sat_list,'1845')),...
                  find(contains(sat_list,'1900')), find(contains(sat_list,'1915')), find(contains(sat_list,'1930')), find(contains(sat_list,'1945')),...
                  find(contains(sat_list,'2000')), find(contains(sat_list,'2015')), find(contains(sat_list,'2030')), find(contains(sat_list,'2045')),...
                  find(contains(sat_list,'2100')), find(contains(sat_list,'2115')), find(contains(sat_list,'2130')), find(contains(sat_list,'2145')),...
                  find(contains(sat_list,'2200')), find(contains(sat_list,'2215')), find(contains(sat_list,'2230')), find(contains(sat_list,'2245')),...
                  find(contains(sat_list,'2300')), find(contains(sat_list,'2315')), find(contains(sat_list,'2330')), find(contains(sat_list,'2345'))    );
    sat_list(KILL) = [];

end
[sa sb] = size(sat_list); clear sb;










for tsat = 1:sa  %9:1:22 %37
 
    % tsat = 8   ;
    % tsat = 19  ;
    % tsat = 16   ;
    % tsat = 11  ;
    
    %%%%% load GOES data
    thisfile = char(sat_list(tsat)) ;
    % ncdisp( thisfile )
    lon = ncread(thisfile,'longitude');
    lat = ncread(thisfile,'latitude');
    CTtemp = ncread(thisfile,'cloud_top_temperature');
    CTht = ncread(thisfile,'cloud_top_height');
    CTeht = ncread(thisfile,'cloud_effective_height');
    refl = ncread(thisfile,'reflectance_vis');
    time = ncread(thisfile,'time');
    basetime = ncread(thisfile,'base_time');

    sat_hr = round( 24*(single(time)) / 86400 , 5 )  

    tind_swath = find(  abs(sat_hr-swath_hr) == min(abs(sat_hr-swath_hr) )   ) ;



    %%% load LASSO cld & make cld top height field
    thiscldfile = char(lassocld_list(tsat))  ;
    %  ncdisp( thiscldfile )
    qcloud = ncread(thiscldfile,'QCLOUD');
    qice= ncread(thiscldfile,'QICE');
    qcloud(qcloud==0) = NaN;
    qice(qice==0) = NaN;
        %qcloud(qcloud<0.00025) = NaN;
    cld_lon = ncread(thiscldfile,'XLONG');
    cld_lat = ncread(thiscldfile,'XLAT');
    cld_hasl = ncread(thiscldfile,'HAMSL');
    cldtime = ncread(thisfile,'time');
    cld_hr = round( 24*(single(cldtime)) / 86400 , 5 )  
    [as df gh] = size(qcloud) ;
    LassoCldTopHt = zeros(as, df);        LassoCldTopHt(:) = NaN;
%     for i = 1:as
%        for j = 1:df
%             ktop = find(  qcloud(i,j,:)  ==  max(qcloud(i,j,:),[],'omitnan') );
%             if( isnan(ktop)==0 )
%                 LassoCldTopHt(i,j) = cld_hasl(ktop(end)); 
%             end
%        end
%     end
    for i = 1:as
       for j = 1:df
            ktop = find(  qcloud(i,j,:)  > 0.0000000 | qice(i,j,:)  > 0.0000000);  
            if( isnan(ktop)==0 )
                ktop = ktop(end);
                LassoCldTopHt(i,j) = cld_hasl(ktop(end)); 
            end
       end
    end



%     % coarsen the lasso cld field to something closer to GOES?
%     dx = 0.1;  %model dx
% 
%     % first pass to kill little ones
%     sample = 0.75 ; %km
%     LassoCldTopHt_coarse = movmean(LassoCldTopHt,        sample/dx, 1);
%     LassoCldTopHt_coarse = movmean(LassoCldTopHt_coarse, sample/dx, 2);
%     % second pass after littles dead
%     sample = 3.0 ; %km
%     LassoCldTopHt_coarse = movmean(LassoCldTopHt_coarse, sample/dx, 1, 'omitnan');
%     LassoCldTopHt_coarse = movmean(LassoCldTopHt_coarse, sample/dx, 2, 'omitnan');
  






%     %%%% plot GOES REFL:
% 
%     %refl(refl<0.05) = NaN;
%     refl(refl<0.1) = NaN;
% 
%     ff = figure
%     ff.Position = [1999,154,442,676];
% 
%     set(gca,'XTick',[])
%     set(gca,'YTick',[])
% 
%     title(   ['REFL GOES time = ', num2str( sat_hr ) ,'hr UTC;   swath time = ', num2str( swath_hr(tind_swath)),'hr UTC']    )
% 
%     ax1 = axes;
%     ax2 = axes;
%     ax3 = axes;
%     ax4 = axes;
%     ax5 = axes;
%     linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');
% 
%     contourf(ax1,lon_terr,lat_terr,terr,[00:250:3000],'LineColor','none')
%     colormap(ax1,demcmap(13))
% 
%     hold on
%     
%     contourf(ax2,lon,lat,refl,20,'LineColor','none') ;
%     colormap(ax2,bone)
%     caxis(ax2,[-0.5 1.1])
% 
%     %contour(ax3,lon,lat,refl,[0.1 0.11],'LineColor','k') ;
% 
%     plot(ax4,-64.7287,-32.1264,'ko','MarkerFaceColor','b','MarkerSize',10)
% 
% %     convobj = double(  permute(Conv_object(:,:,tind_swath,:),[1 4 2 3])   );
% %     convobj(convobj ==0) = NaN;
% %     surf(ax4, Lon, Lat, max(convobj,[],3,'omitnan'),'FaceColor',[0 0 0],'EdgeColor','none' ) ; % [0.7 0.0 0.5],'EdgeColor','none' )
% 
%     wobj = double(  permute(W_object(:,:,tind_swath,:),[1 4 2 3])   );
%     wobj(wobj ==0) = NaN;
%     surf(ax5, Lon, Lat, max(wobj,[],3,'omitnan'),'EdgeColor','r','FaceColor','none' )
% 
%     set(ax2,'Color','None')       %p
%     set(ax2, 'visible', 'off');   %p
% 
%     set(ax3,'Color','None')       %p
%     set(ax3, 'visible', 'off');   %p
% 
%     set(ax4,'Color','None')       %p
%     set(ax4, 'visible', 'off');   %p
% 
%     set(ax5,'Color','None')       %p
%     set(ax5, 'visible', 'off');   %p
% 
%     axis equal
%     axis([-65.2 -64.4 -32.8 -31.4])
% 
% %%%%%%%% image out:
% 
% %saveas(ff,horzcat('/Users/marq789/Downloads/images/Swath_obsgoes_',num2str(tsat),'.png'));
% 
% % outlab = horzcat(imout,'/MPOrigins_AtMCSI_filtLS',num2str(filteroutLS),'.eps');
% % EPSprint = horzcat('print -painters -depsc ',outlab);
% %eval([EPSprint]);















    %%%% plot GOES CTH:

    CTht(CTht<1) = NaN;

    ff = figure
    ff.Position = [1999,154,442,676];

    set(gca,'XTick',[])
    set(gca,'YTick',[])

    title(   ['CldHt GOES time = ', num2str( sat_hr ) ,'hr UTC;   swath time = ', num2str( swath_hr(tind_swath)),'hr UTC']    )

    ax1 = axes;
    ax2 = axes;
    ax3 = axes;
    ax4 = axes;
    ax5 = axes;
    linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');

    contourf(ax1,lon_terr,lat_terr,terr,[00:250:3000],'LineColor','none')
    colormap(ax1,demcmap(13))
    %colormap(ax1,flipud(bone(13)))

    hold on
    
    contourf(ax2,lon,lat,CTht,20,'LineColor','none') ;
    colormap(ax2,bone)
    caxis(ax2,[2 10])

    %contour(ax3,lon,lat,refl,[0.1 0.11],'LineColor','k') ;

    plot(ax4,-64.7287,-32.1264,'ko','MarkerFaceColor','b','MarkerSize',10)

%     convobj = double(  permute(Conv_object(:,:,tind_swath,:),[1 4 2 3])   );
%     convobj(convobj ==0) = NaN;
%     surf(ax4, Lon, Lat, max(convobj,[],3,'omitnan'),'FaceColor',[0 0 0],'EdgeColor','none' ) ; % [0.7 0.0 0.5],'EdgeColor','none' )

    wobj = double(  permute(W_object(:,:,tind_swath,:),[1 4 2 3])   );
    wobj(wobj ==0) = NaN;
    surf(ax5, Lon, Lat, max(wobj,[],3,'omitnan'),'EdgeColor','r','FaceColor','none' )

    set(ax2,'Color','None')       %p
    set(ax2, 'visible', 'off');   %p

    set(ax3,'Color','None')       %p
    set(ax3, 'visible', 'off');   %p

    set(ax4,'Color','None')       %p
    set(ax4, 'visible', 'off');   %p

    set(ax5,'Color','None')       %p
    set(ax5, 'visible', 'off');   %p

    axis equal
    axis([-65.2 -64.4 -32.8 -31.4])

%%%%%%%% image out:

saveas(ff, horzcat('/Users/marq789/Downloads/images/Swath_GOEScth_',RUNS(R,:),'_',num2str(tsat),'.png') );
outlab = horzcat('/Users/marq789/Downloads/images/Swath_GOEScth_',RUNS(R,:),'_',num2str(tsat),'.eps');
EPSprint = horzcat('print -painters -depsc ',outlab);
eval([EPSprint]);




tsat = tsat + 4

    %%%% plot LASSO CTH:  LassoCldTopHt

%     LassoCldTopHt(LassoCldTopHt<1) = NaN;

    ff = figure
    ff.Position = [1999,154,442,676];


    set(gca,'XTick',[])
    set(gca,'YTick',[])

    title(   ['CldHt LASSO time = ', num2str( cld_hr ) ,'hr UTC;   swath time = ', num2str( swath_hr(tind_swath)),'hr UTC']    )

    ax1 = axes;
    ax2 = axes;
    ax3 = axes;
    ax4 = axes;
    ax5 = axes;
    linkaxes([ax1,ax2,ax3,ax4,ax5],'xy');

    contourf(ax1,lon_terr,lat_terr,terr,[00:250:3000],'LineColor','none')
    colormap(ax1,demcmap(13))
    %colormap(ax1,flipud(bone(13)))

    hold on
    
    contourf(ax2,cld_lon,cld_lat,LassoCldTopHt,20,'LineColor','none') ;
    colormap(ax2,bone)
    caxis(ax2,[2000 10000])

    

    %contour(ax3,lon,lat,refl,[0.1 0.11],'LineColor','k') ;

 %   plot(ax4,-64.7287,-32.1264,'ko','MarkerFaceColor','b','MarkerSize',10)

%     convobj = double(  permute(Conv_object(:,:,tind_swath,:),[1 4 2 3])   );
%     convobj(convobj ==0) = NaN;
%     surf(ax4, Lon, Lat, max(convobj,[],3,'omitnan'),'FaceColor',[0 0 0],'EdgeColor','none' ) ; % [0.7 0.0 0.5],'EdgeColor','none' )

    wobj = double(  permute(W_object(:,:,tind_swath,:),[1 4 2 3])   );
    wobj(wobj ==0) = NaN;
    surf(ax5, Lon, Lat, max(wobj,[],3,'omitnan'),'EdgeColor','r','FaceColor','none' )

    set(ax2,'Color','None')       %p
    set(ax2, 'visible', 'off');   %p

    set(ax3,'Color','None')       %p
    set(ax3, 'visible', 'off');   %p

    set(ax4,'Color','None')       %p
    set(ax4, 'visible', 'off');   %p

    set(ax5,'Color','None')       %p
    set(ax5, 'visible', 'off');   %p

    axis equal
    axis([-65.2 -64.4 -32.8 -31.4])

%%%%%%%% image out:

saveas(ff,horzcat('/Users/marq789/Downloads/images/Swath_LASSOcth_',RUNS(R,:),'_',num2str(tsat),'.png'))
outlab = horzcat('/Users/marq789/Downloads/images/Swath_LASSOcth_',RUNS(R,:),'_',num2str(tsat),'.eps');
EPSprint = horzcat('print -painters -depsc ',outlab);
eval([EPSprint]);


end




% for colorbars:

figure
    contourf(cld_lon,cld_lat,LassoCldTopHt,20,'LineColor','none') ;
    colormap(bone)
    caxis([2000 10000])
colorbar


figure
    contourf(lon_terr,lat_terr,terr,[00:250:3000],'LineColor','none')
    colormap(demcmap(13))
colorbar
    axis equal
    axis([-65.2 -64.4 -32.8 -31.4])



















%%%%%%%%%%% just lasso's for the cases I dont have goes downladed


R = 7 ;

%%%%%%  load LASSO swaths:
rootdir =  '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/';
RUNS = ['20181129_gefs09_base_d4';
        '20181204_gefs19_base_d4';
        '20190122_gefs18_base_d4';
        '20190123_gefs18_base_d4';
        '20190125_eda07_base__d4';
        '20190129_gefs11_base_d4';
        '20190208_eda08__base_d4'];

swathmatfilein = [rootdir,'/', RUNS(R,:) ,'/SwathIDthermoOutput_100m_rhcamod_wID_qcldicesnow0.0001_dbz10_run_',RUNS(R,:),'.mat'] ;
load(swathmatfilein,'Terr','Lat','Lon') % ,'terr','W_object','Conv_object','HHMMSS_all','Lat','Lon')
% swathHH = str2num(HHMMSS_all(:,1:2)) ;
% swathMM = str2num(HHMMSS_all(:,3:4)) ;
% swath_hr = swathHH + swathMM/60;
% % ignore the 05,10,20,25,35,40,50,55 min analyses:
% kill = find( swathMM == 5  |  swathMM == 10  |  swathMM == 20  |  swathMM == 25  | ...
%              swathMM == 35  |  swathMM == 40  |  swathMM == 50  |  swathMM == 55 ) ;
% swathHH(kill) = [];
% swathMM(kill) = [];
% swath_hr(kill) = [];
% W_object(:,:,kill,:) = [];
% Conv_object(:,:,kill,:) = [];
% HHMMSS_all(kill,:) = [];


%%%% load LASSO cloud & w fields
lassocldfiles = ls(strcat(rootdir,'/', RUNS(R,:) , '/trimxy_asl_corlasso.cldhaasl.*nc'))
lassocld_list = split(lassocldfiles) ;
lassocld_list(end) = []  ;

[sa sb] = size(lassocld_list) ;   clear sb





thiscldfile = char(lassocld_list(1))   ;
qcloud = ncread(thiscldfile,'QCLOUD')  ;

[as df gh] = size(qcloud) ;
LassoCldTopHt = zeros(as, df,gh);        LassoCldTopHt(:) = NaN;

parfor tsat = 1:sa

    %%% load LASSO cld & make cld top height field
    thiscldfile = char(lassocld_list(tsat))  ;
    %  ncdisp( thiscldfile )

    qcloud = ncread(thiscldfile,'QCLOUD');
    qice= ncread(thiscldfile,'QICE');

    qcloud(qcloud==0) = NaN;
    qice(qice==0) = NaN;

    %qcloud(qcloud<0.00025) = NaN;
    cld_lon = ncread(thiscldfile,'XLONG');
    cld_lat = ncread(thiscldfile,'XLAT');
    cld_hasl = ncread(thiscldfile,'HAMSL');
    cldtime = ncread(thiscldfile,'Time');
    cld_hr = round( 24*(single(cldtime)) / 86400 , 5 )
    [as df gh] = size(qcloud) ;
    LCldTopHt = zeros(as, df);        LCldTopHt(:) = NaN;

    for i = 1:as
        for j = 1:df
            ktop = find(  qcloud(i,j,:)  > 0.0000000 | qice(i,j,:)  > 0.0000000);
            if( isnan(ktop)==0 )
                ktop = ktop(end);
                LCldTopHt(i,j) = cld_hasl(ktop(end));
            end
        end
    end
    
    LassoCldTopHt(:,:,tsat) = LCldTopHt;

end



for tet = 1 : 3 : 67

    ff = figure
    ff.Position = [1999,154,442,676];

    thiscldfile = char(lassocld_list(tet))   ;

    hhmmss = thiscldfile(end-8:end-3)

    set(gca,'XTick',[])
    set(gca,'YTick',[])

    title(   ['  08 Feb CldHt LASSO time = ', hhmmss ,'hr UTC ']    )

    ax1 = axes;
    ax2 = axes;
    linkaxes([ax1,ax2],'xy');

    contourf(ax1,Lon,Lat,Terr,[0:250:3000],'LineColor','none')
    colormap(ax1,demcmap(13))
    %colormap(ax1,flipud(bone(13)))

    hold on

    contourf(ax2,cld_lon,cld_lat,LassoCldTopHt(:,:,tet),20,'LineColor','none') ;
    colormap(ax2,bone)
    caxis(ax2,[2000 10000])

    set(ax2,'Color','None')       %p
    set(ax2, 'visible', 'off');   %p

%     pic = ['/Users/marq789/Downloads/images/Feb08_CldHt_LASSO_', hhmmss ,'UTC.png '];
%     saveas(ff,pic)

end
