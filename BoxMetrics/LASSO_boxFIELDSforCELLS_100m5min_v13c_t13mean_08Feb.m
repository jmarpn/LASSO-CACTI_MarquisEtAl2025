% 100m run processing:
%  W and conv swath detecion in raw (clear + rainy/cloudy) air':
%  W and conv swath detecion in 'clear(ish) air':
%  w-swath: width, depth, peak magnitude per j-grid point
%  w-swath: vertical mass flux, mucape, mucin, mulfc, qvap, theta, thetae 
%  along swath: RH 2.5-5 km asl, ground-rel U,V from bottome of swath (ASL) to to of grid (10km)
%  For non-clear-air swath detection, confines search area for depth to the
%        local LFC [to help mitigate significant cloudy buoyant updrafts from
%        contaminating mesoscale swath width/depth/depth. THough this will
%        make the local LFC a ceiling for swath depth measurements, we can probably assume
%        that if the ascent hits the LFC (which it does a lot postCI) then
%        it's probably less relevant to look at swath depth anymore] % 

%        **** need to finish testing(revisit the conv pattern on the west
%        side of terrain) ****

% v7B -  adds identification for all-cloudy-updrafts (i.e., includes pre-CI); B = parallel version
% v7C -  FIXED BUG in calcualtion of updraft area and non-clear-air-VM from previous v7's. Also, now only keeps low-level updraft profile data within layer relative to LFC 
% v8: this is the 500m v7D but adapted to 100m, including reading in ASL fields from a separate subsetting on cumulus using Bill's code.

% v8c: looks for max(conv_mean) using the conv_agl field rather than asl field (asl in all prev versions)


% v9: separated into 1) swath [ID, width/depth/mag], 2) wID. 
% v9b:   added back in-swath thermo and swath-box metrics to (1)



clear all

% for nersc:


delete(gcp('nocreate'))
pc = parcluster('local')
parpool(pc, 32);
spmd rank = labindex;
    fprintf(1,'Hello from %d\n',rank);
end
disp(' START ')
now





% % % home
% rootdir = '/Volumes/LaCie/LASSO_runs/cacti/lasso/les/';
% %rootdir = '/Users/marq789/Downloads/';

% % nersc
rootdir = '/pscratch/sd/j/jmarquis/cacti/lasso/d4_100m_5min/'

%nersc:
%RUNS = ['22Jan100m/20190122/eda00/morr/les/subset_d4/'];
RUNS = ['08Feb100m/20190208/eda08/base/les/subset_d4/'];

% %%home tester:
% rootdir = ['/Volumes/LaCie/LASSO_runs/cacti/lasso/les/'];
% RUNS = ['20190129/gefs11/base/d4/'];



[ar br] = size( RUNS ) ;  clear br


% % loop through ensemble member runs:
% for R = 1 : ar

%    clearvars -except rootdir RUNS ar R

R = 1

%wrfdir = strcat(rootdir,RUNS(R,:))  ;
wrfdir = strcat(rootdir,RUNS(R,:) ,'/f15min/')  ;

imout = strcat(wrfdir,'/testimages/');
mkdir(imout)

%wrfdir = strcat('/Volumes/LaCie/LASSO_runs/', RUN ,'/') ;
%matoutsave = '/Volumes/LaCie/LASSO_runs/29Jan500m/SwathObjectOutput.mat' ;


%somewhat hardcoded lable for output
wherethe2at = find(wrfdir == '2');
wherethe4at = find(wrfdir == '4');
runlab = wrfdir(wherethe2at(1):wherethe4at(end))   ;
runlab(runlab == '/') = '_'  ;
matoutsave = strcat( wrfdir, 'UVboxesterrainfollowing_Output_100m_qcldicesnow0.0001_dbz10_v13c_run_',runlab,'.mat' )  ;

%CLD_THRESH = 0.0001; DBZ_THRESH = 30;


% want to plot? 1 = yes
plotme = 0;

%%% LASSo grid deltax,y (km)
dx = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% START DOING STUFF NOW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' processing data in : ')
disp( wrfdir )
disp('   ')

%{
%post-processed LASSO AGL cloud data files:
cldlist = ls( horzcat(wrfdir,'/trimxy_corlasso.cldhaagl.*nc') );
cldlist = split(cldlist);
[sa sb] = size(cldlist); clear sb;

%post-processed LASSO lifted parcel metric files:
parcellist = ls( horzcat(wrfdir,'/trimxy_corlasso.liftphagl*nc') );
parcellist = split(parcellist);
[sa sb] = size(parcellist); clear sb;
%}

%post-processed LASSO met files:
metlist = ls( horzcat(wrfdir,'/trimxy_corlasso.methagl*nc') );
metlist = split(metlist);
[sa sb] = size(metlist); clear sb;

%%%%%%%%%%%%% diagnostic tool
%sa = 4;   %tester: reduction of the full sample set
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I looked at the terrain field and chose a select subdomain, so trim away some domain fat early:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  % AGL SUBDOMAIN TO LOAD, dont use these for science
%     i1 = 1;   i2 = 900;    %first/last i index to keep
%     j1 = 1;   j2 = 1901;   %first/last j index to keep
%    i1 = 200;   i2 = 700;    %first/last i index to keep
i1 = 100;   i2 = 600;    %first/last i index to keep
j1 = 201;   j2 = 1801;   %first/last j index to keep

di = i2 - i1;
dj = j2 - j1;

metf = char(metlist(1)) ;
start = [i1 j1 1];
count = [di dj 1];
Terr = ncread(metf,'HGT',start,count);

start = [i1 j1];
count = [di dj];
Lon = ncread(metf,'XLONG',start,count);
Lat = ncread(metf,'XLAT',start,count);
Ht = ncread(metf,'HAGL');


dum = ncread(metf,'UA') ;
[aa bb cc] = size(dum) ;  clear dum
dummy = zeros(di,dj,cc,sa-1) ;
dummy(:) = NaN;


disp('CC')

%     Temp_all = dummy;
%     THETA_all = dummy;
%     THETAE_all = dummy;
%     QVAP_all = dummy;
%     RH_all = dummy;
%     P_all = dummy;

U_all = dummy;
V_all = dummy;
%W_all = dummy;

%newFRN
%     QVAP_all = dummy;
%     THETA_all = dummy;

%     DBZ_all = dummy;
%     Qcld_all = dummy;

clear dummy

MULFC_all = zeros(di,dj,sa-1) ;      MULFC_all(:) = NaN;
MUCIN_all = zeros(di,dj,sa-1) ;      MUCIN_all(:) = NaN;
MUCAPE_all = zeros(di,dj,sa-1) ;     MUCAPE_all(:) = NaN;



disp('   ')
disp( ' loading post-processed AGL LASSO files for: '    )
disp(runlab)
disp('   ')
tic

parfor tt = 1 : sa-1

    disp('   ')
    disp('loading 3d vars for time')
    disp( num2str(tt) )

    metfile = char(metlist(tt));   % ncdisp(metfile)
    %cldfile = char(cldlist(tt));
    %parcelfile = char(parcellist(tt));

    start3d = [i1 j1 1 1];
    count3d = [di dj Inf 1];
    start2d = [i1 j1 1];
    count2d = [di dj 1];

    U_all(:,:,:,tt) = ncread(metfile,'UA',start3d,count3d);
    V_all(:,:,:,tt) = ncread(metfile,'VA',start3d,count3d);
    %W_all(:,:,:,tt) = ncread(metfile,'WA',start3d,count3d);

    %newFRN
    %THETA_all(:,:,:,tt) = ncread(metfile,'THETA',start3d,count3d);
    %QVAP_all(:,:,:,tt) = ncread(metfile,'QVAPOR',start3d,count3d);

    %        Temp_all(:,:,:,tt) = ncread(metfile,'TEMPERATURE',start3d,count3d);       %K
    %         THETA_all(:,:,:,tt) = ncread(metfile,'THETA',start3d,count3d);
    %         THETAE_all(:,:,:,tt) = ncread(metfile,'THETA_E',start3d,count3d);
    %        QVAP_all(:,:,:,tt) = ncread(metfile,'QVAPOR',start3d,count3d);
    %         RH_all(:,:,:,tt) = ncread(metfile,'RH',start3d,count3d);
    %        P_all(:,:,:,tt) = ncread(metfile,'PRESSURE',start3d,count3d);    %mb
    %
    %         DBZ_all(:,:,:,tt) = ncread(cldfile, 'REFL_10CM',start3d,count3d);
    %         Qcld_all(:,:,:,tt) = ncread(cldfile, 'QCLOUD',start3d,count3d);

    %MULFC_all(:,:,tt) = ncread(parcelfile, 'MULFC',start2d,count2d);      % according to bill, this is in ht ASL frame
    %MUCIN_all(:,:,tt) = ncread(parcelfile, 'MUCIN',start2d,count2d);
    %MUCAPE_all(:,:,tt) = ncread(parcelfile, 'MUCAPE',start2d,count2d);


    %   ncdisp(parcelfile)

end
toc
disp('   ')

disp('DD')


THETAv_all  = THETA_all .* (1 + 0.61 * QVAP_all);   %newFRN

% semi-hardcoded time/date labels for LASSO data:
YYMMDD_all = [];
HHMMSS_all = [];
for tt = 1:1:sa-1
    metfile = char(metlist(tt)) ;
    YYMMDD_all = vertcat( YYMMDD_all, metfile(end-17:end-10) );
    HHMMSS_all = vertcat( HHMMSS_all, metfile(end-8:end-3) );
end


%{
% use this later for low-level updraft ID'ing
MULFC_agl_all = MULFC_all;   MULFC_agl_all(:) = NaN;   %  MULFC in AGL frame
for t = 1:sa-1
    MULFC_agl_all(:,:,t) = MULFC_all(:,:,t) - Terr ;
end
%}

% move on to calculate desireable variables
CLD_THRESH = 0.0001; DBZ_THRESH = 10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%     % define cloudy and clear-air updrafts:
%     Wcld_all = W_all; Wcld_all(:,:,:,:) = NaN;
%     cloudy = find(Qcld_all > CLD_THRESH | DBZ_all > DBZ_THRESH);
%     Wcld_all(cloudy) = W_all(cloudy);
%
%     Wca_all = W_all;
%     Wca_all(cloudy) = NaN;  % Wca "clear air" is now background (not cloudy) updraft
%
%     %%%% clear-air qvapor:
%     QVca_all = QVAP_all;
%     QVca_all(cloudy) = NaN;  % QVAPca "clear air" is now background (not cloudy)
%
%     %%%% clear-air theta:
%     THca_all = THETA_all;
%     THca_all(cloudy) = NaN;  % QVAPca "clear air" is now background (not cloudy)
%
%     %%%% clear-air theta-e:
%     THEca_all = THETAE_all;
%     THEca_all(cloudy) = NaN;  % QVAPca "clear air" is now background (not cloudy)
%
%     RHca_all = RH_all;
%     RHca_all(cloudy) = NaN;  % QVAPca "clear air" is now background (not cloudy)


[Xs Ys Zs Ts] = size(U_all(:,:,:,:));


%clear cloudy %QVAP_all THETA_all THETAE_all RH_all


%     %%%%%%%%%%%
%     %%%%%  vert mass flux calculations:
%     %%%%%%%%%%%
%
%     P_all = (P_all*100);            %pa
%     tmp = (QVAP_all ./ 0.622);
%     clear QVAP_all
%
%     vappres = (tmp .* P_all) ./ ( 1 + tmp);       %pa
%     rhom = (    P_all ./ (287 .* Temp_all)  ) .* (   1 -  (vappres./P_all).* ( 1 - 0.622 )    ) ;
%     clear P_all tmp vappres Temp_all
%
%     %VMca = Wca_all .* rhom ;
%     VM = W_all .* rhom ;
%     clear rhom


disp('EE')

%calculate height ASL field
Ht_asl= U_all(:,:,:,1); Ht_asl(:) = NaN;
for i = 1:Xs
    for j = 1:Ys
        Ht_asl(i,j,:) = Ht  + Terr(i,j);
    end
end


%3d-ize height agl field for later conveneience:
Ht_agl = U_all(:,:,:,1); Ht_agl(:) = NaN;
for i = 1:Xs
    for j = 1:Ys
        Ht_agl(i,j,:) = Ht;
    end
end
clear Ht_agl


%%%%%% define the x,y arrays:
X = Lon; X(:) = 0; Y = X ;
xx = dx*[1:Xs]; yy = dx*[1:Ys] ;

%3D-fy horiz grid for later conveneience:
X = U_all(:,:,:,1); X(:) = NaN; Y = X;
for k = 1:Zs
    for j = 1:Ys
        X(:,j,k) = xx ;
    end
    for i = 1:Xs
        Y(i,:,k) = yy ;
    end
end

%regrid X,Y to have 0.0 be at radar position:
AMF_lat = -32.2160 ;  AMF_lon = -64.7284;
dlat = abs(Lat - AMF_lat); dlon = abs(Lon - AMF_lon);
dtot = abs(dlat + dlon);
[AMFi AMFj] = find( dtot == min(min(dtot))    );
AMFx = X(AMFi,AMFj,1);      AMFy = Y(AMFi,AMFj,1);
%corrects to CSAPR2 location
X = X - AMFx;  Y = Y - AMFy;
%then further corrects to typical RELAMPAGO DOW7 position (for dual-Doppler comparisons)
X = X - 1.97;
Y = Y - 16.53;



% disp(' calc AGL convergence')
% 
% % make convergence 'mask' that encompasses all times:
% [dudy dudx dudz dudt] = gradient(U_all,dx*1000,dx*1000,1,1);
% [dvdy dvdx dvdz dvdt] = gradient(V_all,dx*1000,dx*1000,1,1);
% clear U_all V_all
% 
% 
% clear dudy dudz dvdx dvdz
% Conv = -( dudx +  dvdy );
% clear dudx dvdy dvdt dudt dlat dlon dtot

% up to here, the full AGL variables are loaded and have derivatives calculated



% figure; contourf(Conv(:,:,5,5)',20,'LineColor','none'); hold on; contour(Terr'/100000,[1500:200:2500]/100000,'k')





%%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%%%%  ASL grid arrays
%%%%%%%%%%%%
%%%%%%%%%%%%


%%%% do I really need U_all,VCS V_all, V_VCS anymore?
%%%% check hereafter which full vs ca QV, THETA THETAe I really need anymore?
%%%% AGLinASL? Qcld DBZ?

disp(' Reading in subset ASL data ')

%%%% parallel version --- seems to have memory crashes %%%%%
%%%%%%%%%%%% generate height ASL versions by interpolating


%define your own ASL grid
%newz_asl = [0,500:25:4500,4750:250:10000];  %500m original set
newz_asl = [10,500:100:4500,4750,5000:1000:10000] ;

%new x,y,z grids in ASL frame:
[regxasl regyasl regzasl] = meshgrid( X(:,1,1), Y(1,:,1)', newz_asl');
regxasl = permute(regxasl,[2 1 3]); regyasl = permute(regyasl,[2 1 3]); regzasl = permute(regzasl,[2 1 3]);
[xasl yasl zasl] = meshgrid(X(:,1,1), Y(1,:,1)', Ht);  %Ht is just a dummy var, it will be tossed out
zasl = Ht_asl;
xasl = permute(xasl,[2 1 3]); yasl = permute(yasl,[2 1 3]); %zasl = permute(zasl,[2 1 3]);

clear Ht_asl regyasl

% make a ht agl field in asl coord
disp(' calc aglinasl  ')
[e1 e2 e3] = size(regzasl)
AGLinASL = regzasl; AGLinASL(:) = NaN;
for i = 1:e1
    for j = 1:e2
        AGLinASL(i,j,:) = regzasl(i,j,:) - Terr(i,j);
    end
end
AGLinASL(AGLinASL<0) = NaN;

%    save( strcat(wrfdir,'test_aglinasl.mat'), 'AGLinASL', 'Terr', 'regzasl', '-v7.3' )
%    disp('tmp done')

%{
%post-processed LASSO AGL cloud data files:
aslcldlist = ls( horzcat(wrfdir,'/trimxy_asl_corlasso.cldhaasl.*nc') );
aslcldlist = split(aslcldlist);
[sa sb] = size(aslcldlist); clear sb;

%post-processed LASSO lifted parcel metric files:
aslparcellist = ls( horzcat(wrfdir,'/trimxy_asl_corlasso.liftphasl*nc') );
aslparcellist = split(aslparcellist);
[sa sb] = size(aslparcellist); clear sb;
%}

%post-processed LASSO met files:
aslmetlist = ls( horzcat(wrfdir,'/trimxy_asl_corlasso.methasl*nc') );
aslmetlist = split(aslmetlist);
[sa sb] = size(aslmetlist); clear sb;


%sa = 4; % Ts = sa

%%%%%%%%
%%%   readin ASL fields from BIll subsets
%%%%%%%%

%     metf = char(aslmetlist(1)) ;
%
%     dum = ncread(metf,'UA') ;
%     [aa bb cc] = size(dum) ;  clear dum
%     dummy = zeros(di,dj,cc,sa-1) ;
%     dummy(:) = NaN;
%
%     Temp_all = dummy;
%     THETA_all = dummy;
%     THETAE_all = dummy;
%     QVAP_all = dummy;
%     RH_all = dummy;
%     P_all = dummy;
%     U_all = dummy;
%     V_all = dummy;
%     W_all = dummy;
%     DBZ_all = dummy;
%     Qcld_all = dummy;


disp('seeding full ASL vars')
%Conv_ASL 	= zeros( (i2-i1), DJ, length(newz_asl),Ts);
U_ASL2 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
V_ASL2 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%W_ASL2 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%Wca_ASL 	= zeros( (i2-i1), DJ, length(newz_asl),Ts);
%Qcld_ASL 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%DBZ_ASL 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%QVca_ASL 	= zeros( (i2-i1), DJ, length(newz_asl),Ts);
% THca_ASL 	= zeros( (i2-i1), DJ, length(newz_asl),Ts);
% THEca_ASL 	= zeros( (i2-i1), DJ, length(newz_asl),Ts);
%RHca_ASL 	= zeros( (i2-i1), (i2-i1), length(newz_asl),Ts);
%AGLinASL 	= zeros( (i2-i1), DJ, length(newz_asl));
%QVAP_ASL 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%THETA_ASL 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%THETAE_ASL 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%RH_ASL 	    = zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%VMca_ASL 	= zeros( (i2-i1), DJ, length(newz_asl),Ts);
%VM_ASL      = zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%Temp_ASL    = zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
%P_ASL       = zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);

%newFRN
QVAP_ASL2 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);
THETA_ASL2 	= zeros( (i2-i1), (j2-j1), length(newz_asl),Ts);

%disp('size of seeded asl var')
%size(U_ASL)

disp('   ')
disp( ' loading post-processed ASL LASSO files for: '    )
disp(runlab)
disp('   ')
tic

parfor tt = 1 : sa-1

    disp('   ')
    disp('loading 3d asl vars for time')
    disp( num2str(tt) )

    aslmetfile = char(aslmetlist(tt));   % ncdisp(metfile)
    %aslcldfile = char(aslcldlist(tt));
    %aslparcelfile = char(aslparcellist(tt));

    start3d = [i1 j1 1 1];
    count3d = [di dj Inf 1];
    start2d = [i1 j1 1];
    count2d = [di dj 1];

    U_ASL2(:,:,:,tt) = ncread(aslmetfile,'UA',start3d,count3d);
    V_ASL2(:,:,:,tt) = ncread(aslmetfile,'VA',start3d,count3d);
 %   W_ASL2(:,:,:,tt) = ncread(aslmetfile,'WA',start3d,count3d);

    %newFRN
    THETA_ASL2(:,:,:,tt) = ncread(aslmetfile,'THETA',start3d,count3d);
    QVAP_ASL2(:,:,:,tt) = ncread(aslmetfile,'QVAPOR',start3d,count3d);

 %   Temp_ASL(:,:,:,tt) = ncread(aslmetfile,'TEMPERATURE',start3d,count3d);       %K
 %   THETA_ASL(:,:,:,tt) = ncread(aslmetfile,'THETA',start3d,count3d);
 %   THETAE_ASL(:,:,:,tt) = ncread(aslmetfile,'THETA_E',start3d,count3d);
 %   QVAP_ASL(:,:,:,tt) = ncread(aslmetfile,'QVAPOR',start3d,count3d);
 %   RH_ASL(:,:,:,tt) = ncread(aslmetfile,'RH',start3d,count3d);
    %
    %DBZ_ASL(:,:,:,tt) = ncread(aslcldfile, 'REFL_10CM',start3d,count3d);
    %Qcld_ASL(:,:,:,tt) = ncread(aslcldfile, 'QCLOUD',start3d,count3d);
 %   P_ASL(:,:,:,tt) = ncread(aslmetfile,'PRESSURE',start3d,count3d);    %mb
    %Qice_ASL(:,:,:,tt) = ncread(aslcldfile, 'QICE',start3d,count3d);
    %   ncdisp(parcelfile)

end
toc
disp('   ')

disp(' done loading ASLs')


%{
%calculate RHca
RHca_ASL = RH_ASL;
clear RH_ASL

cloudy = find( (Qcld_ASL + Qice_ASL ) > CLD_THRESH   |   DBZ_ASL > DBZ_THRESH   |   RHca_ASL > 90.0  );
RHca_ASL(cloudy) = NaN;
clear Qcld_ASL DBZ_ASL cloudy
%}


% disp(' calc ASL convergence')
% 
% % make convergence 'mask' that encompasses all times:
% [dudy dudx dudz dudt] = gradient(U_ASL2,dx*1000,dx*1000,1,1);
% [dvdy dvdx dvdz dvdt] = gradient(V_ASL2,dx*1000,dx*1000,1,1);
% clear dudy dudz dvdx dvdz dvdt dudt
% 
% Conv_ASL2 = -( dudx +  dvdy );
% clear dudx dvdy
% disp( strcat(' size  conv_asl ') )
% size(Conv_ASL2)

%newFRN
THETAv_ASL2 = THETA_ASL2 .* (1 + 0.61 * QVAP_ASL2);

% %%     %%%%%%%%%%%
% %%     %%%%%  vert mass flux calculations:
% %%     %%%%%%%%%%%
%P_ASL = (P_ASL*100);            %pa
%tmp = (QVAP_ASL ./ 0.622);
%%clear QVAP_ASL

%vappres = (tmp .* P_ASL) ./ ( 1 + tmp);       %pa
%rhom = (    P_ASL ./ (287 .* Temp_ASL)  ) .* (   1 -  (vappres./P_ASL).* ( 1 - 0.622 )    ) ;
%clear P_ASL tmp vappres Temp_ASL

%VMca = Wca_all .* rhom ;
%VM_ASL = W_ASL2 .* rhom ;
%clear P_ASL tmp vappres rhom  Temp_ASL




%black out the stuff below the terrain:
%Conv_ASL2 = Conv_ASL ;
%Conv_ASL2(find(regzasl < Terr)) = -999;
%clear Conv_ASL

%U_ASL2 = U_ASL ;
U_ASL2(find(regzasl < Terr)) = -999;
%clear U_ASL

%V_ASL2 = V_ASL ;
V_ASL2(find(regzasl < Terr)) = -999;
%clear V_ASL

%%W_ASL2 = W_ASL ;
%W_ASL2(find(regzasl < Terr)) = -999;
%%clear W_ASL

%     Wca_ASL2 = Wca_ASL ;
%     Wca_ASL2(find(regzasl < Terr)) = -999;
%     clear Wca_ASL

%     Qcld_ASL2 = Qcld_ASL ;
%     Qcld_ASL2(find(regzasl < Terr)) = -999;
%     clear Qcld_ASL
%
%     DBZ_ASL2 = DBZ_ASL ;
%     DBZ_ASL2(find(regzasl < Terr)) = -999;
%     clear DBZ_ASL

%     QVca_ASL2 = QVca_ASL ;
%     QVca_ASL2(find(regzasl < Terr)) = -999;
%     clear QVca_ASL
%
%     THca_ASL2 = THca_ASL ;
%     THca_ASL2(find(regzasl < Terr)) = -999;
%     clear THca_ASL
%
%     THEca_ASL2 = THEca_ASL ;
%     THEca_ASL2(find(regzasl < Terr)) = -999;
%     clear THEca_ASL

%QVAP_ASL2 = QVAP_ASL;
QVAP_ASL2(find(regzasl < Terr)) = -999;  %newFRN
%clear QVAP_ASL

THETAv_ASL2(find(regzasl < Terr)) = -999;  %newFRN

%THETA_ASL2 = THETA_ASL;
THETA_ASL2(find(regzasl < Terr)) = -999;  %newFRN

%{
QVAP_ASL2 = QVAP_ASL;
QVAP_ASL2(find(regzasl < Terr)) = -999;
clear QVAP_ASL

THETA_ASL2 = THETA_ASL;
THETA_ASL2(find(regzasl < Terr)) = -999;
clear THETA_ASL

THETAE_ASL2 = THETAE_ASL;
THETAE_ASL2(find(regzasl < Terr)) = -999;
clear THETAE_ASL
%}
%
%     RH_ASL2 = RH_ASL;
%     RH_ASL2(find(regzasl < Terr)) = -999;
%     clear RH_ASL

%RHca_ASL2 = RHca_ASL;
%RHca_ASL2(find(regzasl < Terr)) = -999;
%clear RHca_ASL

%
%     VMca_ASL2 = VMca_ASL;
%     VMca_ASL2(find(regzasl < Terr)) = -999;
%     clear VMca_ASL

%{
VM_ASL2 = VM_ASL;
VM_ASL2(find(regzasl < Terr)) = -999;
clear VM_ASL
%}

[Xa Ya Za Ta] = size(U_ASL2);

disp(' permuting ASL into VCS frame ')

%3D cross sections on ASL frame:
A_Ht_asl_vcs = permute(regzasl,[1 3 2]) ;
A_AGL_vcs = permute(AGLinASL,[1 3 2]) ;

A_U_vcs = permute(U_ASL2,[1 3 4 2]) ;
clear U_ASL2
A_V_vcs = permute(V_ASL2,[1 3 4 2]) ;
clear V_ASL2
%A_W_vcs = permute(W_ASL2,[1 3 4 2]) ;
%clear W_ASL2

%     A_Wca_vcs = permute(Wca_ASL2,[1 3 4 2]) ;
%     clear Wca_ASL2
% A_Conv_vcs = permute(Conv_ASL2,[1 3 4 2]) ;
% clear Conv_ASL2
%     A_Qcld_vcs = permute(Qcld_ASL2,[1 3 4 2]) ;
%     clear Qcld_ASL2
%     A_DBZ_vcs = permute(DBZ_ASL2,[1 3 4 2]) ;
%     clear DBZ_ASL2
%     A_RH_vcs = permute(RH_ASL2,[1 3 4 2]) ;
%     clear RH_ASL2

%keep these _ASL2 fields:
%    A_QVca_vcs = permute(QVca_ASL2,[1 3 4 2]) ;
%     A_THca_vcs = permute(THca_ASL2,[1 3 4 2]) ;
%     A_THEca_vcs = permute(THEca_ASL2,[1 3 4 2]) ;
%A_RHca_vcs = permute(RHca_ASL2,[1 3 4 2]) ;
%clear RHca_ASL2

%newFRN
A_THETAv_vcs = permute(THETAv_ASL2,[1 3 4 2]) ;

%     A_QVAP_vcs = permute(QVAP_ASL2,[1 3 4 2]) ;
%     A_THETA_vcs = permute(THETA_ASL2,[1 3 4 2]) ;
%     A_THETAE_vcs = permute(THETAE_ASL2,[1 3 4 2]) ;
%     A_RH_vcs = permute(RH_ASL2,[1 3 4 2]) ;
%    A_VMca_vcs = permute(VMca_ASL2,[1 3 4 2]) ;
% A_VM_vcs = permute(VM_ASL2,[1 3 4 2]) ;
% A_X_vcs = permute(regxasl,[1 3 2]) ;
clear regcasl

[Xvcsa Zvcsa Tvcsa Yvcsa] = size(A_U_vcs)  ;



% %%%%%%%%%%%%%%%%%%%%
% %%%i-j smooth the conv_asl field
% %%%%%%%%%%%%%%%%%%%%
% 
% % %   preserve VCS fields
% Conv_unsmoothed = Conv;
% A_Conv_vcs_unsmoothed = A_Conv_vcs;
% A_W_vcs_unsmoothed = A_W_vcs;  %preserved
% 
% %     % BASE
% %     A_Conv_vcs = A_Conv_vcs_unsmoothed;
% %     A_W_vcs = A_W_vcs_unsmoothed;
% %     Conv = Conv_unsmoothed;
% 
% 
% %     % smoothers:
% %1st pass
% A_Conv_vcs = movmean(A_Conv_vcs_unsmoothed,11,4,'omitnan');   %11 = 1-km window, 5 = 400m window
% A_Conv_vcs = movmean(A_Conv_vcs,11,1,'omitnan');
% clear A_Conv_vcs_unsmoothed
% %2pass
% %A_Conv_vcs = movmean(A_Conv_vcs,11,4,'omitnan');
% %A_Conv_vcs = movmean(A_Conv_vcs,11,1,'omitnan');
% %3pass
% %A_Conv_vcs = movmean(A_Conv_vcs,11,4,'omitnan');
% %A_Conv_vcs = movmean(A_Conv_vcs,11,1,'omitnan');
% %4pass
% %A_Conv_vcs = movmean(A_Conv_vcs,11,4,'omitnan');
% %A_Conv_vcs = movmean(A_Conv_vcs,11,1,'omitnan');
% 
% 
% %1st pass
% A_W_vcs = movmean(A_W_vcs_unsmoothed,11,4,'omitnan');   %11 = 1-km window, 5 = 400m window
% A_W_vcs = movmean(A_W_vcs,11,1,'omitnan');
% %2pass
% %A_W_vcs = movmean(A_W_vcs,11,4,'omitnan');
% %A_W_vcs = movmean(A_W_vcs,11,1,'omitnan');
% %3pass
% %A_W_vcs = movmean(A_W_vcs,11,4,'omitnan');
% %A_W_vcs = movmean(A_W_vcs,11,1,'omitnan');
% %4pass
% %A_W_vcs = movmean(A_W_vcs,11,4,'omitnan');
% %A_W_vcs = movmean(A_W_vcs,11,1,'omitnan');
% 
% 
% %1st pass
% Conv = movmean(Conv_unsmoothed,11,1,'omitnan');   %11 = 1-km window
% clear Conv_unsmoothed
% Conv = movmean(Conv,11,2,'omitnan');
% %2pass
% %Conv = movmean(Conv,11,4,'omitnan');
% %Conv = movmean(Conv,11,1,'omitnan');
% %3pass
% %Conv = movmean(Conv,11,4,'omitnan');
% %Conv = movmean(Conv,11,1,'omitnan');
% %4pass
% %Conv = movmean(Conv,11,4,'omitnan');
% %Conv = movmean(Conv,11,1,'omitnan');




% 
% %%%%%%%%%%%%%%%%%%%%%%%   Window means in VCSs:
% 
% disp('    ')
% disp(' mean magic on conv field to find the mean position of the swath  ')
% disp('    ')
% tic
% 
% %     %1 option for FOR AGL MEAN DECTTION: all time conv
% %     Convagl_alltmean = zeros(Xvcsa,38,Yvcsa,Tvcsa);
% %     cagltmp = permute(Conv,[1 3 4 2]) ;
% %     alltmean = mean( cagltmp ,3,'omitnan');      % all time mean
% %     for t = 1:Tvcsa
% %         Convagl_alltmean(:,:,:,t) = alltmean;
% %     end
% %     Convagl_alltmean = permute(Convagl_alltmean,[1 2 4 3]);
% %     clear alltmean
% %     Convagl_alltmean = movmean(Convagl_alltmean,11,1,'omitnan');
% %     Convagl_alltmean = movmean(Convagl_alltmean,11,4,'omitnan');
% %
% %
% %
% %     %all time conv asl frame
% %     Conv_alltmean = zeros(Xvcsa,Zvcsa,Yvcsa,Tvcsa);
% %     alltmean = mean(A_Conv_vcs,3,'omitnan');      % all time mean
% %     for t = 1:Tvcsa
% %         Conv_alltmean(:,:,:,t) = alltmean;
% %     end
% %     Conv_alltmean = permute(Conv_alltmean,[1 2 4 3]);
% %     clear alltmean
% 
% 
% 
% disp(' time running mean conv fields ')
% 
% %     % running time means:
% %      Conv_5tmean = movmean(A_Conv_vcs,5,3,'omitnan');   % running time means (window = 5)  @ f:15, t5 = 1 hr window
% %      Conv_9tmean = movmean(A_Conv_vcs,9,3,'omitnan');   % running time means (window = 9) @ f:15, t9 = 2 hr window
% Conv_mean = movmean(A_Conv_vcs,13,3,'omitnan');   % running time means (window = 13) @ f:15, t13  = 3 hr window   @5minf 13 = 1 hr
% %Conv_25tmean = movmean(A_Conv_vcs,25,3,'omitnan');   %  @5minf 25 = 2 hr
% 
% 
% 
% %cagltmp = permute(Conv,[1 3 4 2]) ;
% Convagl_mean = movmean(permute(Conv,[1 3 4 2]),13,3,'omitnan');   % running time means (window = 13) @ f:15, t13  = 3 hr window   @5minf 13 = 1 hr
% %Convagl_25tmean = movmean(Convvcs,25,3,'omitnan');   %  @5minf 25 = 2 hr
% 
% clear Conv
% 
% %     % spatial smoother:
% %     Conv_mean_11imean_11jmean = movmean(Conv_alltmean,11,4,'omitnan');   %11 = 1-km window
% %     Conv_mean_11imean_11jmean = movmean(Conv_mean_11imean_11jmean,11,1,'omitnan');
% 
% % %2d smoothing:
% % Conv_mean_ismooth = movmean(Conv_5tmean,5,1,'omitnan') ;
% % Conv_mean_2dsmooth = movmean(Conv_mean_ismooth,5,2,'omitnan') ;
% 
% %Conv_mean = Conv_mean_9jmean_9tmean;   %running along-y mean
% % A_U_mean = mean(A_U_vcs,3,'omitnan');
% 
% 
% 
% 
% 
% %{
%     %try time-weighted mean:
%     
%     twt_kernel = zeros(Tvcsa);  twt_kernel = twt_kernel(:,1);
%     % manual weight toward first time
%     nt = Tvcsa;
%     for t = 1:1:Tvcsa
%         twt_kernel(t) = nt;
%         nt = nt-1;
%     end
%     %normalize t-weight kernel
%     twtF_kernel = twt_kernel/sum(twt_kernel);
% 
%     % manual weight toward last time
%     for t = 1:1:Tvcsa
%         twt_kernel(t) = t;
%     end
%     %normalize t-weight kernel
%     twtL_kernel = twt_kernel/sum(twt_kernel);
% 
%     % manual weight toward middle time, linear d(weight)/dt
%     for t = 1:1:floor(Tvcsa/2)
%         twt_kernel(t) = t;
%     end
%     nt = floor(Tvcsa/2)+1;
%     for t = floor(Tvcsa/2)+1:Tvcsa
%         twt_kernel(t) = nt;
%         nt = nt-1;
%     end
%     %normalize t-weight kernel
%     twtM_kernel = twt_kernel/sum(twt_kernel);
% 
% 
% 
%     % pick your kernel:
%     twt_kernel = twtM_kernel;
% 
%     %apply kernel-weighted conv t-mean [AGL]
%     conv_tweight = zeros(Xvcsa,Yvcsa,38,Tvcsa);
%     for t = 1:Tvcsa
%         conv_tweight(:,:,:,t) = twt_kernel(t) .* Conv(:,:,:,t) ;
%     end
%     conv_tweight = sum(conv_tweight,4,'omitnan');
%     %build it out for all times
%     Conv_tweight = [];
%     for t = 1:Tvcsa
%         Conv_tweight = cat(4,Conv_tweight,conv_tweight);
%     end
%     Convagl_tweight = permute(Conv_tweight,[1 3 4 2]);
%     clear Conv_tweight 
% 
% 
%     %apply kernel-weighted conv t-mean [ASL]
%     conv_tweight = zeros(Xvcsa,Yvcsa,Zvcsa,Tvcsa);
%     casl = permute(A_Conv_vcs,[1 4 2 3]);
%     for t = 1:Tvcsa
%         conv_tweight(:,:,:,t) = twt_kernel(t) .* casl(:,:,:,t) ;
%     end
%     conv_tweight = sum(conv_tweight,4,'omitnan');
%     %build it out for all times
%     Conv_tweight = [];
%     for t = 1:Tvcsa
%         Conv_tweight = cat(4,Conv_tweight,conv_tweight);
%     end
%     Convasl_tweight = permute(Conv_tweight,[1 3 4 2]);
%     clear Conv_tweight casl
% %}
% 
% 
% 
% toc
% disp('   ')
% 
% % the field you pick to use in swath identification
% %Conv_mean = Conv_alltmean;
% %Convagl_mean = Convagl_alltmean;
% 
% %Conv_mean = Conv_9tmean;
% 
% %Conv_mean = Convasl_tweight;
% %Convagl_mean = Convagl_tweight;
% 
% %     Conv_mean = Conv_13tmean;
% %     Convagl_mean = Convagl_13tmean;
% 
% disp(' conv mean sizes - asl then agl')
% size(Conv_mean)
% size(Convagl_mean)
% 
% 
% %FOR AGL MEAN DECTTION: kill agl where terrain < 1000m:
% parfor i = 1:Xvcsa
%     for j = 1:Yvcsa
%         if(Terr(i,j) < 500.)
%             Convagl_mean(i,:,:,j) = NaN;
%         end
%     end
% end
% 
% 
% %%%%%%%% manual domain restricting:
% %blanket thresholding:
% %Conv_mean(find(Conv_mean < 0.0001)) = NaN;
% 
% %need to add additional smothering (kill stuff "far" from ridgeline):
% Conv_mean(1:60,:,:,:) = nan;
% Conv_mean(350:end,:,:,:) = nan;
% %  kill whats below Z meters asl
% Conv_mean(:,1:9,:,:) = nan;   %5 = 800m;   7 = 1km
% 
% 
% 
% 
% % find the index of maxima in conv_mean ASL field that is very near the ground and is the largest in targeted subdomain:
% 
% disp('    ')
% disp(' locating the time mean-ed conv max near the terrain  ')
% disp('    ')
% tic
% 
% 
% %temp conv field that is just near the surface:
% dummy = permute(Conv_mean,[1 2 4 3]);
% %abovegrnd =  find( A_AGL_vcs  > 200);
% abovegrnd =  find( A_AGL_vcs  > 300);  %raise it a little from d3 version because of d4's 100m asl DZ
% Conv_mean_nearsfc = dummy; Conv_mean_nearsfc(:) = NaN;
% for t = 1:Tvcsa
%     Cmn = dummy(:,:,:,t);
%     Cmn(abovegrnd) = nan;
%     Conv_mean_nearsfc(:,:,:,t) = Cmn;
% end
% Conv_mean_nearsfc = permute(Conv_mean_nearsfc,[1 2 4 3]);
% clear Cmn dummy Conv_mean
% 
% 
% 
% % figure; contourf(Conv_mean_nearsfc(:,:,15,200)',[-0.01:0.0001:0.01])
% 
% i_max_nearsfc_meanfield = zeros(Yvcsa,Tvcsa); i_max_nearsfc(:) = NaN;
% k_max_nearsfc_meanfield = i_max_nearsfc;
% for t = 1:Tvcsa   %keep 1:1 when dealing with alltime mean field
%     for j = 1:Yvcsa
%         [imaxgrnd kmaxgrnd] = find(   Conv_mean_nearsfc(:,:,t,j) == max(max(Conv_mean_nearsfc(:,:,t,j)))  ) ;
%         if( isempty(imaxgrnd)==0 )
%             i_max_nearsfc_meanfield(j,t) = imaxgrnd(1);
%             k_max_nearsfc_meanfield(j,t) = kmaxgrnd(1);
%         end
%     end
% end
% clear Conv_mean_nearsfc
% 
% %   figure; contourf(Conv_mean_nearsfc(:,:,15,200)',[-0.01:0.0001:0.01],'LineColor','none'); hold on; plot(i_max_nearsfc_meanfield(200,15),k_max_nearsfc_meanfield(200,15),'ok')
% %   figure; plot(i_max_nearsfc_meanfield(:,15),[1:1600],'ok');
% 
% 
% 
% 
% 
% 
% 
% %FOR AGL MEAN DECTION: find the index of maxima in conv_mean AGL field that is very near the ground and is the largest in targeted subdomain:
% 
% disp('    ')
% disp(' locating the time mean-ed conv max near the terrain  ')
% disp('    ')
% tic
% 
% %temp conv field that is just near the surface:
% Conv_meanagl_nearsfc = permute(Convagl_mean,[1 2 3 4]);
% clear Convagl_mean
% Conv_meanagl_nearsfc(:,6:end,:,:)  = NaN;
% 
% i_max_nearsfcagl_meanfield = zeros(Yvcsa,Tvcsa); i_max_nearsfc(:) = NaN;
% k_max_nearsfcagl_meanfield = i_max_nearsfc;
% for t = 1:Tvcsa   %keep 1:1 when dealing with alltime mean field
%     for j = 1:Yvcsa
%         [imaxgrnd kmaxgrnd] = find(   Conv_meanagl_nearsfc(:,:,t,j) == max(max(Conv_meanagl_nearsfc(:,:,t,j)))  ) ;
%         if( isempty(imaxgrnd)==0 )
%             i_max_nearsfcagl_meanfield(j,t) = imaxgrnd(1);
%             k_max_nearsfcagl_meanfield(j,t) = kmaxgrnd(1);
%         end
%     end
% end
% clear Conv_meanagl_nearsfc
% 
% %convert k agl ito k asl via interpolation:
% i_max_nearsfc_agl2asl_meanfield =i_max_nearsfcagl_meanfield;
% k_max_nearsfc_agl2asl_meanfield = k_max_nearsfcagl_meanfield; k_max_nearsfc_agl2asl_meanfield(:) = NaN;
% for j = 1:Yvcsa
%     for t = 1:Tvcsa
%         %k_max_nearsfc_agl2asl_meanfield(j,t)  =  interp1(Ht,k_max_nearsfcagl_meanfield(j,t), AGLinASL(i_max_nearsfc_agl2asl_meanfield(j,t),j,:));
%         aslhp = AGLinASL(i_max_nearsfc_agl2asl_meanfield(j,t),j,:);
%         aglh = Ht(k_max_nearsfcagl_meanfield(j,t));  %m
%         k_max_nearsfc_agl2asl_meanfield(j,t)  =  find(  abs(aslhp - aglh) == min(min(abs(aslhp - aglh)))        );
%     end
% end
% % do this if you want to detect mean swath position in AGL frame, but
% % detect in instant conv/w fields in ASL frame
% i_max_nearsfc_meanfield = i_max_nearsfc_agl2asl_meanfield;
% k_max_nearsfc_meanfield = k_max_nearsfc_agl2asl_meanfield;
% 
% 
% % do this if you want to multi-pass movmean smooth the conv_mean locations along j
% %note, this may put some points underground?
% i_max_nearsfc_meanfield_unsmooth = i_max_nearsfc_meanfield;
% i_max_nearsfc_meanfield_sm = movmean(i_max_nearsfc_meanfield,11,1);
% k_max_nearsfc_meanfield_sm = movmean(k_max_nearsfc_meanfield,11,1);
% for n = 1:400
%     i_max_nearsfc_meanfield_sm = movmean(i_max_nearsfc_meanfield_sm,11,1);
%     k_max_nearsfc_meanfield_sm = movmean(k_max_nearsfc_meanfield_sm,11,1);
% end
% i_max_nearsfc_meanfield = floor(i_max_nearsfc_meanfield_sm);
% k_max_nearsfc_meanfield = floor(k_max_nearsfc_meanfield_sm);
% %reset k_max to the 3rd pt above teh ground if the smoothed one is below ground
% for j = 1:Yvcsa
%     for t = 1:Tvcsa
%         % t = 1; j = 1332;
%         aslhp = AGLinASL(i_max_nearsfc_meanfield(j,t),j,k_max_nearsfc_meanfield(j,t));  %asl height of point
%         localaglprof = A_AGL_vcs(i_max_nearsfc_meanfield(j,t),:,j);
%         aboveground = find(localaglprof>=0);
%         if(  isnan(aslhp)  )
%             k_max_nearsfc_meanfield(j,t) = aboveground(3);
%         end
%     end
% end
% 
% clear dummy Conv_mean_nearsfc AGLinASL
% 
% % figure;
% % caglm = permute(Convagl_mean,[1 4 2 3]);
% % contourf( max(caglm(:,:,:,6),[],3,'omitnan')', 20 ,'LineColor','none')
% % hold on
% % plot(i_max_nearsfc_meanfield_sm,[1:1600],'ok')
% % plot(i_max_nearsfc_meanfield_unsmooth,[1:1600],'or')
% 
% 
% 
% %   caglm = permute(Convagl_alltmean,[1 4 2 3]);
% %   figure; contourf( caglm(:,:,3,5)',[-0.01:0.0001:0.01],'LineColor','none'); hold on;  plot(i_max_nearsfcagl_meanfield(:,5),[1:1600],'ok');
% 
% 
% 
% % %     %diagnostic compare conv_agl vs conv_asl:
% % %
% % %     %agl
% % caglm = permute(Convagl_alltmean,[1 4 2 3]);   caglm(caglm <= -1)=NaN;
% %     asdf = Conv;
% %     asdf = movmean(asdf,11,1,'omitnan');
% %     asdf = movmean(asdf,11,2,'omitnan');
% % %     %figure; contourf( caslm(:,:,15,5)',[-0.01:0.0001:0.01],'LineColor','none'); hold on;  plot(i_max_nearsfc_meanfield(:,5),[1:1600],'ok');
% % % for t = 1:sa-1
% % %     fig1 = figure('position',[220,387,464,549]);
% % %     %contourf(caslm(:,:,10,t)',[-0.003:0.0003:0.003],'LineColor','none'); caxis([-0.003 0.003]);
% % %     contourf(asdf(:,:,5,t)',[-0.003:0.0003:0.003],'LineColor','none'); caxis([-0.003 0.003]);
% % %     hold on
% % %     contour(Terr'/10000,[1500:250:2500]/10000,'k')
% % %     plot(i_max_nearsfcagl_meanfield(:,5),[1:1600],'ok');
% % %     title('InstConvagl at t',num2str(t))
% % % end
% % %
% %
% %     %asl
% %     caslm = permute(Conv_alltmean,[1 4 2 3]);   caslm(caslm <= -1)=NaN;
% % %     asdf = permute(A_Conv_vcs,[1 4 2 3]);
% % %     asdf = movmean(asdf,11,1,'omitnan');
% % %     asdf = movmean(asdf,11,2,'omitnan');
% %     %figure; contourf( caslm(:,:,15,5)',[-0.01:0.0001:0.01],'LineColor','none'); hold on;  plot(i_max_nearsfc_meanfield(:,5),[1:1600],'ok');
% % for t = 1:sa-1
% %     fig1 = figure('position',[220,387,464,549]);
% %     contourf(caglm(:,:,5,t)',[-0.003:0.0003:0.003],'LineColor','none'); caxis([-0.003 0.003]);
% %     %contourf(asdf(:,:,5,t)',[-0.003:0.0003:0.003],'LineColor','none'); caxis([-0.003 0.003]);
% %     hold on
% %     contour(Terr'/10000,[1500:250:2500]/10000,'k')
% %     plot(i_max_nearsfc_meanfield(:,t),[1:1600],'ok');
% %     title('InstConvasl at t',num2str(t))
% % end
% 
% 
% 
% 
% 
% 
% %%% reset W & conv to raw for inst swath ID
% %A_Conv_vcs = A_Conv_vcs_unsmoothed;
% A_W_vcs = A_W_vcs_unsmoothed;
% %Conv = Conv_unsmoothed;
% 
% clear A_W_vcs_unsmoothed
% 
% %%%%%%ASL FRAME Now start looking in the raw field within a window near the
% %%%%%% conv_mean sfc max, and w-max a little above it
% 
% disp('    ')
% disp(' locating raw near sfc conv/w max in prescribed box  ')
% disp('    ')
% tic
% 
% %seed i,k max raw con fields
% i_max_nearsfc_rawfield = i_max_nearsfc_meanfield;  i_max_nearsfc_rawfield(:) = nan;
% k_max_nearsfc_rawfield = k_max_nearsfc_meanfield;  k_max_nearsfc_rawfield(:) = nan;
% 
% %prescribed search box size from the nearsfc conv_mean max to look within the instantaneous sfc conv field to find its max (to subsequently flood fill upon)
% %hconwindow = 15;
% %vconwindow = 15;
% hrconwindow = 40;  % these are i,k space, not physical space. i = 15 was for d3, so try just mult by 5 to get to d4. (15*5)
% hlconwindow = 40;
% vdconwindow = 15;  %15 was for d3 when i used newdz = 25m. I guess divide by for for d4's newdz=100m
% vuconwindow = 10;  %later changed this to a large number to compensate for the iterative smoothing of i_max,k_max(j) above. should be okay to
% %make this vertical window largeish as long as abovegrn below is limited to just a few 100 m agl.
% 
% 
% 
% %define near-sfc raw conv field
% dummy = permute(A_Conv_vcs,[1 2 4 3]);
% abovegrnd =  find( A_AGL_vcs  > 500 | A_AGL_vcs  < 25);  %physical space thresholds to use to kill
% Conv_raw_nearsfc = dummy; Conv_raw_nearsfc(:) = NaN;
% for t = 1:Tvcsa
%     Cmn = dummy(:,:,:,t);
%     Cmn(abovegrnd) = nan;
%     Conv_raw_nearsfc(:,:,:,t) = Cmn;
% end
% Conv_raw_nearsfc = permute(Conv_raw_nearsfc,[1 2 4 3]);
% Conv_raw_nearsfc(find(Conv_raw_nearsfc < -1)) = NaN;
% clear abovegrnd %A_AGL_vcs
% 
% 
% crn = Conv_raw_nearsfc;
% %   figure; contourf(crn(:,:,15,200)',[-0.01:0.0001:0.01],'LineColor','none'); hold on; ; hold on; plot(i_max_nearsfc_meanfield(200,15),k_max_nearsfc_meanfield(200,15),'ok')
% %disp('bbb')
% 
% clear dummy Conv_raw_nearsfc
% 
% %seed i,k max raw w fields
% i_max_w_rawfield = i_max_nearsfc_rawfield ;
% k_max_w_rawfield = k_max_nearsfc_rawfield ;
% 
% W_dummy = A_W_vcs;
% %Wca_dummy = A_Wca_vcs;
% 
% %prescribed w search box
% %hwwindow = 5;  %i,k space. 5 & 15 were for d3 (dx=500m,dz=25m)
% %vwwindow = 15;
% hwwindow = 15;
% vwwindow = 5;
% 
% 
% for t = 1:Tvcsa   %keep 1:1 when dealing with alltime mean field
%     for j = 1:Yvcsa
% 
%         % j = 1334; t = 1;
% 
%         %CONVERGENCE:
%         %nan-out stuff outside of wanted box:
%         crn(1:i_max_nearsfc_meanfield(j,t)-hlconwindow,:,t,j) = nan;
%         crn(i_max_nearsfc_meanfield(j,t)+hrconwindow:Xvcsa,:,t,j) = nan;
% 
%         crn(:,1:k_max_nearsfc_meanfield(j,t)-vdconwindow,t,j) = nan;
%         crn(:,k_max_nearsfc_meanfield(j,t)+vuconwindow:Zvcsa,t,j) = nan;
% 
%         [imaxgrnd kmaxgrnd] = find( crn(:,:,t,j) == max( max( crn(:,:,t,j),[],'omitnan' ),[],'omitnan')  );
% 
%         %{
%             if(j == 868)
% 		    disp(['j868 t=',num2str(t)])
%                     disp(['Xvcsa=',num2str(Xvcsa),' Zvcsa=',num2str(Zvcsa)])
%                     disp('i_max_nearsfc_meanfield(j,t)' )
% 		    i_max_nearsfc_meanfield(j,t)
% 		    disp('i_max_nearsfc_meanfield(j,t)-hlconwindow' )
% 		    i_max_nearsfc_meanfield(j,t)-hlconwindow
%                     disp('i_max_nearsfc_meanfield(j,t)+hrconwindow')
% 		    i_max_nearsfc_meanfield(j,t)+hrconwindow
%                     disp(' k_max_nearsfc_meanfield(j,t)-vdconwindow ')
%                     k_max_nearsfc_meanfield(j,t)-vdconwindow 
% 		    disp('k_max_nearsfc_meanfield(j,t)+vuconwindow')
% 		    k_max_nearsfc_meanfield(j,t)+vuconwindow
% 		    disp('  ')
% 		    
% 		    disp('imaxgrnd')
% 		    imaxgrnd
% 		    disp('  ')
% 		    
% 		    disp('kmaxgrnd')
% 		    kmaxgrnd
% 		    disp('  ')		    
% 
%                     disp(' size crn(:,:,t,j)  ')
% 		    size(crn(:,:,t,j))
% 		    disp('   ')
% 
% 
%                     disp('max max  crn(:,:,t,j) ')
%                     max( max( crn(:,:,t,j),[],'omitnan' ),[],'omitnan') )
%                     disp('  ')
% 
% 
%                     disp(' crn(:,:,t,j) ')
%                     crn(:,:,t,j)
%                     disp('  ')  
% 		    
% 	    end
%         %}
% 
% 
% 
%         if(isempty(imaxgrnd)==0)
%             i_max_nearsfc_rawfield(j,t) = imaxgrnd(1);
%             k_max_nearsfc_rawfield(j,t) = kmaxgrnd(1);
%         end
% 
%         %UPDRAFT:
%         %             %nan-out stuff outside of wanted box - mostly located above the near-sfc conv_mean max:
%         %             W_dummy(1:i_max_nearsfc_meanfield(j,t)-hwwindow,:,t,j) = nan;
%         %             W_dummy(i_max_nearsfc_meanfield(j,t)+hwwindow:Xvcsa,:,t,j) = nan;
%         %             W_dummy(:,1:k_max_nearsfc_meanfield(j,t)-1,t,j) = nan;
%         %             W_dummy(:,k_max_nearsfc_meanfield(j,t)+vwwindow:Zvcsa,t,j) = nan;
% 
%         %       ....OR....
% 
% 
%         %nan-out stuff outside of wanted box - mostly located above the near-sfc conv_mean max:
% 
%         %disp(['i=',num2str(i),' j=',num2str(j),' t=',num2str(t)])
%         %disp(['i_max_nearsfc_rawfield(j,t) = ', num2str(i_max_nearsfc_rawfield(j,t)) ])
%         %disp(['hwwindow',num2str(hwwindow)])
% 
%         W_dummy(1:i_max_nearsfc_rawfield(j,t)-hwwindow,:,t,j) = nan;
%         W_dummy(i_max_nearsfc_rawfield(j,t)+hwwindow:Xvcsa,:,t,j) = nan;
%         W_dummy(:,1:k_max_nearsfc_rawfield(j,t)-1,t,j) = nan;
%         W_dummy(:,k_max_nearsfc_meanfield(j,t)+vwwindow:Zvcsa,t,j) = nan;
%         % % % figure; contourf(cnv(:,:,t,j)',30)
% 
%         [imaxgrnd kmaxgrnd] = find( W_dummy(:,:,t,j) == max(max( W_dummy(:,:,t,j),[],'omitnan' ))) ;
% 
%         if(isempty(imaxgrnd)==0)
%             i_max_w_rawfield(j,t) = imaxgrnd(1);
%             k_max_w_rawfield(j,t) = kmaxgrnd(1);
%         end
% 
%     end  %j
% end  %t
% clear crn Cmn W_dummy Wca_dummy
% 
% toc
% disp('   ')
% 
% 
% 
% 
% disp('    ')
% disp('  flood-filling on instantaneous conv & w fields to ID swath objects  ')
% disp('    ')
% tic
% 
% % now try flood-filling on the instantaneous full CONV field near sfc max to mask the Conv-swath - where conv > prescribed amount  (Conv_object):
% Conv_bin = A_Conv_vcs;
% Conv_bin(A_Conv_vcs > 0.00075) = 1;
% Conv_bin(A_Conv_vcs < 0.00075) = 0;
% Conv_bin(isnan(A_Conv_vcs)) = 0;
% clear A_Conv_vcs
% Conv_bin = logical(Conv_bin);
% Conv_object = Conv_bin; Conv_object(:) = 0;
% 
% for t = 1:Tvcsa
%     for j = 1:Yvcsa
%         if( isnan(k_max_nearsfc_rawfield(j,t))==0 )
%             Conv_object(:,:,t,j) = bwselect(Conv_bin(:,:,t,j),k_max_nearsfc_rawfield(j,t),i_max_nearsfc_rawfield(j,t),4) ; % result = mask of 0,1
%         end
%     end
% end
% 
% clear Conv_bin
% 
% % now try flood-filling on the instantaneous full w field near sfc max to mask the W-swath - where w > prescribed amount  (w_object).
% % as of v6, the vertical component of flood filling only goes up to local MULFC
% W_bin = A_W_vcs;
% W_bin(A_W_vcs > 0.5) = 1;
% W_bin(A_W_vcs < 0.5) = 0;
% W_bin(isnan(A_W_vcs)) = 0;
% W_bin = logical(W_bin);
% W_object = W_bin; W_object(:) = 0;
% 
% for t = 1:Tvcsa
%     for j = 1:Yvcsa
%         if( isnan(k_max_w_rawfield(j,t)) == 0 )
% 
%             %trim out stuff above MU LFC
%             for i = 1:Xvcsa
%                 ceiling  = find(      abs(A_Ht_asl_vcs(i,:,j) - MULFC_all(i,j,t)) == min( abs(A_Ht_asl_vcs(i,:,j) - MULFC_all(i,j,t) ))         )  ;
%                 W_bin(i,ceiling:end,t,j) = 0;
%             end
% 
%             W_object(:,:,t,j) = bwselect(W_bin(:,:,t,j),k_max_w_rawfield(j,t),i_max_w_rawfield(j,t),4) ; % result = mask of 0,1
%             %   figure;   contourf(W_object(:,:,t,j)',4);    hold on;    plot( [1:Xvcsa],ceilingk,'ow')
%         end
%     end
% end
% 
% clear W_bin ceiling
% 
% 
% 
% % % diagnostic plots to test 100m conv/w detection:
% wobj = permute(W_object,[1 4 2 3]);
% %casl = permute(A_Conv_vcs,[1 4 2 3]);
% %cobj = permute(Conv_object,[1 4 2 3]);
% %wasl = permute(A_W_vcs,[1 4 2 3]);
% 
% %CMallt = permute(Conv_mean,[1 4 2 3]);
% % CM11ijallt = permute(Conv_mean_11imean_11jmean,[1 4 2 3]);
% % CM5t = permute(Conv_5tmean,[1 4 2 3]);
% % CM9t = permute(Conv_9tmean,[1 4 2 3]);
% % CM13t = permute(Conv_13tmean,[1 4 2 3]);
% 
% %CMallt(CMallt == -999) = NaN;
% % CM11ijallt(CM11ijallt < -1) = NaN;
% % CM5t(CM5t == -999) = NaN;
% % CM9t(CM9t == -999) = NaN;
% % CM13t(CM13t == -999) = NaN;
% 
% 
% 
% 
% 
% %{
% % % swath diagnostic plots kagl = 5;
% 
% %casl = permute(A_Conv_vcs,[1 4 2 3]); for t = 1:sa-1; fig1 =
% figure('position',[220,387,464,549]);
% contourf(max(Conv(:,:,1:kagl+4,t),[],3)',20,'LineColor','none');
% caxis([-0.01 0.01]);
% %contourf(max(casl(:,:,20,t),[],3)',20,'LineColor','none'); caxis([-0.01
% 0.01]); hold on contour(Terr'/10000,[1500:250:2500]/10000,'k') wom =
% max(wobj(:,:,:,t),[],3); contour(wom'/10000,'r','LineWidth',0.5) %   com =
% max(cobj(:,:,:,t),[],3); %   contour(com'/10000,'g')
% 
% plot(i_max_nearsfc_meanfield(:,t),[1:1600],'k*')
% plot(i_max_nearsfc_meanfield(:,t)+hrconwindow,[1:1600],'k<','MarkerSize',0.5)
% plot(i_max_nearsfc_meanfield(:,t)-hlconwindow,[1:1600],'k>','MarkerSize',0.5)
% plot(i_max_nearsfc_rawfield(:,t),[1:1600],'w*','MarkerSize',0.75)
% 
% title( [ 'agl2aslcmean 11ijsmbothC t13mean1hr cmeansm400j rawWinst
% hcwinl40r40 vcwinu5d5 hwwin5 vwwin5 cthr0.00075 wthr0.5 Convagl wobj at
% ',YYMMDD_all(1,:),' ',HHMMSS_all(t,:) ] ) pic =
% strcat(imout,YYMMDD_all(1,:),'_',HHMMSS_all(t,:),'gefs18base_B_agl2aslcmean_11ijsmbothC_t13mean1hr_cmeansm400j_rawWinst_hcwinl40r40_vcwinu5d5_hwwin5_vwwin5_cthr0.00075_wthr0.5_Convagl_wobj',num2str(kagl),'t',num2str(t),'.png');
% saveas(fig1,pic) close(fig1)
% 
% end
% 
% %%%% plot swath kagl = 5;
% 
% %casl = permute(A_Conv_vcs,[1 4 2 3]); for t = 1:sa-1; fig1 =
% figure('position',[220,387,464,549]);
% contourf(max(Conv(:,:,1:kagl+4,t),[],3)',20,'LineColor','none');
% caxis([-0.01 0.01]);
% %contourf(max(casl(:,:,20,t),[],3)',20,'LineColor','none'); caxis([-0.01
% 0.01]); hold on contour(Terr'/10000,[1500:250:2500]/10000,'k') %    wom =
% max(wobj(:,:,:,t),[],3); %    contour(wom'/10000,'r','LineWidth',0.5) %
% com = max(cobj(:,:,:,t),[],3); %   contour(com'/10000,'g')
% 
% plot(i_max_nearsfc_meanfield(:,t),[1:1600],'k*')
% plot(i_max_nearsfc_meanfield(:,t)+hrconwindow,[1:1600],'k<','MarkerSize',0.5)
% plot(i_max_nearsfc_meanfield(:,t)-hlconwindow,[1:1600],'k>','MarkerSize',0.5)
% plot(i_max_nearsfc_rawfield(:,t),[1:1600],'wo','MarkerSize',0.75)
% 
% title( [ 'agl2aslcmean 11ijsmbothC t13mean1hr cmeansm400j rawWinst
% hcwinl40r40 vcwinu5d5 hwwin5 vwwin5 cthr0.00075 wthr0.5 Convagl at
% ',YYMMDD_all(1,:),' ',HHMMSS_all(t,:) ] ) pic =
% strcat(imout,YYMMDD_all(1,:),'_',HHMMSS_all(t,:),'gefs18base_agl2aslcmean_11ijsmbothC_t13mean1hr_cmeansm400j_rawWinst_hcwinl40r40_vcwinu5d5_hwwin5_vwwin5_cthr0.00075_wthr0.5_Convagl',num2str(kagl),'t',num2str(t),'.png');
% saveas(fig1,pic) close(fig1)
% 
% end
% %}
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %{
% 
% %%% VCS plots
% t = 3;
% j = 857;
% fig1 = figure; %('position',[220,387,464,549]);
% contourf(A_Conv_vcs(:,:,t,j)',20,'LineColor','none'); caxis([-0.01 0.01]);
% hold on
% contour(W_object(:,:,t,j)','r')
% plot(i_max_nearsfc_meanfield(j,t),k_max_nearsfc_meanfield(j,t),'ks','MarkerSize',10)
% plot(i_max_nearsfc_meanfield(j,t)-hlconwindow,k_max_nearsfc_meanfield(j,t),'ks','MarkerSize',5)
% plot(i_max_nearsfc_meanfield(j,t)+hrconwindow,k_max_nearsfc_meanfield(j,t),'ks','MarkerSize',5)
% plot(i_max_nearsfc_rawfield(j,t),k_max_nearsfc_rawfield(j,t),'g*','MarkerSize',10)
% plot(i_max_w_rawfield(j,t),k_max_w_rawfield(j,t),'r^','MarkerSize',10)
% 
% 
% 
% kasl = 20;
% fig2 = figure('position',[220,387,464,549]);
% contourf(max(CMallt(:,:,:,t),[],3)',20,'LineColor','none'); caxis([-0.008 0.008]);
% hold on
% contour(Terr'/1000000,[1500:250:2500]/1000000,'k')
% % wom = max(wobj(:,:,:,t),[],3);
% % contour(wom'/10000,'r')
%     plot(i_max_nearsfc_meanfield,[1:1600],'k*')
% title('Convasl_alltmean')
% % pic = strcat(wrfdir,'Convasl_alltmean',num2str(kasl),'t',num2str(t),'.png');
% % saveas(fig2,pic)
% 
% 
% fig3 = figure('position',[220,387,464,549]);
% contourf(CM11ijallt(:,:,kasl,t)',20,'LineColor','none'); caxis([-0.005 0.005]);
% hold on
% contour(Terr'/1000000,[1500:250:2500]/1000000,'k')
% wom = max(wobj(:,:,:,t),[],3);
% contour(wom'/10000,'r')
% title('Convasl_ij11smooth_alltmean')
% % pic = strcat(wrfdir,'Convasl_ij11smooth_alltmean',num2str(kasl),'t',num2str(t),'.png');
% % saveas(fig3,pic)
% 
% 
% fig4 = figure('position',[220,387,464,549]);
% contourf(CM9t(:,:,kasl,t)',20,'LineColor','none'); caxis([-0.005 0.005]);
% hold on
% contour(Terr'/1000000,[1500:250:2500]/1000000,'k')
% wom = max(wobj(:,:,:,t),[],3);
% contour(wom'/10000,'r')
% title('Convasl_9tmean')
% % pic = strcat(wrfdir,'Convasl_9tmean',num2str(kasl),'t',num2str(t),'.png');
% % saveas(fig4,pic)
% 
% 
% 
% %}
% 
% 
% 
% %{
%     disp(' saving objects  ')
% 
%     %Wob1 = W_object(:,:,1,:);
%     %Wob6 = W_object(:,:,6,:);    
%     Wob11 = W_object(:,:,11,:);
% 
%     %Wagl1 = W_all(:,:,:,1);
%     %Wagl6 = W_all(:,:,:,6);
%     Wagl11 = W_all(:,:,:,11);
% 
%     %Wasl1 = A_W_vcs(:,:,1,:);
%     %Wasl6 = A_W_vcs(:,:,6,:);
%     Wasl11 = A_W_vcs(:,:,11,:); 
% 
%     %Cagl1 = Conv(:,:,:,1);
%     %Cagl6 = Conv(:,:,:,6);
%     Cagl11 = Conv(:,:,:,11);
% 
%     %Casl1 = A_Conv_vcs(:,:,1,:);
%     %Casl6 = A_Conv_vcs(:,:,6,:);
%     Casl11 = A_Conv_vcs(:,:,11,:);
% 
%     %Cob1  = Conv_object(:,:,1,:);
%     %Cob6  = Conv_object(:,:,6,:);
%     Cob11 = Conv_object(:,:,11,:);
% 
%     Cma = Conv_alltmean(:,:,11,:);
%     Cmt5 = Conv_5tmean(:,:,11,:);
%     Cmt9 = Conv_9tmean(:,:,11,:);
%     Cmt13 = Conv_13tmean(:,:,11,:);
% 
%     save( strcat(wrfdir,'testobj.mat'), 'Lat', 'Lon', 'Wagl11', 'Wasl11', 'Cagl11', 'Casl11',...
%      	    'Terr', 'AGLinASL', 'Wob11','Cob11', 'Cma', 'Cmt5', 'Cmt9', 'Cmt13',  '-v7.3' ) 
% %}
% 
% 
% 
% 
% 
% 
% % %     % now try flood-filling on the instantaneous clear-air w field near sfc max to mask the W-swath - where w > prescribed amount (wca_object):
% % %     Wca_bin = A_Wca_vcs;
% % %     Wca_bin(A_Wca_vcs > 0.5) = 1;
% % %     Wca_bin(A_Wca_vcs < 0.5) = 0;
% % %     Wca_bin(isnan(A_Wca_vcs)) = 0;
% % %     Wca_bin = logical(Wca_bin);
% % %     Wca_object = Wca_bin; Wca_object(:) = 0;
% % %
% % %     for t = 1:Tvcsa
% % %         for j = 1:Yvcsa
% % %             if( isnan(k_max_w_rawfield(j,t)) == 0 )
% % %                 Wca_object(:,:,t,j) = bwselect(Wca_bin(:,:,t,j),k_max_w_rawfield(j,t),i_max_w_rawfield(j,t),4) ; % result = mask of 0,1
% % %             end
% % %         end
% % %     end
% % %
% % %     toc
% % %     disp('    ')
% 
% 
% 
% 
% 
% 
% 
% 
% %     % % diagnostic plots of W & Conv objects
% %     %   figure; contourf(regxasl(:,1,1),newz_asl/1000,Wca_object(:,:,15,300)',20)
% %
% %     %  debug clear of unneeded vars for memory:
% %     %clear A_QVca_vcs ATHEca_vcs A_THca_vcs A_DBZ_vcs A_Qcld_vcs abovegrnd %MUCAPE_all MUCIN_all
% %     %clear Qcld_all QV_3Dswath RH_ASL2 RHca_ASL2 THETA_ASL2 THETAE_ASL2  A_RHca_vcs A_THEca_VCS %THEca_ASL2
% 
% % %     % x-y plane diagnositc plots:
% % %     convinsta_3dtmp = permute(A_Conv_vcs,[1 4 2 3]);
% % %     wobj_2dtmp = permute(W_object,[1 4 2 3]);
% % %     wcaobj_2dtmp = permute(Wca_object,[1 4 2 3]);
% % %     cobj_2dtmp = permute(Conv_object,[1 4 2 3]);
% % %
% % %     k = 100;
% % %     t = 11;
% % %
% % %
% % %     figure;
% % %     contourf(convinsta_3dtmp(:,:,k,t)',[-0.01:0.001:0.01],'LineColor','none'); caxis([-0.003 0.003])
% % %     hold on
% % %     mxcobj = cobj_2dtmp(:,:,:,t)/100;
% % %     mxwobj = wobj_2dtmp(:,:,:,t)/100;
% % %     mxcobj = max(mxcobj,[],3);
% % %     mxwobj = max(mxwobj,[],3);
% % %     contour( mxcobj', [0.99 0.99]/100,  'LineColor', 'm', 'LineWidth', 1.5);
% % %     contour( mxwobj', [0.99 0.99]/100,  'LineColor', 'r', 'LineWidth', 1.0);
% % %     contour(Terr'/100000,[1500:250:2500]/100000,'k')
% % %     title(['Insta conv at k',num2str(k),' at t ',num2str(t)])
% % %
% % %
% % %     figure;
% % %     convvcomp_2dtmp = convinsta_3dtmp(:,:,:,t);
% % %     convvcomp_2dtmp = max(convvcomp_2dtmp,[],3);
% % %     contourf(convvcomp_2dtmp',[-0.01:0.001:0.01],'LineColor','none'); caxis([-0.003 0.003])
% % %     hold on
% % %     mxcobj = cobj_2dtmp(:,:,:,t)/100;
% % %     mxwobj = wobj_2dtmp(:,:,:,t)/100;
% % %     mxcobj = max(mxcobj,[],3);
% % %     mxwobj = max(mxwobj,[],3);
% % %     contour( mxcobj', [0.99 0.99]/100,  'LineColor', 'm', 'LineWidth', 1.5);
% % %     contour( mxwobj', [0.99 0.99]/100,  'LineColor', 'r', 'LineWidth', 1.0);
% % %     contour(Terr'/100000,[1500:250:2500]/100000,'k')
% % %     title(['Insta conv at k',num2str(k),' at t ',num2str(t)])
% % %
% % %
% % %     for j = 120:10:240
% % %
% % %         height = A_Ht_asl_vcs(:,:,j) ;
% % %         Xx = A_X_vcs(:,:,j) ;
% % %         lfc = MULFC_all(:,j,t)/1000;
% % %
% % %
% % %         figure
% % %         contourf(Xx,height/1000, A_Conv_vcs(:,:,t,j),[-0.01:0.001:0.01],'LineColor','none'); caxis([-0.003 0.003])
% % %         hold on
% % %         contour(Xx,height/1000, permute(wobj_2dtmp(:,j,:,t)/100,[1 3 2]), [0.99 0.99]/100,  'LineColor', 'r', 'LineWidth', 1);
% % %         contour(Xx,height/1000, permute(wcaobj_2dtmp(:,j,:,t)/100,[1 3 2]), [0.99 0.99]/100,  'LineColor', 'c', 'LineWidth', 1, 'LineStyle','--');
% % %         contour(Xx,height/1000, permute(cobj_2dtmp(:,j,:,t)/100,[1 3 2]), [0.99 0.99]/100,  'LineColor', 'm', 'LineWidth', 1.5);
% % %         plot(Xx(:,1),lfc,':w','LineWidth',1.5)
% % %         %contour(Terr'/100000,[1500:250:2500]/100000,'k')
% % %         title(['Insta conv at j',num2str(j),' at t ',num2str(t)])
% % %         %axis([80 140 20 160])
% % %         axis([-30 10 0 5])
% % %
% % %     end
% 
% 
% 
% %%%%%%%%%%%%%
% %%%% Identify the 'central' longitude of the swath at each lat (full field, not clear air):
% %%%%%%%%%%%%%
% 
% size(W_object)   ;
% 
% W_object_vertcomp = max( W_object, [], 2 );   size(W_object_vertcomp);
% W_object_vertcomp  = permute(W_object_vertcomp, [1 4 3 2]);   size(W_object_vertcomp);
% 
% swath_centlongitude = zeros(Yvcsa,Tvcsa) ;   swath_centlongitude(:) = NaN;
% swath_centlatitude = zeros(Yvcsa,Tvcsa)  ;   swath_centlatitude(:) = NaN;
% swath_j = zeros(Yvcsa,Tvcsa)  ;   swath_j(:) = NaN;
% 
% %call the central lat/lon the mean within this j of the swath's span.
% for j = 1:Yvcsa
%     for t = 1:Tvcsa
%         swathi = find(W_object_vertcomp(:,j,t)) ;
%         swath_centlongitude(j,t) =  mean( Lon(swathi,j) ,'omitnan')  ;
%         swath_centlatitude(j,t)  =  mean( Lat(swathi,j) ,'omitnan')  ;
%         swath_j(j,t)  =  j  ;
%     end
% end
% clear W_object_vertcomp
% 
% 
% % interpolate values to fill in across the full j=[1:Yvcsa] span:
% swath_centlon_interp = zeros(Yvcsa,Tvcsa); swath_centlon_interp(:) = NaN;
% swath_centlat_interp = zeros(Yvcsa,Tvcsa); swath_centlat_interp(:) = NaN;
% for t = 1:Tvcsa
%     tmpj = swath_j(:,t);
%     tmplon = swath_centlongitude(:,t);
%     tmplat = swath_centlatitude(:,t);
%     tmpj(isnan(swath_j(:,t))) = [];  tmplon(isnan(swath_j(:,t))) = [];   tmplat(isnan(swath_j(:,t))) = [];
% 
%     swath_centlon_interp(:,t) = interp1(tmpj,tmplon,[1:Yvcsa]);
%     swath_centlat_interp(:,t) = interp1(tmpj,tmplat,[1:Yvcsa]);
% end
% clear swath_centlongitude swath_centlatitude swath_j
% 
% 
% % smooth the swath central lat/lon field:
% dummy = swath_centlon_interp; dummy(:) = NaN;
% Wswath_centlat_sm1 = dummy; Wswath_centlon_sm1 = dummy;
% Wswath_centlat_sm2 = dummy; Wswath_centlon_sm2 = dummy;
% WIN = 9;  %movmean window length
% for t = 1:1:Tvcsa
%     Wswath_centlat_sm1(:,t)    = movmean(swath_centlat_interp(:,t),WIN,'omitnan');
%     Wswath_centlat_sm2(:,t)    = movmean(Wswath_centlat_sm1(:,t),WIN,'omitnan');
% 
%     Wswath_centlon_sm1(:,t)    = movmean(swath_centlon_interp(:,t),WIN,'omitnan');
%     Wswath_centlon_sm2(:,t)    = movmean(Wswath_centlon_sm1(:,t),WIN,'omitnan');
% end
% clear Wswath_centlon_sm1 Wswath_centlat_sm1
% 
% 
% 
% %%%%%%%%%% define some thermo quantities within the 3d w-swath mask here:
% 
% disp('    ')
% disp('  define some thermo quantities within the 3d w-swath mask  ')
% disp('    ')
% tic
% 
% % %     %clear-air version:
% % %     Wca_3Dmask = single(Wca_object); Wca_3Dmask(find(Wca_3Dmask==0)) = NaN;
% % %     QVca_3Dswath = permute(QVca_ASL2,[1 3 4 2]) .* Wca_3Dmask;   % 4D thermo field nan'ed outside of the 3D swath mask
% % %     THca_3Dswath = permute(THca_ASL2,[1 3 4 2]) .* Wca_3Dmask;
% % %     THEca_3Dswath = permute(THEca_ASL2,[1 3 4 2]) .* Wca_3Dmask;
% 
% %     %full air version:
% W_3Dmask = single(W_object); W_3Dmask(find(W_3Dmask==0)) = NaN;
% QV_3Dswath = permute(QVAP_ASL2,[1 3 4 2]) .* W_3Dmask;   % 4D thermo field nan'ed outside of the 3D swath mask
% TH_3Dswath = permute(THETA_ASL2,[1 3 4 2]) .* W_3Dmask;
% THE_3Dswath = permute(THETAE_ASL2,[1 3 4 2]) .* W_3Dmask;
% 
% %VMca_3Dswath = permute(VMca_ASL2,[1 3 4 2]) .* W_3Dmask;
% VM_3Dswath = permute(VM_ASL2,[1 3 4 2]) .* W_3Dmask;
% 
% save('/global/cfs/projectdirs/m1657/jmarquis/LASSO/matlab/checkpt1.mat','W_3Dmask','QV_3Dswath')
% 
% 
% toc
% disp('   ')
% 
% clear W_3Dmask VM_ASL2 QVAP_ASL2 THETA_ASL2 THETAE_ASL2
% 
% 
% % now start cataloging depth, width, strength, and thermo within w swath:
% 
% % convergence (temporary punt):
% % Convswath_width = Wswath_width;  Convswath_width(:) = nan;
% % Convswath_depth = Wswath_depth;  Convswath_depth(:) = nan;
% 
% 
% %seed width/depth fields using swath of full w field:
% Wswath_width = zeros(Yvcsa,Zvcsa,Tvcsa); Wswath_width(:) = nan;
% Wswath_depth = zeros(Yvcsa,Tvcsa);  Wswath_depth(:) = nan;
% Wswath_maxmag = zeros(Yvcsa,Tvcsa);  Wswath_maxmag(:) = nan;
% Wswath_toplfcagl = zeros(Yvcsa,Tvcsa);  Wswath_toplfcagl(:) = nan;
% Wswath_toplfcasl = zeros(Yvcsa,Tvcsa);  Wswath_toplfcasl(:) = nan;
% 
% %vertical mean of 'full-air' thermo in each swath along-lat slice  % PUTBACKINLATER
% Wswath_thetamean = zeros(Yvcsa,Tvcsa); Wswath_thetamean(:) = nan;
% Wswath_thetaemean = zeros(Yvcsa,Tvcsa);  Wswath_thetaemean(:) = nan;
% Wswath_vapormean = zeros(Yvcsa,Tvcsa);  Wswath_vapormean(:) = nan;
% 
% 
% % %     %seed width/depth fields using swath of clear air  w field:
% % %     Wcaswath_width = zeros(Yvcsa,Zvcsa,Tvcsa); Wcaswath_width(:) = nan;
% % %     Wcaswath_depth = zeros(Yvcsa,Tvcsa);  Wcaswath_depth(:) = nan;
% % %     Wcaswath_maxmag = zeros(Yvcsa,Tvcsa);  Wcaswath_maxmag(:) = nan;
% % %     %vertical mean of 'clear-air' thermo in each swath along-lat slice
% % %     Wcaswath_thetamean = zeros(Yvcsa,Tvcsa); Wcaswath_thetamean(:) = nan;
% % %     Wcaswath_thetaemean = zeros(Yvcsa,Tvcsa);  Wcaswath_thetaemean(:) = nan;
% % %     Wcaswath_vapormean = zeros(Yvcsa,Tvcsa);  Wcaswath_vapormean(:) = nan;
% %
% %     % top & bot of the w swath along each lat:
% % %     Wcaswath_top_asl = zeros(Yvcsa,Tvcsa);  Wcaswath_top_asl(:) = nan;
% % %     Wcaswath_top_agl = zeros(Yvcsa,Tvcsa);  Wcaswath_top_agl(:) = nan;
% Wswath_top_asl = zeros(Yvcsa,Tvcsa);  Wswath_top_asl(:) = nan;
% Wswath_top_agl = zeros(Yvcsa,Tvcsa);  Wswath_top_agl(:) = nan;
% 
% % %     Wcaswath_bot_asl = zeros(Yvcsa,Tvcsa);  Wcaswath_bot_asl(:) = nan;
% % %     Wcaswath_bot_agl = zeros(Yvcsa,Tvcsa);  Wcaswath_bot_agl(:) = nan;
% Wswath_bot_asl = zeros(Yvcsa,Tvcsa);  Wswath_bot_asl(:) = nan;
% Wswath_bot_agl = zeros(Yvcsa,Tvcsa);  Wswath_bot_agl(:) = nan;
% 
% Wswath_ktop_asl = zeros(Yvcsa,Tvcsa);  Wswath_ktop_asl(:) = nan;
% Wswath_kbot_asl = zeros(Yvcsa,Tvcsa);  Wswath_kbot_asl(:) = nan;
% % %     Wcaswath_ktop_asl = zeros(Yvcsa,Tvcsa);  Wcaswath_ktop_asl(:) = nan;
% % %     Wcaswath_kbot_asl = zeros(Yvcsa,Tvcsa);  Wcaswath_kbot_asl(:) = nan;
% %
% % %     WMcaflux_sum = zeros(Yvcsa,Tvcsa);  WMcaflux_sum(:) = nan;
% WMflux_sum = zeros(Yvcsa,Tvcsa);  WMflux_sum(:) = nan;
% 
% disp('    ')
% disp('  cataloging depth, width, strength, and thermo within w swath  ')
% disp('    ')
% tic
% %populate them:
% for t = 1:Tvcsa
%     for j = 1:Yvcsa
% 
%         %%%%%% full-air W swath stuff:
% 
%         % find width of full air w swath along j and t and k:
%         %             VWcaflux_perz = [1:Zvcsa]' ;  VWcaflux_perz(:) = NaN;
%         VWflux_perz = [1:Zvcsa]' ;  VWflux_perz(:) = NaN;
% 
%         for k = 1:Zvcsa
%             % k = 80;
%             wobj = W_object(:, k, t, j) ; % e-w slice of w swath at this height and time
%             wswath = find( wobj == 1 );
%             if( isempty(wswath) == 0  &   length(wswath) > 0)
%                 Wswath_width(j,k,t) =  dx * ( wswath(end) - wswath(1));   %km
% 
%                 %%% do ca vertical mass flux per slice along the SDC:
%                 %                      %vcaminswath = mean( A_VMca_vcs(wswath(1):wswath(end), k, t, j) , 'omitnan' ) ;
%                 %                      vcaminswath = mean( VMca_3Dswath(wswath(1):wswath(end), k, t, j) , 'omitnan' ) ;
%                 %                      VWcaflux_perz(k)   =   vcaminswath  .*   (dx*1000)   .*   ( dx*1000 * ( wswath(end) - wswath(1))) ;  %  = <W*rhom> * swath width * model dy  (where <> denote mean within the swath at current j,k,t)
% 
%                 %%% do full vertical mass flux per slice along the SDC:
%                 %vminswath = mean( A_VM_vcs(wswath(1):wswath(end), k, t, j) , 'omitnan' ) ;
%                 vminswath = mean( VM_3Dswath(wswath(1):wswath(end), k, t, j) , 'omitnan' ) ;
%                 VWflux_perz(k)   =   vminswath  .*   (dx*1000)   .*   ( dx*1000 * ( wswath(end) - wswath(1))) ;  %  = <W*rhom> * swath width * model dy  (where <> denote mean within the swath at current j,k,t)
%             end
%         end
% 
%         %finish the vert mass flux via summing throughout depth of the current j,t slice.
%         %             WMcaflux_sum(j,t) = sum(VWcaflux_perz,'omitnan');
%         WMflux_sum(j,t) = sum(VWflux_perz,'omitnan');
%         clear VWflux_perz VWcaflux_perz
% 
%         % find depth of w column along j and t:
%         wvobj = double(W_object(:, :, t, j)) ;
%         wvobj( find(wvobj == 0) ) = nan;
% 
% 
%         check = find(isnan(wvobj)); 
%         if( isempty(check)==0  )  %if there is a w obj at this j,t
% 
%             aslswath = A_Ht_asl_vcs(:, :, j) .* wvobj ; % height asl of the current slice of the w swath
%             aglswath = A_AGL_vcs(:, :, j) .* wvobj ; % height asl of the current slice of the w swath
% 
% 	    check2 = find(isnan(aslswath)==0 );
% 	    if( isempty(check2) == 0 )
%             	
% 		swathtop = max( max( aslswath ,[],'omitnan') ) ;
%             	swathbot = min( min( aslswath ,[],'omitnan') ) ;
% 
%                 Wswath_depth(j,t) = (swathtop - swathbot)/1000;   % km
%                 Wswath_top_asl(j,t) = swathtop/1000;
%                 Wswath_bot_asl(j,t) = swathbot/1000;
% 
%                 %log the MULFC_asl/agl @ top_asl/agl because you later find
%                 %that ratio of Wtopasl/agl :  LFCasl/agl yields values >>
%                 %1. Hoping using this lfc in the top/lfc ratio fixes it!
%                 [itop, ktop] = find(aslswath == swathtop) ;
%                 itop = itop(1); ktop = ktop(1);
%                 Wswath_toplfcasl(j,t) = MULFC_all(itop,j,t);
%                 clear itop ktop
% 
%                 swathtopagl = max( max( aglswath ,[],'omitnan') ) ;
%                 swathbotagl = min( min( aglswath ,[],'omitnan') ) ;
%                 %log the MULFC_asl/agl @ top_asl/agl
%                 [itop, ktop] = find(aglswath == swathtopagl) ;
%                 itop = itop(1); ktop = ktop(1);
%                 Wswath_toplfcagl(j,t) = MULFC_agl_all(itop,j,t);
%                 clear itop ktop
% 
%                 %find just the k indices to maybe use later:
%                 [we sd] = find(isnan(aslswath)==0) ;
%                 if(isempty(sd)==0)
%                     maxk = max(sd); maxk = maxk(1);
%                     mink = min(sd); mink = mink(1);
%                     Wswath_ktop_asl(j,t) =  maxk;
%                     Wswath_kbot_asl(j,t) =  mink;
%                 end
% 
% 
%                 if (  isnan(swathtop)==0 & isnan(swathbot)==0  )
%                   % % %calcuate the AGL versions:
%                   %locate the top&bot points along the swath in the domain:
%                     [itop ktop] = find(    swathtop ==  aslswath    ) ; itop= itop(1);
%                     [ibot kbot] = find(    swathbot ==  aslswath    ) ; ibot= ibot(1);
%                     %terr elev at top/bot points:
%                     terrtop = Terr(itop,j);  terrbot = Terr(ibot,j);
%                     %final answer:
%                     Wswath_top_agl(j,t) = (swathtop - terrtop)/1000;  %km
%                     Wswath_bot_agl(j,t) = (swathbot - terrbot)/1000;  %km
%                 end
% 
%             end
%             clear swathtop swathbot itop ktop ibot kbot terrtop terrbot aslswath
% 
% 	end% wvobj nan check
% 
%         % find max of w within column along j and t:
%         wvobj = double(W_object(:, :, t, j)) ;
%         wvobj( find(wvobj == 0) ) = nan;
%         winswath = A_W_vcs(:, :, t, j) .* wvobj ; % height asl of the current slice of the w swath
%         if( isempty(winswath) == 0 )
%             Wswath_maxmag(j,t) = max( max( winswath, [],'omitnan') ) ;
%         end
% 
%         % calculate thermo quantities within each E-W slice of the w swath:
%         Wswath_thetamean(j,t)  = mean(  mean( TH_3Dswath(:,:,t,j) ,'omitnan'),'omitnan' );
%         Wswath_thetaemean(j,t) = mean(  mean( THE_3Dswath(:,:,t,j) ,'omitnan'),'omitnan' );
%         Wswath_vapormean(j,t)  = mean(  mean( QV_3Dswath(:,:,t,j) ,'omitnan'),'omitnan' );
% 
%         %             %%%%%%%%%%%%%%%%%%
%         %             %%%%%% clear-air W swath stuff:
%         %             %%%%%%%%%%%%%%%%%%
%         %
%         %             % find width of clear air w  swath along j and t and k:
%         %             for k = 1:Zvcsa
%         %                 wobj = Wca_object(:, k, t, j) ;
%         %                 wswath = find( wobj == 1 );
%         %                 if( isempty(wswath)==0 & length(wswath) > 0)
%         %                     Wcaswath_width(j,k,t) =  dx * ( wswath(end) - wswath(1));   %km
%         %                 end
%         %             end
%         %
%         %             % j = 100;  %diag tester index
%         %             % find depth of wca column along j and t:
%         %             wvobj = double(Wca_object(:, :, t, j)) ;
%         %             wvobj( find(wvobj == 0) ) = nan;
%         %             aslswath = A_Ht_asl_vcs(:, :, j) .* wvobj ; % height asl of the current slice of the w swath
%         %
%         %             if( isempty(aslswath) == 0 )
%         %                 swathtop = max( max( aslswath ,[],'omitnan') ) ;
%         %                 swathbot = min( min( aslswath ,[],'omitnan') ) ;
%         %
%         %                 Wcaswath_depth(j,t) = (swathtop - swathbot)/1000;   % km
%         %                 Wcaswath_top_asl(j,t) = swathtop/1000;
%         %                 Wcaswath_bot_asl(j,t) = swathbot/1000;
%         %
%         %                 %find just the k indices to maybe use later:
%         %                 [we sd] = find(isnan(aslswath)==0) ;
%         %                 if(isempty(sd)==0)
%         %                     maxk = max(sd); maxk = maxk(1);
%         %                     mink = min(sd); mink = mink(1);
%         %                     Wcaswath_ktop_asl(j,t) =  maxk;
%         %                     Wcaswath_kbot_asl(j,t) =  mink;
%         %                 end
%         %
%         %                 if (  isnan(swathtop)==0 & isnan(swathbot)==0  )
%         %                     % % %calcuate the AGL versions:
%         %                     %locate the top&bot points along the swath in the domain:
%         %                     [itop ktop] = find(    swathtop ==  aslswath    ) ; itop= itop(1);
%         %                     [ibot kbot] = find(    swathbot ==  aslswath    ) ; ibot= ibot(1);
%         %                     %terr elev at top/bot points:
%         %                     terrtop = Terr(itop,j);  terrbot = Terr(ibot,j);
%         %                     %final answer:
%         %                     Wcaswath_top_agl(j,t) = (swathtop - terrtop)/1000;  %km
%         %                     Wcaswath_bot_agl(j,t) = (swathbot - terrbot)/1000;  %km
%         %                 end
%         %
%         %             end
%         %             clear swathtop swathbot itop ktop ibot kbot terrtop terrbot aslswath
%         %
%         %             % find max of w within column along j and t:
%         %             wvobj = double(Wca_object(:, :, t, j)) ;
%         %             wvobj( find(wvobj == 0) ) = nan;
%         %             winswath = A_W_vcs(:, :, t, j) .* wvobj ; % height asl of the current slice of the w swath
%         %             if( isempty(winswath) == 0 )
%         %                 Wcaswath_maxmag(j,t) = max( max( winswath, [],'omitnan') ) ;
%         %             end
%         %
%         %             % calculate clear-air thermo quantities within each E-W slice of the w swath:
%         %             Wcaswath_thetamean(j,t) = mean(  mean( THca_3Dswath(:,:,t,j) ,'omitnan'),'omitnan' );
%         %             Wcaswath_thetaemean(j,t) = mean(  mean( THEca_3Dswath(:,:,t,j) ,'omitnan'),'omitnan' );
%         %             Wcaswath_vapormean(j,t) = mean(  mean( QVca_3Dswath(:,:,t,j) ,'omitnan'),'omitnan' );
% 
%     end
% end
% toc
% disp('     ')
% 
% clear VM_3Dswath winswath TH_3Dswath THE_3Dswath QV_3Dswath A_W_vcs
% 
% %clear  THca_3Dswath QVca_3Dswath THEca_3Dswath winswath we sd aslswath
% 
% % % diagnostic plots:
% %    figure; plot(WMflux_sum(:,2),'ok')
% 
% % % diagnositc plots:
% % figure; contourf(Wswath_width(:,:,15),20)
% % figure; plot(Wcaswath_vapormean(:,15)*1000)
% 
% save('/global/cfs/projectdirs/m1657/jmarquis/LASSO/matlab/checkpt2.mat','Wswath_vapormean','Wswath_maxmag','WMflux_sum','Wswath_width')
% 
% % do some along-SDC smoothing of the w-swath metrics:
% disp('    ')
% disp('  smoothing W swath along SDC  ')
% disp('    ')
% tic
% 
% %1-pass movmean::
% 
% %kinematics:
% Wswath_depth_sm1 = Wswath_depth;   Wswath_depth_sm1(:) = NaN;  dummy = Wswath_depth_sm1;
% %Wcaswath_depth_sm1 = dummy;
% Wswath_maxmag_sm1 = dummy;
% %Wcaswath_maxmag_sm1 = dummy;
% Wswath_width_sm1 = dummy;
% %Wcaswath_width_sm1 = dummy;
% Wswath_VMcasum_sm1 = dummy;
% Wswath_VMsum_sm1 = dummy;
% %thermo:
% Wswath_thetamean_sm1 = dummy;
% Wswath_thetaemean_sm1 = dummy;
% Wswath_vapormean_sm1 = dummy;
% %Wcaswath_thetamean_sm1 = dummy;
% %Wcaswath_thetaemean_sm1 = dummy;
% %Wcaswath_vapormean_sm1 = dummy;
% 
% 
% WIN = 9;  %movmean window length
% 
% for t = 1:1:Tvcsa
% 
%     % %swath depth:
%     Wswath_depth_sm1(:,t)    = movmean(Wswath_depth(:,t),WIN,'omitnan');
%     %Wcaswath_depth_sm1(:,t)  = movmean(Wcaswath_depth(:,t),WIN,'omitnan');
%     % %swath width:
%     ws = mean(Wswath_width(:,:,t),2,'omitnan');
%     %wsca = mean(Wcaswath_width(:,:,t),2,'omitnan');
%     Wswath_width_sm1(:,t)    = movmean(ws,WIN,1,'omitnan');
%     %Wcaswath_width_sm1(:,t)  = movmean(wsca,WIN,1,'omitnan');
%     % %swath strength:
%     Wswath_maxmag_sm1(:,t)   = movmean(Wswath_maxmag(:,t),WIN,'omitnan');
%     %Wcaswath_maxmag_sm1(:,t) = movmean(Wcaswath_maxmag(:,t),WIN,'omitnan');
%     % % vert mass flux
%     Wswath_VMsum_sm1(:,t)   = movmean(WMflux_sum(:,t),WIN,'omitnan');
%     %Wswath_VMcasum_sm1(:,t)   = movmean(WMcaflux_sum(:,t),WIN,'omitnan');
% 
%     % %thermo:
%     Wswath_thetamean_sm1(:,t) = movmean(Wswath_thetamean(:,t),WIN,'omitnan');
%     Wswath_thetaemean_sm1(:,t) = movmean(Wswath_thetaemean(:,t),WIN,'omitnan');
%     Wswath_vapormean_sm1(:,t) = movmean(Wswath_vapormean(:,t),WIN,'omitnan');
%     %         Wcaswath_thetamean_sm1(:,t) = movmean(Wcaswath_thetamean(:,t),WIN,'omitnan');
%     %         Wcaswath_thetaemean_sm1(:,t) = movmean(Wcaswath_thetaemean(:,t),WIN,'omitnan');
%     %         Wcaswath_vapormean_sm1(:,t) = movmean(Wcaswath_vapormean(:,t),WIN,'omitnan');
% 
% end
% clear WMflux_sum
% 
% save('/global/cfs/projectdirs/m1657/jmarquis/LASSO/matlab/checkpt3.mat','Wswath_vapormean_sm1','Wswath_width_sm1','Wswath_maxmag_sm1','Wswath_VMsum_sm1')
% 
% % % 2-pass movmean
% Wswath_depth_sm2 = dummy;
% %Wcaswath_depth_sm2 = dummy;
% Wswath_width_sm2 = dummy;
% %Wcaswath_width_sm2 = dummy;
% Wswath_maxmag_sm2 = dummy;
% %Wcaswath_maxmag_sm2 = dummy;
% %Wswath_VMcasum_sm2 = dummy;
% Wswath_VMsum_sm2 = dummy;
% 
% %     Wcaswath_thetamean_sm2 = dummy;
% %     Wcaswath_thetaemean_sm2 = dummy;
% %     Wcaswath_vapormean_sm2 = dummy;
% Wswath_thetamean_sm2 = dummy;
% Wswath_thetaemean_sm2 = dummy;
% Wswath_vapormean_sm2 = dummy;
% 
% for t = 1:1:Tvcsa
%     Wswath_depth_sm2(:,t)    = movmean(Wswath_depth_sm1(:,t),WIN,'omitnan');
%     %Wcaswath_depth_sm2(:,t)  = movmean(Wcaswath_depth_sm1(:,t),WIN,'omitnan');
% 
%     Wswath_width_sm2(:,t)    = movmean(Wswath_width_sm1(:,t),WIN,'omitnan');
%     %Wcaswath_width_sm2(:,t)  = movmean(Wcaswath_width_sm1(:,t),WIN,'omitnan');
% 
%     Wswath_maxmag_sm2(:,t)   = movmean(Wswath_maxmag_sm1(:,t),WIN,'omitnan');
%     %Wcaswath_maxmag_sm2(:,t) = movmean(Wcaswath_maxmag_sm1(:,t),WIN,'omitnan');
% 
%     %Wswath_VMcasum_sm2(:,t)   = movmean(Wswath_VMcasum_sm1(:,t),WIN,'omitnan');
%     Wswath_VMsum_sm2(:,t)   = movmean(Wswath_VMsum_sm1(:,t),WIN,'omitnan');
% 
%     %         Wcaswath_thetamean_sm2(:,t) = movmean(Wcaswath_thetamean_sm1(:,t),WIN,'omitnan');
%     %         Wcaswath_thetaemean_sm2(:,t) = movmean(Wcaswath_thetaemean_sm1(:,t),WIN,'omitnan');
%     %         Wcaswath_vapormean_sm2(:,t) = movmean(Wcaswath_vapormean_sm1(:,t),WIN,'omitnan');
%     Wswath_thetamean_sm2(:,t) = movmean(Wswath_thetamean_sm1(:,t),WIN,'omitnan');
%     Wswath_thetaemean_sm2(:,t) = movmean(Wswath_thetaemean_sm1(:,t),WIN,'omitnan');
%     Wswath_vapormean_sm2(:,t) = movmean(Wswath_vapormean_sm1(:,t),WIN,'omitnan');
% 
% end
% toc
% disp('   ')
% 
% 
% 
% save('/global/cfs/projectdirs/m1657/jmarquis/LASSO/matlab/checkpt4.mat')
% 
% 
% % % diagnostic plots:
% % figure; plot(Wcaswath_vapormean(:,15)*1000,'k'); hold on; plot(Wcaswath_vapormean_sm1(:,15)*1000,'r'); plot(Wcaswath_vapormean_sm2(:,15)*1000,'b');
% % figure; plot( WMflux_sum(:,1),'*k' ); hold on; plot( Wswath_VMcasum_sm1(:,1),'or' ); hold on;  plot( Wswath_VMcasum_sm2(:,1),'dg' )
% 
% 
% 
% % works to here
% 
% 
% 
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    %%%%% catalog in-swath 2D sounding metrics as function of y,t:
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     size(W_object);
% 
%     W_object_vertcomp = max( W_object, [], 2 );   size(W_object_vertcomp);
%     W_object_vertcomp  = permute(W_object_vertcomp, [1 4 3 2]);   size(W_object_vertcomp);
% 
%     Wswath_mulfcagl = zeros(Yvcsa,Tvcsa) ;    Wswath_mulfcagl(:) = NaN; %JNM-27Nov
%     Wswath_mulfcasl = zeros(Yvcsa,Tvcsa) ;    Wswath_mulfcasl(:) = NaN; %JNM-27Nov
%     Wswath_mucin = zeros(Yvcsa,Tvcsa)  ;   Wswath_mucin(:) = NaN;
%     Wswath_mucape = zeros(Yvcsa,Tvcsa)  ;  Wswath_mucape(:) = NaN;
%     %   swath_j = zeros(Yvcsa,Tvcsa)  ;   swath_j(:) = NaN;
% 
%     for j = 1:Yvcsa
%         for t = 1:Tvcsa
%             %   j = 250; t = 10;   %diagnostic coordinates
% 
%             %locate the swath in this j:
%             swathi = find(W_object_vertcomp(:,j,t)) ;
% 
%             Wswath_mulfcagl(j,t)   =  mean( MULFC_agl_all(swathi,j,t) ,'omitnan')/1000  ;  %km  
%             Wswath_mulfcasl(j,t)   =  mean( MULFC_all(swathi,j,t) ,'omitnan')/1000  ;  %km      
%             Wswath_mucin(j,t)   =  mean( MUCIN_all(swathi,j,t) ,'omitnan')  ;
%             Wswath_mucape(j,t)  =  mean( MUCAPE_all(swathi,j,t) ,'omitnan')  ;
% 
%         end
%     end
% 
%     clear MULFC_agl_all MULFC_all MUCIN_all MUCAPE_all
% 
% %     %clear Wca_object_vertcomp swathi
% %
% %
% %     % works to here
% %
% % %     %quick clear-out for memory:
% % %     clear Conv_3tmean Conv_5tmean Conv_9tmean Conv_alltmean Conv_mean_nearsfc Conv_raw_nearsfc Conv_mean_9jmean_9tmean Conv_mean
% % %     clear Ht_agl Ht_asl abovegrnd A_THETA_vcs A_THETAE_vcs  A_QVAP_vcs alltmean  % A_V_vcs
% % %     %clear QVAP_all QVca_all RH_all RHca_all sU_ASL sV_ASL sW_ASL THca_all THEca_all THETA_all THETAE_all U_all V_all W_all Qcld_all Wcld_all
% % %     clear Wca_all W_bin Wca_bin X Y Wca_3Dmask W_3Dmask
% % %     clear THEca_3Dswath THE_3Dswath TH_3Dswath THca_3Dswath sAGLinASL
% % %     clear THca_ASL2 THEca_ASL2 THETA_ASL2 THETAE_ASL2 sConv_ASL RH_ASL2 RHca_ASL2 QVca_3Dswath QV_3Dswath QVAP_ASL2 QVca_ASL2
% % %     clear dum  Conv_bin
% %
% %
% %
% % %     %%% SDC BOX HERE:    ASL version
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %     %%%%% catalog above and downstream-swath env metric profiles as function of y,z,t:
% % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %
% % %     size(A_RHca_vcs)   ;
% % %
% % %     %layer (bot/top hts asl) you want to analyze
% % %     bot = find( newz_asl == 2500 );  % meters
% % %     top = find( newz_asl == 5000 );
% % %
% % %     swathboxaloft_RH = zeros(Yvcsa,length(newz_asl(bot:top)),Tvcsa) ;    swathboxaloft_RH(:) = NaN;
% % %     swathboxaloft_Zasl = newz_asl(bot:top);
% % %
% % %     %     swathaloft_mucin = zeros(Yvcsa,Tvcsa)  ;   swathaloft_mucin(:) = NaN;
% % %     %     swathaloft_mucape = zeros(Yvcsa,Tvcsa)  ;  swathaloft_mucape(:) = NaN;
% % %
% % %     tic
% % %     for t = 1:Tvcsa
% % %         for j = 1:Yvcsa
% % %             %for k = bot:top
% % %             %   j = 250; t = 1;   %diagnostic coordinates
% % %
% % %             %locate the swath in this j:
% % %             swathi = find(Wca_object_vertcomp(:,j,t)) ;
% % %
% % %             if(isempty(swathi)==0)
% % %                 % define box width in e-w direction along each j, going west 5 km from
% % %                 % west edge of swath to 15 km east of east edge of swath
% % %
% % %                 boxwidth = [ swathi(1) - 5/dx : 1 : swathi(end) + 15/dx ] ;
% % %                 swathboxaloft_RH(j,:,t)   =  mean(  A_RHca_vcs( boxwidth, bot:top, t, j ), 1, 'omitnan' )  ;  %percent   "1" is for meaning along the k dimension - checked manually by length(bot:top) - (even tho it feels like it should be a "2")
% % %             end
% % %
% % %             clear boxwidth swathi
% % %
% % %             %end
% % %         end
% % %     end
% % %     toc   %~ 10 sec
% % %
% % %     clear Wca_object_vertcomp swathi
% % %
% %
% %     % diagnostic plot:
% %     % figure; plot(swathboxaloft_RH(:,2,1),'ok')
% %
% %
% %
% %
% %
% %
% 
% 
%     %%% SDC BOX HERE:  AGL FRAME:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%% catalog above and downstream-swath env metric profiles as function of y,z,t:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % size(A_RHca_vcs)   ;
% 
%     %top/bottom of layer (m agl) you want to evaluate:
%     botagl = 2500.;
%     topagl = 5000.;
%     aglspan = [botagl:200:topagl];
% 
%     Wswathboxaloft_RH = zeros(Yvcsa,length(aglspan),Tvcsa) ;    Wswathboxaloft_RH(:) = NaN;
%     Wswathboxaloft_Zasl = aglspan;
% 
%     %swathboxsubcloud_U = zeros(Yvcsa,length(aglspan),Tvcsa) ;    swathboxaloft_RH(:) = NaN;
% 
%     %     swathaloft_mucin = zeros(Yvcsa,Tvcsa)  ;   swathaloft_mucin(:) = NaN;
%     %     swathaloft_mucape = zeros(Yvcsa,Tvcsa)  ;  swathaloft_mucape(:) = NaN;
%     %tmpnewzasl = newz_asl';  tmpnewzasl(tmpnewzasl == -999)=NaN;
%     tic
%     for t = 1:Tvcsa
%         for j = 1:Yvcsa
%             %for k = bot:top
%             %   j = 250; t = 1;   %diagnostic coordinates
% 
%             %locate the swath in this j:
%             swathi = find(W_object_vertcomp(:,j,t)) ;  %i span of swath width in LASSO frame
% 
%             if(isempty(swathi)==0)
%                 % define box width in e-w direction along each j, going west 5 km from
%                 % west edge of swath to 15 km east of east edge of swath
% 
%                 ileft = swathi(1) - 5/dx;
%                 iright = swathi(end) + 15/dx;
%                 
%                 %if box is too wide for the domain, set it to domain edge:
%                 if(ileft < 1);  ileft =1;  end
%                 if(iright > Xvcsa);  iright = Xvcsa;  end                
%                 
%                 boxwidth = [ ileft : 1 : iright ] ;
% 
%                 RHprofs= zeros( length(boxwidth), length(aglspan) );
%                 ic = 0;
%                 for i =  boxwidth(1) : boxwidth(end)
%                     ic = ic+1;
% 
% %                     bot = find( A_AGL_vcs(i,:,j) <= botagl );     bot = bot(end);
% %                     top = find( A_AGL_vcs(i,:,j) >= topagl );     top = top(1);     %  A_AGL_vcs(i,top,j)
% 
%                     %gerenrate agl profile by interp asl to steady prescribed AGL layer:
%                     RHprofs(ic,:) = interp1( newz_asl',    A_RHca_vcs( i, :, t, j ),   aglspan  ) ;
% 
%                 end %i
% 
%                 Wswathboxaloft_RH(j,:,t)   =  mean(  RHprofs, 1, 'omitnan' )  ;  %percent   "1" is for meaning along the k dimension - checked manually by length(bot:top) - (even tho it feels like it should be a "2")
% 
%             end %j
% 
%             clear boxwidth swathi
% 
%         end %
%     end
%     toc   %~ 10 sec
% 
%     % manual QC:
%     Wswathboxaloft_RH(Wswathboxaloft_RH < 0) = NaN;
% 
%     %  figure; contourf(swathboxaloft_RH(:,:,27),20)
% 
%     clear  swathi  RHprofs top bot boxwidth
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     %%% SDC BOX HERE:  AGL FRAME:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%% catalog above and downstream-swath env wind profiles as function of y,local-swath-bottom to top of domain (asl),t:
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [as df gh jk] = size(A_U_vcs) ;
%     Wswathbox_Uswathdepth = zeros(Yvcsa,df,Tvcsa);   Wswathbox_Uswathdepth(:) = NaN;
%     Wswathbox_Vswathdepth = zeros(Yvcsa,df,Tvcsa);   Wswathbox_Vswathdepth(:) = NaN;
% 
%     tic
%     for t = 1:Tvcsa
%         for j = 1:Yvcsa
% 
% 
%             % t = 10;  j = 150;  %tester set
% 
%             %locate the swath in this j:
%             swathi = find(W_object_vertcomp(:,j,t)) ;  %i span of swath width in LASSO frame
% 
%             if(isempty(swathi)==0)
% 
%                 % define box width in e-w direction along each j, going west 5 km from
%                 % west edge of swath to 15 km east of east edge of swath
% 
%                 ileft = swathi(1) - 5/dx;
%                 iright = swathi(end) + 15/dx;
%                 
%                 %if box is too wide for the domain, set it to domain edge:
%                 if(ileft < 1);  ileft =1;  end
%                 if(iright > Xvcsa);  iright = Xvcsa;  end                
%                 
%                 boxwidth = [ ileft : 1 : iright ] ;
%                 
%                 
%                 %depth of wswath (num of vert indices on ASL grid)
%                 dep = length([Wswath_kbot_asl(j,t):df]) ;
% 
%                 Uprofs = zeros( length(boxwidth),  dep ) ;
%                 Vprofs = zeros( length(boxwidth),  dep ) ;
%                 ic = 0;
%                 for i =  boxwidth(1) : boxwidth(end)
%                     ic = ic+1;
%                     Uprofs(ic,:) = A_U_vcs( i, Wswath_kbot_asl(j,t):df, t, j ) ;
%                     Vprofs(ic,:) = A_V_vcs( i, Wswath_kbot_asl(j,t):df, t, j ) ;
%                 end %i
%                 Uprofs(Uprofs == -999) = NaN;
%                 Vprofs(Vprofs == -999) = NaN;
%                 Wswathbox_Uswathdepth(j,(df-dep)+1:df,t)  =  mean(  Uprofs, 1, 'omitnan' )  ;
%                 Wswathbox_Vswathdepth(j,(df-dep)+1:df,t)  =  mean(  Vprofs, 1, 'omitnan' )  ;
% %                  swathbox_Uswathdepth(j,1:dep,t)   =  mean(  Uprofs, 1, 'omitnan' )  ;     %percent   "1" is for meaning along the k dimension - checked manually by length(bot:top) - (even tho it feels like it should be a "2")
% %                 swathbox_Vswathdepth(j,1:dep,t)   =  mean(  Vprofs, 1, 'omitnan' )  ;
% 
%             end %if
% 
%             clear boxwidth swathi dep
% 
%         end % j
%     end %t
%     toc   %~ 1 sec
% 
%     clear  swathi Uprofs Vprofs top bot boxwidth %Wca_object_vertcomp
%


% figure; contourf(Lon,Lat,Terr,20)
% 
% figure; contourf(Lon,20)

disp('FF')

%%% SDC BOX HERE:  ASL FRAME:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% catalog above and downstream-swath env wind profiles as function of y,local-swath-bottom to top of domain (asl),t:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imax kmax tmax jmax] = size(A_U_vcs) 
westSDCbox_Uprofs  = zeros(jmax,kmax,tmax);   westSDCbox_Uprofs(:) = NaN;
westSDCbox_Vprofs  = zeros(jmax,kmax,tmax);   westSDCbox_Vprofs(:) = NaN;
eastSDCbox_Uprofs  = zeros(jmax,kmax,tmax);   eastSDCbox_Uprofs(:) = NaN;
eastSDCbox_Vprofs  = zeros(jmax,kmax,tmax);   eastSDCbox_Vprofs(:) = NaN;
%newFRN
eastSDCbox_THVprofs  = zeros(jmax,kmax,tmax);   eastSDCbox_THVprofs(:) = NaN;
westSDCbox_THVprofs  = zeros(jmax,kmax,tmax);   eastSDCbox_THVprofs(:) = NaN;

% east and west i index limits to the east/west averging boxes
% westi1 = 60;   westi2= 120;
% easti1= 440;   easti2= 500;

% wloni1 = -65.079;   wloni2 = -65.028; 
% westi1 = find( abs((Lon(:,800) - wloni1)) ==  min(abs((Lon(:,800) - wloni1))) )
% westi2 = find( abs((Lon(:,800) - wloni2)) ==  min(abs((Lon(:,800) - wloni2))) )
% eloni1 = -64.726;   eloni2 = -64.682; 
% easti1 = find( abs((Lon(:,800) - eloni1)) ==  min(abs((Lon(:,800) - eloni1))) )
% easti2 = find( abs((Lon(:,800) - eloni2)) ==  min(abs((Lon(:,800) - eloni2))) )

westi1 = zeros(jmax,1); %easternmost box
westi2 = zeros(jmax,1);
easti1 = zeros(jmax,1);
easti2 = zeros(jmax,1);
for j = 1:jmax
    itermax = find(  Terr(:,j) == max( Terr(:,j) )  ) ; itermax = itermax(1);
    wi2 = find( Terr(1:itermax,j) <= 1100. ) ;   westi2(j) = wi2(end);        %newFRN
    westi1(j) = westi2(j) ;  %- 20;  %box width of di = 20      %newFRN
    
    ei1 = find( Terr(itermax:end,j) >= 1100. ) ;   easti1(j) = ei1(end) + itermax;     %newFRN
    easti2(j) = easti1(j) ;%  + 20;                             %newFRN
end
easti2(easti2>imax)=imax;
westi1(westi1<1)=1;

% figure; contourf(Terr',20); hold on; plot(westi1,[1:1600],'ow'); plot(westi2,[1:1600],'or');  plot(easti1,[1:1600],'mo'); plot(easti2,[1:1600],'og')


A_U_vcs(A_U_vcs == -999) = NaN;
A_V_vcs(A_V_vcs == -999) = NaN;
A_THETAv_vcs(A_THETAv_vcs == -999) = NaN;   %newFRN

%tabulate profiles of i-horizontal mean asl wind profiles
tic
for t = 1:tmax
    for j = 1:jmax
      
        westSDCbox_Uprofs(j,:,t)  =  permute(mean(A_U_vcs(westi1(j):westi2(j),:,t,j),1,'omitnan'),[2,1]);
        westSDCbox_Vprofs(j,:,t)  =  permute(mean(A_V_vcs(westi1(j):westi2(j),:,t,j),1,'omitnan'),[2,1]);

        eastSDCbox_Uprofs(j,:,t)  =  permute(mean(A_U_vcs(easti1(j):easti2(j),:,t,j),1,'omitnan'),[2,1]);
        eastSDCbox_Vprofs(j,:,t)  =  permute(mean(A_V_vcs(easti1(j):easti2(j),:,t,j),1,'omitnan'),[2,1]);

        %newFRN
        eastSDCbox_THVprofs(j,:,t)  = permute(mean(A_THETAv_vcs(easti1(j):easti2(j),:,t,j),1,'omitnan'),[2,1]);
        westSDCbox_THVprofs(j,:,t)  = permute(mean(A_THETAv_vcs(westi1(j):westi2(j),:,t,j),1,'omitnan'),[2,1]);

    end % j
end %t
toc   %~ 1 sec





% %         %%% SDC BOX HERE:  AGL FRAME:
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%% catalog above and downstream-swath env wind profiles as function of y,z,t:
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% %
% %     %top/bottom of layer (m agl) you want to evaluate:
% %     botagl = 10.;
% %     topagl = 3010.;
% %     aglspan = [botagl:200:topagl];
% %
% %     swathboxLL_U = zeros(Yvcsa,length(aglspan),Tvcsa) ;    swathboxLL_U(:) = NaN;
% %     swathboxLL_V = zeros(Yvcsa,length(aglspan),Tvcsa) ;    swathboxLL_V(:) = NaN;
% %     swathboxLL_Zasl = aglspan;
% %
% %     tic
% %     for t = 1:Tvcsa
% %         for j = 1:Yvcsa
% %             %for k = bot:top
% %             %   j = 250; t = 1;   %diagnostic coordinates
% %
% %             %locate the swath in this j:
% %             swathi = find(Wca_object_vertcomp(:,j,t)) ;  %i span of swath width in LASSO frame
% %
% %             if(isempty(swathi)==0)
% %                 % define box width in e-w direction along each j, going west 5 km from
% %                 % west edge of swath to 15 km east of east edge of swath
% %
% %                 boxwidth = [ swathi(1) - 5/dx : 1 : swathi(end) + 15/dx ] ;
% %
% %                 Uprofs= zeros( length(boxwidth), length(aglspan) );
% %                 Vprofs= zeros( length(boxwidth), length(aglspan) );
% %                 ic = 0;
% %                 for i =  boxwidth(1) : boxwidth(end)
% %                     ic = ic+1;
% %
% %                     %gerenrate agl profile by interp asl to steady prescribed AGL layer:
% %                     Uprofs(ic,:) = interp1( newz_asl',    A_U_vcs( i, :, t, j ),   aglspan  ) ;
% %                     Vprofs(ic,:) = interp1( newz_asl',    A_V_vcs( i, :, t, j ),   aglspan  ) ;
% %
% %                 end %i
% %
% %                 swathboxLL_U(j,:,t)   =  mean(  Uprofs, 1, 'omitnan' )  ;  %percent   "1" is for meaning along the k dimension - checked manually by length(bot:top) - (even tho it feels like it should be a "2")
% %                 swathboxLL_V(j,:,t)   =  mean(  Vprofs, 1, 'omitnan' )  ;
% %
% %             end %if
% %
% %             clear boxwidth swathi
% %
% %         end % j
% %     end %t
% %
% %     toc   %~ 10 sec
% %
% %
% %     %  figure; contourf(swathboxaloft_RH(:,:,27),20)
% %
% %     clear Wca_object_vertcomp swathi Uprofs Vprofs top bot boxwidth



%newFRN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  v13c-ii:   dz(j,t) =  U>?(j,t) -  peak(j) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dzpeakU5  = zeros(Ys,Ts);    dzpeakU5(:) = NaN;
dzpeakU10 = zeros(Ys,Ts);   dzpeakU10(:) = NaN;
dzpeakU15 = zeros(Ys,Ts);   dzpeakU15(:) = NaN;
dzpeakU20 = zeros(Ys,Ts);   dzpeakU20(:) = NaN;
dzpeakU8 = zeros(Ys,Ts);   dzpeakU8(:) = NaN;
ridge  = max(Terr,[],1,'omitnan') ;
Iridge = zeros(Ys,1);
Z5= zeros(Ys,1);  Z5(:) = NaN;
Z8= zeros(Ys,1); Z8(:) = NaN;
Z10= zeros(Ys,1); Z10(:) = NaN;
Z15= zeros(Ys,1); Z15(:) = NaN;
Z20= zeros(Ys,1); Z20(:) = NaN;
for t = 1:Ts
    for j = 1:Ys
       %  t = 16;  j = 500;
       iridge = find( Terr(:,j)  == ridge(j)  ) ;    iridge = iridge(1) ;
       Iridge(j) = iridge ;
       ku = find(  A_U_vcs(iridge,:,t,j) >= 5.0 ) ;
       if( isempty(ku)==0)
           ku = ku(1);  %pick the first height asl
           dzpeakU5(j,t)  =  regzasl(iridge,j,ku) - ridge(j) ;
           Z5(j,t) = regzasl(iridge,j,ku);
       end
       ku = find(  A_U_vcs(iridge,:,t,j) >= 10.0 ) ;
       if( isempty(ku)==0)
           ku = ku(1);  %pick the first height asl
           dzpeakU10(j,t)  =  regzasl(iridge,j,ku) - ridge(j) ;
           Z10(j,t) = regzasl(iridge,j,ku);
       end
       ku = find(  A_U_vcs(iridge,:,t,j) >= 7.5 ) ;
       if( isempty(ku)==0)
           ku = ku(1);  %pick the first height asl
           dzpeakU8(j,t)  =  regzasl(iridge,j,ku) - ridge(j) ;
           Z8(j,t) = regzasl(iridge,j,ku);
       end
       ku = find(  A_U_vcs(iridge,:,t,j) >= 15.0 ) ;
       if( isempty(ku)==0)
           ku = ku(1);  %pick the first height asl
           dzpeakU15(j,t)  =  regzasl(iridge,j,ku) - ridge(j) ;
           Z15(j,t) = regzasl(iridge,j,ku);
       end
       ku = find(  A_U_vcs(iridge,:,t,j) >= 20.0 ) ;
       if( isempty(ku)==0)
           ku = ku(1);  %pick the first height asl
           dzpeakU20(j,t)  =  regzasl(iridge,j,ku) - ridge(j) ;
           Z20(j,t) = regzasl(iridge,j,ku);
       end
    end
end
%newFRN  end
%figure; plot(Z8(:,16),'o');  axis([1 1600 1000 10000])



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%  v13c-i:  convert U,V,ThetaV as f(j,k,t) to a Mtn FroudeNumber; 
% use markowski formula FRM = Uo/N*HM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rzasl = permute(regzasl(1,1,:),[3 1 2] );
dzasl = gradient(rzasl);


% %newFRN
HM = ridge - 1100  ; 
% find upslope e/w between terr=1100km and local (along-i) peak:
k1100 = find(rzasl == 1100) ;
for j = 1:Ys
    kpeak(j) = find( abs(ridge(j) - rzasl) == min(abs(ridge(j) - rzasl)) ) ;
end


Nwest = zeros(Ys,Zvcsa,Ts);  Nwest(:) = NaN;
west_N_mean = zeros(Ys,Ts);  west_N_mean(:) = NaN;
west_upslope_mean = zeros(Ys,Ts);   west_upslope_mean(:) = NaN;
FroudeW = zeros(Ys,Ts);   FroudeW(:) = NaN;
west_upslope_depth = zeros(Ys,Ts);   west_upslope_depth(:) = NaN;          %NEWD

Neast = zeros(Ys,Zvcsa,Ts);  Neast(:) = NaN;
east_N_mean = zeros(Ys,Ts);  east_N_mean(:) = NaN;
east_upslope_mean = zeros(Ys,Ts);   east_upslope_mean(:) = NaN;
FroudeE = zeros(Ys,Ts);   FroudeE(:) = NaN;
east_upslope_depth = zeros(Ys,Ts);   east_upslope_depth(:) = NaN;          %NEWD

%{
%level by level method that yields imaginary N: 
for t = 1:Ts
    for j = 1:Ys
        %  t = 15; j = 500;
        eTHVprof = permute(eastSDCbox_THVprofs(j,:,t),[2 1 3]) ; 
        dTedz = gradient(eTHVprof) ./ dzasl ; 
        wTHVprof = permute(westSDCbox_THVprofs(j,:,t),[2 1 3]) ;
        dTwdz = gradient(wTHVprof) ./ dzasl ;
        blah1 = ( 9.81 ./ eTHVprof) ;
        Neast(j,:,t) = (   blah1 .* dTedz   ).^0.5 ;
        blah1 = ( 9.81 ./ wTHVprof) ;
        Nwest(j,:,t) = (   blah1 .* dTwdz   ).^0.5 ;
    end
end
%}

%N calcuated by discretizing over dz = 1100 to ridge:
Neast = zeros(Ys,Ts);  Neast(:) = NaN;
Nwest = zeros(Ys,Ts);  Nwest(:) = NaN;
for t = 1:Ts
    for j = 1:Ys
        %  t = 12; j = 171;
        eTHVprof = permute(eastSDCbox_THVprofs(j,:,t),[2 1 3]) ; 
        pp = find( isnan(eTHVprof)==0 ) ;
        dTedz =  ( eTHVprof( kpeak(j) ) -  eTHVprof( pp(1) ) ) ./ HM(j) ; 
       
        wTHVprof = permute(westSDCbox_THVprofs(j,:,t),[2 1 3]) ; 
        pp = find( isnan(wTHVprof)==0 ) ;
        dTwdz =  ( wTHVprof( kpeak(j) ) -  wTHVprof( pp(1) ) ) ./ HM(j) ; 

        blah1 = ( 9.81 ./ mean(eTHVprof(1:kpeak(j)),'omitnan') ) ;
        Neast(j,t) = (   blah1 .* dTedz   ).^0.5 ;
       
        blah1 = ( 9.81 ./ mean(wTHVprof(1:kpeak(j)),'omitnan') ) ;
        Nwest(j,t) = (   blah1 .* dTwdz   ).^0.5 ;

        if(dTedz <= 0)
            Neast(j,t) = 9.99;
        end
        if(dTwdz <= 0)
            Nwest(j,t) = 9.99;
        end
    end
end

% find top of upslope layer or kpeak (whichever lower) WEST SIDE
for t = 1:Ts
    for j = 1:Ys
        %   t = 15; j = 369;
        west_u = permute(westSDCbox_Uprofs(j,:,t), [2 3 1]) ; %west side profile
        %west_u = west_u(k1100:kpeak(j)) ;  %cut to just sub local terr peak
        west_u = west_u(1:kpeak(j)) ;  %cut to just sub local terr peak
        kwtop = kpeak(j) ;  %start with terr peak height
        for k = 2:length(west_u)  %loop to lower the top k if it stops being upslope
            if(  west_u(k) < 0   )
                kwtop = k-1 ;
                break
            end
        end
        west_u = west_u(1:kwtop)  ;
        west_upslope_mean(j,t) = mean( west_u,'omitnan') ;
        %NEWD
        west_upslope_depth(j,t) =  abs( rzasl( kwtop )  -  rzasl( k1100 ) ) ./ ( rzasl(kpeak(j)) - rzasl( k1100 ) )   ;  %abs() because some discretized heights are the next one below 1100	
    end
end

% find top of upslope layer or kpeak (whichever lower) EAST SIDE
for t = 1:Ts
    for j = 1:Ys
        %   t = 13; j = 1308;
        east_u = permute(eastSDCbox_Uprofs(j,:,t), [2 3 1]) ; %east side profile
        east_u = east_u(1:kpeak(j)) ;  %cut to just sub local terr peak
        ketop = kpeak(j) ;  %start with terr peak height
        for k = 2:length(east_u)  %loop to lower the top k if it stops being upslope
            if(  east_u(k) > 0   )
                ketop = k-1 ;
                break
            end
        end
        east_u = east_u(1:ketop)  ;
        east_upslope_mean(j,t) = mean( east_u,'omitnan') ;
        east_upslope_depth(j,t) =  abs( rzasl( ketop )  -  rzasl( k1100 ) ) ./ ( rzasl(kpeak(j)) - rzasl( k1100 ) )  ; %NEWD
    end
end
east_upslope_mean = abs(east_upslope_mean);

%evaluate a N (mean from ground to kpeak) for the FrN calc with 1100-ridge
%discretization
for t = 1:Ts
    for j = 1:Yvcsa
        %  j = 600; t = 15;
        west_N_mean(j,t) = mean(  Nwest(j,t),'omitnan') ;
        east_N_mean(j,t) = mean(  Neast(j,t),'omitnan') ;
    end
end

% west side froude num: 
for j = 1:Yvcsa
    for t = 1:Ts

        %set frn to 1 is neutral/unstable
        if(  Nwest(j,t) == 9.99  &  isnan(west_upslope_mean(j,t)) == 0    )
            FroudeW(j,t) = 999.9 ;
        elseif(  Nwest(j,t) ~= 9.99  )
            FroudeW(j,t) = west_upslope_mean(j,t) / ( HM(j) * west_N_mean(j,t) ) ;
        elseif(  Nwest(j,t) == 9.99 &  isnan(west_upslope_mean(j,t)) )
            FroudeW(j,t) = west_upslope_mean(j,t) / ( HM(j) * west_N_mean(j,t) ) ;
        end

        %set frn to 1 is neutral/unstable
        if(  Neast(j,t) == 9.99  &  isnan(east_upslope_mean(j,t)) == 0    )
            FroudeE(j,t) = 999.9 ;
        elseif(  Neast(j,t) ~= 9.99  )
            FroudeE(j,t) = east_upslope_mean(j,t) / ( HM(j) * east_N_mean(j,t) ) ;
        elseif(  Neast(j,t) == 9.99 &  isnan(east_upslope_mean(j,t)) )
            FroudeE(j,t) = east_upslope_mean(j,t) / ( HM(j) * east_N_mean(j,t) ) ;
        end

    end
end
FroudeEjt = FroudeE;
FroudeWjt = FroudeW;

%rename:
east_N_mean_Fr = east_N_mean;
west_N_mean_Fr = west_N_mean;
east_upslope_mean_Fr = east_upslope_mean;
west_upslope_mean_Fr = west_upslope_mean;
HM_Fr = HM;
%
dzpeakU5_vent = dzpeakU5;     Z5_vent = Z5;
dzpeakU8_vent = dzpeakU8;     Z8_vent = Z8;
dzpeakU10_vent = dzpeakU10;   Z10_vent = Z10;
dzpeakU15_vent = dzpeakU15;   Z15_vent = Z15;
dzpeakU20_vent = dzpeakU20;   Z20_vent = Z20;




%diagnostic plots
%{
% mean of Fr with stable conditions and upslope flow found from ground up:
FRE_chop = FroudeEjt;  
FRE_chop(find(FRE_chop==999.9)) = NaN;
meanFRe = mean( FRE_chop, 1,'omitnan'  ) ;
% find fraction of unstable/neutral stable conditions & upslope: 
frac_stableFRE = zeros(Ts,1) ; frac_stableFRE(:) = NaN; 
for t = 1:Ts
    frac_stableFRE(t) =  length(find(  FroudeEjt(:,t) == 999.9   )) / Yvcsa ;
end

FRW_chop = FroudeWjt;  
FRW_chop(find(FRW_chop==999.9)) = NaN;
meanFRw = mean( FRW_chop, 1,'omitnan'  ) ;
% find fraction of unstable/neutral stable conditions & upslope: 
frac_stableFRW = zeros(Ts,1) ; frac_stableFRW(:) = NaN; 
for t = 1:Ts
    frac_stableFRW(t) =  length(find(  FroudeWjt(:,t) == 999.9   )) / Yvcsa ;
end

%plot them
figure;   
yyaxis right
plot(meanFRe,'r');   hold on;   plot(meanFRw,'r--'); 
yyaxis left
plot(frac_stableFRE,'>b');  hold on;   plot(frac_stableFRW,'<b');
%}




%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Try a very along-SDC bulk version:   INCOMPLETE!
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


bulk_west_N_mean = zeros(Ts,1);        bulk_west_N_mean(:) = NaN;
bulk_west_upslope_mean = zeros(Ts,1);  bulk_west_upslope_mean(:) = NaN;
bulk_FroudeW = zeros(Ts,1);            bulk_FroudeW(:) = NaN;

% meankpeak = mean(kpeak); meankpeak = floor(meankpeak);
% HM = rzasl(meankpeak) - 1100 ;
%     eTHVprof = permute(eastSDCbox_THVprofs(:,:,t),[2 1 3]) ;
%     eTHVprof = mean(eTHVprof,2,'omitnan') ;
%     pp = find( isnan(eTHVprof)==0 ) ;
%     dTedz =  ( eTHVprof( meankpeak ) -  eTHVprof( pp(1) ) ) ./ HM(j) ;

%N calcuated by discretizing over dz = 1100 to ridge:
Neast = zeros(Ys,Ts);  Neast(:) = NaN;
Nwest = zeros(Ys,Ts);  Nwest(:) = NaN;
for t = 1:Ts
    for j = 1:Ys
        %  t = 12; j = 171;
        eTHVprof = permute(eastSDCbox_THVprofs(j,:,t),[2 1 3]) ; 
        pp = find( isnan(eTHVprof)==0 ) ;
        dTedz =  ( eTHVprof( kpeak(j) ) -  eTHVprof( pp(1) ) ) ./ HM(j) ; 
       
        wTHVprof = permute(westSDCbox_THVprofs(j,:,t),[2 1 3]) ; 
        pp = find( isnan(wTHVprof)==0 ) ;
        dTwdz =  ( wTHVprof( kpeak(j) ) -  wTHVprof( pp(1) ) ) ./ HM(j) ; 

        blah1 = ( 9.81 ./ mean(eTHVprof(1:kpeak(j)),'omitnan') ) ;
        Neast(j,t) = (   blah1 .* dTedz   ).^0.5 ;
       
        blah1 = ( 9.81 ./ mean(wTHVprof(1:kpeak(j)),'omitnan') ) ;
        Nwest(j,t) = (   blah1 .* dTwdz   ).^0.5 ;

        if(dTedz <= 0)
            Neast(j,t) = 9.99;
        end
        if(dTwdz <= 0)
            Nwest(j,t) = 9.99;
        end
    end
end
%}

%when you did N layer by layer and got imaginary N's: 
%{
%evaluate a N (mean from ground to kpeak) for the FrN calc
for t = 1:Ts
    for j = 1:Yvcsa
        %  j = 600; t = 15;
        west_N_mean(j,t) = mean(  Nwest(j,1:kpeak(j),t),'omitnan') ;
    end
end
% west side froude num: 
for j = 1:Yvcsa
    for t = 1:Ts
        FroudeW(j,t) = west_upslope_mean(j,t) / ( HM(j) * west_N_mean(j,t) ) ;
    end
end
%}
%  newFRN  end





%     %%%%%%%%%% plot'ification:
%
%     if(plotme == 1)
%
%         dualpol_colmap
%
%         % diagnositc plots of W swath and refl properites along the SDC :
%         for t = 1:1:Tvcsa
%
%             %    t = 15
%
%             refl = A_DBZ_vcs(:,:,t,:); refl = permute(refl,[4 2 1 3]);
%             refl_maxalongi = max(refl,[],3,'omitnan') ;
%
%             fig = figure('position',[1196,118,1141,729])
%             Termax = max(Terr,[],1);
%
%             subplot(4,1,1)
%             contourf(yasl(1,:,1)',newz_asl/1000',refl_maxalongi',[0:5:70],'LineColor','none')
%             colormap(zhmap)
%             caxis([0 70])
%             hold on
%             plot( yasl(1,:,1),  Termax'/1000,'Color',[0.8 0.6 0],'LineWidth',3)
%             axis([-110 80 0 8])
%             xlabel( 'y (km)' )
%             ylabel( 'Depth of W swath (km)' )
%
%             subplot(4,1,2)
%             hold on
%             %     plot(yasl(1,:,1),Wswath_depth(:,t),'or')
%             %     plot(yasl(1,:,1),Wcaswath_depth(:,t),'og')
%
%             plot(yasl(1,:,1),Wswath_depth(:,t),'r','LineWidth',0.5)
%             plot(yasl(1,:,1),Wcaswath_depth(:,t),'g','LineWidth',0.5)
%
%             plot(yasl(1,:,1),Wswath_depth_sm1(:,t),'r-')
%             plot(yasl(1,:,1),Wcaswath_depth_sm1(:,t),'g-')
%
%             plot(yasl(1,:,1),Wswath_depth_sm2(:,t),'r-','LineWidth',1.5)
%             plot(yasl(1,:,1),Wcaswath_depth_sm2(:,t),'g-','LineWidth',1.5)
%             axis([-110 80 0 5])
%             xlabel( 'y (km)' )
%             ylabel( 'Depth of W swath (km)' )
%
%             subplot(4,1,3)
%             hold on
%             %     plot(yasl(1,:,1), mean(Wswath_width(:,:,t),2,'omitnan'),'ob')
%             %     plot(yasl(1,:,1), mean(Wcaswath_width(:,:,t),2,'omitnan'),'og')
%             plot(yasl(1,:,1), mean(Wswath_width(:,:,t),2,'omitnan'),'b','LineWidth',0.5)
%             plot(yasl(1,:,1), mean(Wcaswath_width(:,:,t),2,'omitnan'),'g','LineWidth',0.5)
%
%             plot(yasl(1,:,1), Wswath_width_sm1(:,t),'b-')
%             plot(yasl(1,:,1), Wcaswath_width_sm1(:,t),'g-')
%
%             plot(yasl(1,:,1), Wswath_width_sm2(:,t),'b-','LineWidth',1.5)
%             plot(yasl(1,:,1), Wcaswath_width_sm2(:,t),'g-','LineWidth',1.5)
%
%             axis([-110 80 0 5])
%             xlabel( 'y (km)' )
%             ylabel( 'Mean width of W swath (km)' )
%
%             subplot(4,1,4)
%             hold on
%             %     plot(yasl(1,:,1), Wswath_maxmag(:,t),'ok')
%             %     plot(yasl(1,:,1), Wcaswath_maxmag(:,t),'og')
%             plot(yasl(1,:,1),Wswath_maxmag(:,t),'k','LineWidth',0.5)
%             plot(yasl(1,:,1),Wcaswath_maxmag(:,t),'g','LineWidth',0.5)
%
%             plot(yasl(1,:,1),Wswath_maxmag_sm1(:,t),'k-')
%             plot(yasl(1,:,1),Wcaswath_maxmag_sm1(:,t),'g-')
%
%             plot(yasl(1,:,1),Wswath_maxmag_sm2(:,t),'k-','LineWidth',1.5)
%             plot(yasl(1,:,1),Wcaswath_maxmag_sm2(:,t),'g-','LineWidth',1.5)
%             axis([-110 80 0 8])
%             xlabel( 'y (km)' )
%             ylabel( 'Max magnitude of W swath (km)' )
%
%             title([YYMMDD_all(t,:),' - ',HHMMSS_all(t,:),' (',num2str(t),')'])
%
%             pic = strcat(wrfdir,'AlongSDC_Wstr_width_depth_refl_t',num2str(t),'.png');
%             saveas(fig,pic)
%             close(fig)
%         end
%     end

terrfollow_westSDCbox_Uprofs = westSDCbox_Uprofs;
terrfollow_westSDCbox_Vprofs = westSDCbox_Vprofs;
terrfollow_eastSDCbox_Uprofs = eastSDCbox_Uprofs;
terrfollow_eastSDCbox_Vprofs = eastSDCbox_Vprofs;


disp('     ')
disp(' saving output into mat file: ' )
disp(matoutsave)
save(matoutsave,'wrfdir','runlab','Lon','Lat',...
    'Yvcsa','Xvcsa','Zvcsa','Tvcsa','A_Ht_asl_vcs','yasl','Terr','dx','HHMMSS_all','YYMMDD_all','newz_asl',...
    ...
    'FroudeEjt', 'FroudeWjt', 'east_N_mean_Fr', 'west_N_mean_Fr', 'east_upslope_mean_Fr', 'west_upslope_mean_Fr', 'HM_Fr', 'west_upslope_depth', 'east_upslope_depth',...
    ...
    'dzpeakU5_vent','dzpeakU8_vent','dzpeakU10_vent','dzpeakU15_vent','dzpeakU20_vent',...
    'Z5_vent','Z8_vent','Z10_vent','Z15_vent','Z20_vent',...
    ...
    'terrfollow_westSDCbox_Uprofs','terrfollow_westSDCbox_Vprofs','terrfollow_eastSDCbox_Uprofs','terrfollow_eastSDCbox_Vprofs','newz_asl',... 
    '-v7.3')


%swathbox:
 








% toc
% disp('    ')
% disp( 'DONEski with case: ' )
% disp(runlab)


%end %RUN loop



disp('        Fin. The end. The end has ended. What are you still doing here? Go home!           '         )


