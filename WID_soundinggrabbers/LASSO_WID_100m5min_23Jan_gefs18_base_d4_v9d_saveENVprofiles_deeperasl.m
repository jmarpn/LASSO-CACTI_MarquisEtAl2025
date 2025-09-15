

clear all


%{
 delete(gcp('nocreate'))
 pc = parcluster('local')
 parpool(pc, 128);
 spmd rank = labindex;
     fprintf(1,'Hello from %d\n',rank);
 end
 disp(' START')
 now
%}




% % home
%rootdir = '/Volumes/LaCie/LASSO_runs/cacti/lasso/les/';
%rootdir = '/Users/marq789/Downloads/';

% % nersc
rootdir = '/pscratch/sd/j/jmarquis/cacti/lasso/d4_100m_5min/'
RUNS = ['23Jan100m/20190123/gefs18/base/les/subset_d4'];
[ar br] = size( RUNS ) ;  clear br
R = 1;
wrfdir = strcat(rootdir,RUNS(R,:) ,'/f15min/')  ;

deepmetdir = strcat(rootdir,RUNS(R,:) ,'/deepermet/')  ;



%{
% wID, the home game:
wrfdir = '/Users/marq789/Downloads/tester/';
deepmetdir = '/Users/marq789/Downloads/tester/';
%}

% %home tester:
%wrfdir = ['/Users/marq789/Downloads/wenvsave/'];
 





WIDfilelist = ls(strcat(wrfdir,'/wID_ObjectOutput_NSj_100m5min_qcldice0.0001_dbz10_*mat'));
widlist = split(WIDfilelist);   widlist(end)=[] ;
[wa wb] = size(widlist); clear wb;

% %post-processed LASSO AGL cloud data files:
% cldlist = ls( horzcat(wrfdir,'/trimxy_corlasso.cldhaagl.*nc') );
% cldlist = split(cldlist);
% [sa sb] = size(cldlist); clear sb;
% 
% %post-processed LASSO lifted parcel metric files:
% parcellist = ls( horzcat(wrfdir,'/trimxy_corlasso.liftphagl*nc') );
% parcellist = split(parcellist);
% [sa sb] = size(parcellist); clear sb;

%post-processed LASSO met files:
metlist = ls( horzcat(deepmetdir,'/trimxy_asl_corlasso_methamsl_*nc') );
metlist = split(metlist);  metlist(end)=[];
[sa sb] = size(metlist); clear sb;


P_met = ncread(char(metlist(1)),'PRESSURE');
[am bm zm] = size(P_met);   
clear P_met


% 
% 
% %somewhat hardcoded lable for output
% wherethe2at = find(wrfdir == '2');
% wherethe4at = find(wrfdir == '4');
% runlab = wrfdir(wherethe2at(1):wherethe4at(end))   ;
% runlab(runlab == '/') = '_'  ;
% %matoutsave = strcat( wrfdir, 'wID_ObjectOutput_100m_qcld0.0001_dbz10_run_',runlab,'.mat' )  ;



%%% LASSo grid deltax,y (km)
dx = 0.1;

angles = [0 45 90 90+45 180 180+45 270 270+45] ;  % sampling of profile angles around circle in degrees

samplingKM = 2.0;  %km distance N, S, E, W from centroid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% START DOING STUFF NOW:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('   ')
disp(' processing data in : ')
disp( wrfdir )
disp('   ')



tback = 3;   %number of time indices to look backward drom wid t_ind for met data
% 3 = 15min,  6 = 30 min, 9 = 45 min, 12 = 60 min


for tt = 1 : wa
    
    %   tt = 5
    
    % preceeding updrafts outside of look-back window:
    if( (tt - tback) < 1)


        disp('   ')
        disp([ 'loading wIDs at time index ',num2str(tt) ] )

        widfile = char(widlist(tt));
        load(widfile,'wIDprofs_centlat','wIDprofs_centlon','wIDprofs_equivarea');

        [numwid nz]  = size(wIDprofs_centlat);  clear nz;

        %define file name stuff:
        segfile = split(widfile,"/");  segfile(1) = [];
        filestamp = segfile(16);  filestamp = split(filestamp,"_");
        DATE = filestamp(9);
        UTC = filestamp(8);
        ENS = filestamp(10);
        MP = filestamp(11);
        %DOM = filestamp(15);
        DOM = filestamp(14);

        WID_3Dtemper = zeros(numwid,length(angles)+1,zm);    WID_3Dtemper(:) = NaN;
        WID_3Dqvapor = zeros(numwid,length(angles)+1,zm);    WID_3Dqvapor(:) = NaN;
        WID_3Dpress  = zeros(numwid,length(angles)+1,zm);    WID_3Dpress(:) = NaN;
        WID_3Drh     = zeros(numwid,length(angles)+1,zm);    WID_3Drh(:) = NaN;
	    WID_3Du      = zeros(numwid,length(angles)+1,zm);    WID_3Du(:) = NaN;
        WID_3Dv      = zeros(numwid,length(angles)+1,zm);    WID_3Dv(:) = NaN;
        WID_3Dhtagl  = zeros(numwid,length(angles)+1,zm);    WID_3Dhtagl(:) = NaN;
        WID_3Dhtasl  = zeros(numwid,length(angles)+1,zm);    WID_3Dhtasl(:) = NaN;
        WID_3Dterr   = zeros(numwid,length(angles)+1);       WID_3Dterr(:) = NaN;



    else

        %%%%%%%%%%%%%%%%%%%%%%%%

        disp('   ')
        disp([ 'loading wIDs at time index ',num2str(tt) ] )

        widfile = char(widlist(tt));
        load(widfile,'wIDprofs_centlat','wIDprofs_centlon','wIDprofs_equivarea');

        [numwid nz]  = size(wIDprofs_centlat);  clear nz;

        %define file name stuff:
        disp(' segfile')
	segfile = split(widfile,"/");  segfile(1) = []

	disp(' filestamp ')
        filestamp = segfile(16)  
	filestamp = split(filestamp,"_")
        disp(' date ')
 	DATE = filestamp(9)
        
	disp(' utc ')
	UTC = filestamp(8)
        
	disp(' ensemble ');
	ENS = filestamp(10)
        
	disp(' mp ');
	MP = filestamp(11)
        
	disp(' DOM ');
	%DOM = filestamp(15);
        DOM = filestamp(14)

        %%%%%%%%%%%%%%%%%%%%%%%%

        %%%% loop tback here if you want multiples.

        disp('   ')
        disp('loading 3d met at time')
        disp( num2str(tt - tback) )

        metfile = char(metlist(tt - tback) )   % look back in time for met

        % metfile = char(metlist(end)) 

        Temp_met = ncread(metfile,'TEMPERATURE');       %K
        QVAP_met = ncread(metfile,'QVAPOR');
	    RH_met = ncread(metfile,'RH');
        P_met = ncread(metfile,'PRESSURE');    %mb
        U_met = ncread(metfile,'UA');
        V_met = ncread(metfile,'VA');
        Lon_met = ncread(metfile,'XLONG');
        Lat_met = ncread(metfile,'XLAT');
        Terr = ncread(metfile,'HGT');



        %%%%  [9Dec2024] some mods because now I'm reading in asl files instead of the agl
        % files I was using before: 

        Ht = ncread(metfile,'HAMSL');

        %calculate height 3 ASL field
        Ht_asl = P_met(:,:,:,1); Ht_asl(:) = NaN;
        for i = 1:am
            for j = 1:bm
                Ht_asl(i,j,:) = Ht ;
            end
        end

        %3d-ize height agl field for later conveneience:
        Ht_agl = P_met(:,:,:,1); Ht_agl(:) = NaN;
        for i = 1:am
            for j = 1:bm
                Ht_agl(i,j,:) =  Ht_asl(i,j,:) - Terr(i,j);
            end
        end
        Ht_agl( Ht_agl < 0.0 ) = NaN;



        %{
        %%%% when you used to read in AGL-frame data:
        %calculate height ASL field
        Ht_asl= P_met(:,:,:,1); Ht_asl(:) = NaN;
        for i = 1:am
            for j = 1:bm
                Ht_asl(i,j,:) = Ht  + Terr(i,j);
            end
        end

        %3d-ize height agl field for later conveneience:
        Ht_agl = P_met(:,:,:,1); Ht_agl(:) = NaN;
        for i = 1:am
            for j = 1:bm
                Ht_agl(i,j,:) = Ht;
            end
        end
        %}

        % set up a surface value
        P_sfc = ncread(metfile,'PSFC');    %pa
        %P_sfc = P_sfc/100;  %mb
        U_10m = ncread(metfile,'U10');    %m/s
        V_10m = ncread(metfile,'V10');
        T_2m = ncread(metfile,'T2');      %K
        Qvap_2m = ncread(metfile,'Q2');   %kg/kg     

        tdry_2m = T_2m - 273.15;  % degC
        esat_2m =  6.112*exp( (17.67*tdry_2m)./(243.5+tdry_2m) ); % mb
        e_2m    = ( ( (Qvap_2m .* P_sfc)/0.622 ) ./ ( 1 - Qvap_2m./0.622 ) )/100 ;   %mb
        RH_2m   = 100*( e_2m./esat_2m );

        disp('resetting bottom height')

        tic
        %add the surface values to the grid:
        for i = 1 : am
            for j = 1 : bm 
           
                %  i = 459; j = 1569;


                k_firstdirt = find(  isnan( Ht_agl(i,j,:) )  )  ;
                k_firstdirt = k_firstdirt(end)  ;   %the k index that is the first one under the surface

                %found a few points where the first hagl point is below 10m
                if(  Ht_agl(i,j, k_firstdirt+1)  > 10 )  %when first level above dirt is > 10m

                    %rebrand this first dirt as the 10m height:

                    Ht_agl(i,j,k_firstdirt) = 10.0 ;
                    Ht_asl(i,j,k_firstdirt) = Terr(i,j) + 10.0 ;   %not sure if I use this henceforth?

                    U_met(i,j,k_firstdirt)  = U_10m(i,j) ;
                    V_met(i,j,k_firstdirt)  = V_10m(i,j) ;

                    %interpolate to get a 10m T & Q & P:
                    P_2mtemporary   = [  P_sfc(i,j)/100, P_met(i,j,k_firstdirt+1), P_met(i,j,k_firstdirt+2)    ] ;  %mb
                    T_2mtemporary   = [  T_2m(i,j), Temp_met(i,j,k_firstdirt+1), Temp_met(i,j,k_firstdirt+2)    ] ;
                    Q_2mtemporary   = [  Qvap_2m(i,j), QVAP_met(i,j,k_firstdirt+1), QVAP_met(i,j,k_firstdirt+2)  ] ;
                    AGL_2mtemporary = [  2.0, Ht_agl(i,j,k_firstdirt+1), Ht_agl(i,j,k_firstdirt+2) ];

                    T_10m = interp1( AGL_2mtemporary ,T_2mtemporary, 10. ) ;
                    Q_10m = interp1( AGL_2mtemporary ,Q_2mtemporary, 10. ) ;
                    P_10m = interp1( AGL_2mtemporary ,P_2mtemporary, 10. ) ;

                    %insert that T10, Q10 into full T,Q array
                    Temp_met(i,j,k_firstdirt)  = T_10m ;
                    QVAP_met(i,j,k_firstdirt)  = Q_10m ;
                    P_met(i,j,k_firstdirt)     = P_10m ;

                elseif( Ht_agl(i,j, k_firstdirt+1)  <= 10  )  % hagl is between 0-10m, so some sounding values missing. let's replace them

                    %rebrand this sub-10m first above dirt as the 10m height:

                    Ht_agl(i,j,k_firstdirt+1) = 10.0 ;
                    Ht_asl(i,j,k_firstdirt+1) = Terr(i,j) + 10.0 ;   %not sure if I use this henceforth?

                    U_met(i,j,k_firstdirt+1)  = U_10m(i,j) ;
                    V_met(i,j,k_firstdirt+1)  = V_10m(i,j) ;

                    %interpolate to get a 10m T & Q & P:
                    P_2mtemporary   = [  P_sfc(i,j)/100, P_met(i,j,k_firstdirt+2), P_met(i,j,k_firstdirt+3)    ] ;  %mb
                    T_2mtemporary   = [  T_2m(i,j), Temp_met(i,j,k_firstdirt+2), Temp_met(i,j,k_firstdirt+3)    ] ;
                    Q_2mtemporary   = [  Qvap_2m(i,j), QVAP_met(i,j,k_firstdirt+2), QVAP_met(i,j,k_firstdirt+3)  ] ;
                    AGL_2mtemporary = [  2.0, Ht_agl(i,j,k_firstdirt+2), Ht_agl(i,j,k_firstdirt+3) ];

                    T_10m = interp1( AGL_2mtemporary ,T_2mtemporary, 10. ) ;
                    Q_10m = interp1( AGL_2mtemporary ,Q_2mtemporary, 10. ) ;
                    P_10m = interp1( AGL_2mtemporary ,P_2mtemporary, 10. ) ;

                    %insert that T10, Q10 into full T,Q array
                    Temp_met(i,j,k_firstdirt+1)  = T_10m ;
                    QVAP_met(i,j,k_firstdirt+1)  = Q_10m ;
                    P_met(i,j,k_firstdirt+1)     = P_10m ;


                end




            end
        end
        toc

        %%%%% done with agl-asl input mods

        %%%%%% seed 3D environmental arrays

        WID_3Dtemper = zeros(numwid,length(angles)+1,zm);    WID_3Dtemper(:) = NaN;
        WID_3Dqvapor = zeros(numwid,length(angles)+1,zm);    WID_3Dqvapor(:) = NaN;
        WID_3Drh     = zeros(numwid,length(angles)+1,zm);    WID_3Drh(:) = NaN;	
        WID_3Dpress  = zeros(numwid,length(angles)+1,zm);    WID_3Dpress(:) = NaN;
        WID_3Du      = zeros(numwid,length(angles)+1,zm);    WID_3Du(:) = NaN;
        WID_3Dv      = zeros(numwid,length(angles)+1,zm);    WID_3Dv(:) = NaN;
        WID_3Dhtagl  = zeros(numwid,length(angles)+1,zm);    WID_3Dhtagl(:) = NaN;
        WID_3Dhtasl  = zeros(numwid,length(angles)+1,zm);    WID_3Dhtasl(:) = NaN;
        WID_3Dterr   = zeros(numwid,length(angles)+1);       WID_3Dterr(:) = NaN;

        disp('grabbing wid envs now')

        for n = 1:numwid
             
            n

            % n = 100

            latcent_wid = mean( wIDprofs_centlat(n,:),'omitnan' );
            loncent_wid = mean( wIDprofs_centlon(n,:),'omitnan' );
            circradM_wid =  ( max( wIDprofs_equivarea(n,:), [], 2, 'omitnan' ) ./ 3.14159) .^ 0.5 ; 
            circradKM_wid =  circradM_wid ./ 1000.0 ;  

            circradKM_wid =  circradKM_wid * 3 ;


	    %disp('AA')

            latdif =  abs(Lat_met - latcent_wid);   londif  =  abs(Lon_met - loncent_wid);
            totdif = latdif + londif;
            [loncent_met latcent_met] = find( totdif == min(min(totdif))  );  %nearest met grid indices of wid centroid


            % figure; contourf(Terr,30); hold on; contour(totdif ,40, 'r'); plot(latcent_met,loncent_met,'xr'); plot(latcent_met,loncent_met+50,'ko');
            % figure; contourf(Lon_met,Lat_met,Terr,30)
            %         figure; contourf(Terr,30); hold on; contour(totdif ,40, 'r'); plot(latcent_met,loncent_met,'xr'); plot(latcent_met,loncent_met+50,'ko');

            %         correct resulting index referncing format:
            %         Lat_met(loncent_met,latcent_met)
            %         Lon_met(loncent_met,latcent_met)

            % %sometimes tehre are multiple nearest neighbors
            %loncent_met = loncent_met(1);
            %latcent_met = latcent_met(1);

            %disp('BB')

            
            %if( isnan(latcent_wid)==0 &  isempty(loncent_met)==0  &  isempty(latcent_met)==0  &   ...
            %        latcent_met + floor(circradKM_wid/dx) + samplingKM/dx < bm   &   latcent_met - floor(circradKM_wid/dx) - samplingKM/dx > 0  & ...
            %        loncent_met + floor(circradKM_wid/dx) + samplingKM/dx < am   &   loncent_met - floor(circradKM_wid/dx) - samplingKM/dx > 0  )
 
	    if( isnan(latcent_wid)==0 &  isempty(loncent_met)==0  &  isempty(latcent_met)==0  &   ...
                 latcent_met + floor(circradKM_wid/dx) < bm   &   latcent_met - floor(circradKM_wid/dx) > 0 & ...
                 loncent_met + floor(circradKM_wid/dx) < am   &   loncent_met - floor(circradKM_wid/dx) > 0 )

                 %sometimes there are multiple nearest neighbors
                 loncent_met = loncent_met(1);
                 latcent_met = latcent_met(1);

                %disp('CC')

                %wR = floor(circradKM_wid/dx) + samplingKM/dx ;
                wR = floor(circradKM_wid/dx)  ;
                
	        %disp('DD')	
                %index distance from centroid going around circle

%                 i_r = latcent_met + wR *  cosd( 180 + 270 - angles ) ;
%                 j_r = loncent_met + wR *  sind( 180 + 270 - angles ) ;
                
                %Future JIm: after the usual tedious matlab-i-j and cos-sin
                %annoying diagnostics, I'm pretty sure the code here and
                %below is all consistient in N -> NE -> E -> SE -> S - ...
                i_r = latcent_met + wR *  cosd(  angles ) ;
                j_r = loncent_met + wR *  sind(  angles ) ;

                %disp('EE')

                %discretize the indices to nearest:
                i_r = round(i_r)  ;
                j_r = round(j_r)  ;

                %{
                %%% diagnostic plots
                figure;
                contourf( Lon_met, Lat_met, T_2m, 20 )
                hold on
                plot(Lon_met(loncent_met, latcent_met), Lat_met(loncent_met, latcent_met) , 'o', 'MarkerFaceColor', 'r', 'MarkerSize', 10)
                plot(Lon_met(j_r(1), i_r(1)), Lat_met(j_r(1), i_r(1)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 15)
                plot(Lon_met(j_r(2), i_r(2)), Lat_met(j_r(2), i_r(2)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10)
                plot(Lon_met(j_r(3), i_r(3)), Lat_met(j_r(3), i_r(3)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
                plot(Lon_met(j_r(4), i_r(4)), Lat_met(j_r(4), i_r(4)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 6)
                plot(Lon_met(j_r(5), i_r(5)), Lat_met(j_r(5), i_r(5)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
                plot(Lon_met(j_r(6), i_r(6)), Lat_met(j_r(6), i_r(6)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
                plot(Lon_met(j_r(7), i_r(7)), Lat_met(j_r(7), i_r(7)) , 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
                plot(Lon_met(j_r(8), i_r(8)), Lat_met(j_r(8), i_r(8)) , 'd', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
                axis equal
                axis([-65 -64.88 -32.65 -32.55])
                %}

                % at w centroid location  (profile index 1):
                WID_3Dtemper(n,1,:) =  Temp_met(loncent_met,latcent_met,:);
                WID_3Dqvapor(n,1,:) =  QVAP_met(loncent_met,latcent_met,:);
                WID_3Drh(n,1,:)     =  RH_met(loncent_met,latcent_met,:);		
                WID_3Dpress(n,1,:)  =  P_met(loncent_met,latcent_met,:);
                WID_3Du(n,1,:)      =  U_met(loncent_met,latcent_met,:);
                WID_3Dv(n,1,:)      =  V_met(loncent_met,latcent_met,:);
                WID_3Dhtagl(n,1,:)  =  Ht_agl(loncent_met,latcent_met,:);
                WID_3Dhtasl(n,1,:)  =  Ht_asl(loncent_met,latcent_met,:);
                WID_3Dterr(n,1)     =  Terr(loncent_met,latcent_met);


		%disp('FF')

                % now go around the circle. it SHOULD BE (acrroding to the above diagnostic plots)
                % go from N(1) -> NE(2) -> E(3) -> SE -> S -> SW -> W -> NW(8) 
                for an = 1 : length(angles)

                    WID_3Dtemper(n,an+1,:) =  Temp_met(j_r(an),i_r(an),:);
                    WID_3Dqvapor(n,an+1,:) =  QVAP_met(j_r(an),i_r(an),:);
                    WID_3Drh(n,an+1,:)     =  RH_met(j_r(an),i_r(an),:);
                    WID_3Dpress(n,an+1,:)  =  P_met(j_r(an),i_r(an),:);
                    WID_3Du(n,an+1,:)      =  U_met(j_r(an),i_r(an),:);
                    WID_3Dv(n,an+1,:)      =  V_met(j_r(an),i_r(an),:);
                    WID_3Dhtagl(n,an+1,:)  =  Ht_agl(j_r(an),i_r(an),:);
                    WID_3Dhtasl(n,an+1,:)  =  Ht_asl(j_r(an),i_r(an),:);
                    WID_3Dterr(n,an+1)     =  Terr(j_r(an),i_r(an));

                end

                %disp('GG')

                %{
%%% more diagnostics plots to make sure that i's and j's are where you
%%% intend them to be. Goes along with the diag plot above 
% %center point:
%tem = WID_3Dtemper(n,1,:) ; 
%sfc = find( isnan(tem) ==0) ;
%sfc = sfc(1) ; 
%tsfc = tem(sfc)
%
% % circle point 1 (should be north?)
%tem = WID_3Dtemper(n,2,:) ; 
%sfc = find( isnan(tem) ==0) ;
%sfc = sfc(1) ; 
%tsfc = tem(sfc)
%
% circle point 3 (should be east )
%tem = WID_3Dtemper(n,4,:) ; 
%sfc = find( isnan(tem) ==0) ;
%sfc = sfc(1) ; 
%tsfc = tem(sfc)
%
% % circle point 5 (should be due south)
%tem = WID_3Dtemper(n,6,:) ; 
%sfc = find( isnan(tem) ==0) ;
%sfc = sfc(1) ; 
%tsfc = tem(sfc)
%
% % circle point 7 (should be due west)
%tem = WID_3Dtemper(n,8,:) ; 
%sfc = find( isnan(tem) ==0) ;
%sfc = sfc(1) ; 
%tsfc = tem(sfc)
                %}


%                 % at w centroid location  (profile index 1):
%                 WID_3Dtemper(n,1,:) =  Temp_met(loncent_met,latcent_met,:);
%                 WID_3Dqvapor(n,1,:) =  QVAP_met(loncent_met,latcent_met,:);
%                 WID_3Drh(n,1,:)     =  RH_met(loncent_met,latcent_met,:);		
%                 WID_3Dpress(n,1,:)  =  P_met(loncent_met,latcent_met,:);
%                 WID_3Du(n,1,:)      =  U_met(loncent_met,latcent_met,:);
%                 WID_3Dv(n,1,:)      =  V_met(loncent_met,latcent_met,:);
%                 WID_3Dhtagl(n,1,:)  =  Ht_agl(loncent_met,latcent_met,:);
%                 WID_3Dhtasl(n,1,:)  =  Ht_asl(loncent_met,latcent_met,:);
%                 WID_3Dterr(n,1)     =  Terr(loncent_met,latcent_met);
% 
%                 %now go around the circle from N -> NE -> E -> SE -> S -> SW -> W -> NW -> N 
% 
% 
%                 % north of w centroid location  (profile index 2):
%                 WID_3Dtemper(n,2,:) =  Temp_met(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Dqvapor(n,2,:) =  QVAP_met(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Drh(n,2,:)     =  RH_met(loncent_met,latcent_met + samplingKM/dx,:);		
%                 WID_3Dpress(n,2,:)  =  P_met(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Du(n,2,:)      =  U_met(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Dv(n,2,:)      =  V_met(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Dhtagl(n,2,:)  =  Ht_agl(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Dhtasl(n,2,:)  =  Ht_asl(loncent_met,latcent_met + samplingKM/dx,:);
%                 WID_3Dterr(n,2)     =  Terr(loncent_met,latcent_met + samplingKM/dx);
% 
%                 % south of w centroid location  (profile index 3):
%                 WID_3Dtemper(n,3,:) =  Temp_met(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Dqvapor(n,3,:) =  QVAP_met(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Drh(n,3,:)     =  RH_met(loncent_met,latcent_met - samplingKM/dx,:);		
%                 WID_3Dpress(n,3,:)  =  P_met(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Du(n,3,:)      =  U_met(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Dv(n,3,:)      =  V_met(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Dhtagl(n,3,:)  =  Ht_agl(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Dhtasl(n,3,:)  =  Ht_asl(loncent_met,latcent_met - samplingKM/dx,:);
%                 WID_3Dterr(n,3)     =  Terr(loncent_met,latcent_met - samplingKM/dx);
% 
%                 % east of w centroid location  (profile index 4):
%                 WID_3Dtemper(n,4,:) =  Temp_met(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Dqvapor(n,4,:) =  QVAP_met(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Drh(n,4,:)     =  RH_met(loncent_met + samplingKM/dx,latcent_met,:);                
% 		        WID_3Dpress(n,4,:)  =  P_met(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Du(n,4,:)      =  U_met(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Dv(n,4,:)      =  V_met(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Dhtagl(n,4,:)  =  Ht_agl(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Dhtasl(n,4,:)  =  Ht_asl(loncent_met + samplingKM/dx,latcent_met,:);
%                 WID_3Dterr(n,4)     =  Terr(loncent_met + samplingKM/dx,latcent_met);
% 
%                 % west of w centroid location  (profile index 4):
%                 WID_3Dtemper(n,5,:) =  Temp_met(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Dqvapor(n,5,:) =  QVAP_met(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Drh(n,5,:)     =  RH_met(loncent_met - samplingKM/dx,latcent_met,:);
%  		        WID_3Dpress(n,5,:)  =  P_met(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Du(n,5,:)      =  U_met(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Dv(n,5,:)      =  V_met(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Dhtagl(n,5,:)  =  Ht_agl(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Dhtasl(n,5,:)  =  Ht_asl(loncent_met - samplingKM/dx,latcent_met,:);
%                 WID_3Dterr(n,5)     =  Terr(loncent_met - samplingKM/dx,latcent_met);

            end  %position checks

        end   %wid loop

    end % tlookback horizon

    %     DATE = filestamp(10)
    %     UTC = filestamp(8)
    %     ENS = filestamp(11)
    %     MP = filestamp(12)
    %     DOM = filestamp(15)

    matout = strcat( wrfdir, 'wID_3Denv_3R_',string(DATE),'_',string(UTC),'_tindlookback',num2str(tback),'_',string(ENS),'_',MP,'_',DOM  )

    save( matout,'WID_3Dtemper','WID_3Dqvapor','WID_3Drh','WID_3Dpress','WID_3Du','WID_3Dv','WID_3Dhtagl','WID_3Dhtasl','WID_3Dterr',...
        'DATE','UTC','ENS','MP','DOM','tback','samplingKM','WIDfilelist','metlist')



end  %time loop






% clear all
% load('/Users/marq789/Downloads/wenvsave/wID_3Denv_20181129_151000_tindlookback2_gefs09_base_d4.mat')










