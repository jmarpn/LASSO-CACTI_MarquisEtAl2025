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


% v9: separated into 1) swath [ID, log thermo, swathbox calcs], 2) wID

clear all



 delete(gcp('nocreate'))
 pc = parcluster('local')
 parpool(pc, 128);
 spmd rank = labindex;
     fprintf(1,'Hello from %d\n',rank);
 end
 disp(' START')
 now





% % home
%rootdir = '/Volumes/LaCie/LASSO_runs/cacti/lasso/les/';
%rootdir = '/Users/marq789/Downloads/';

% % nersc
rootdir = '/pscratch/sd/j/jmarquis/cacti/lasso/d4_100m_5min/'

%nersc:
RUNS = ['/25Jan100m/20190125/eda07/base/les/subset_d4'];


%home tester:
  %rootdir = ['/Volumes/LaCie/LASSO_runs/cacti/lasso/les/'];
  %RUNS = ['20190129/gefs11/base/d4/'];

%home:
%RUNS = ['20190129/gefs11/base/d4/'];  %  "best" of the day?


[ar br] = size( RUNS ) ;  clear br


%loop through ensemble member runs:
%for R = 1 : ar

%    clearvars -except rootdir RUNS ar R

R = 1

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
%matoutsave = strcat( wrfdir, 'wID_ObjectOutput_100m_qcld0.0001_dbz10_run_',runlab,'.mat' )  ;

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


%post-processed LASSO AGL cloud data files:
cldlist = ls( horzcat(wrfdir,'/trimxy_corlasso.cldhaagl.*nc') );
cldlist = split(cldlist);
[sa sb] = size(cldlist); clear sb;

%post-processed LASSO lifted parcel metric files:
parcellist = ls( horzcat(wrfdir,'/trimxy_corlasso.liftphagl*nc') );
parcellist = split(parcellist);
[sa sb] = size(parcellist); clear sb;

%post-processed LASSO met files:
metlist = ls( horzcat(wrfdir,'/trimxy_corlasso.methagl*nc') );
metlist = split(metlist);  metlist(end)=[];
[sa sb] = size(metlist); clear sb;

%%%%%%%%%%%%% diagnostic tool
% sa = 3;   %tester: reduction of the full sample set
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I looked at the terrain field and chose a select subdomain, so trim away some domain fat early:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  % AGL SUBDOMAIN TO LOAD
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


% dum = ncread(metf,'UA') ;
% [aa bb cc] = size(dum) ;  clear dum
% dummy = zeros(di,dj,cc,sa-1) ;
% dummy(:) = NaN;

% Temp_all = dummy;
%     THETA_all = dummy;
%     THETAE_all = dummy;
%     QVAP_all = dummy;
%     RH_all = dummy;
%     P_all = dummy;

%     U_all = dummy;
%     V_all = dummy;
% W_all = dummy;

%     DBZ_all = dummy;
% Qcld_all = dummy;

%clear dummy

%MULFC_all = zeros(di,dj,sa-1) ;      MULFC_all(:) = NaN;
%     MUCIN_all = zeros(di,dj,sa-1) ;      MUCIN_all(:) = NaN;
%     MUCAPE_all = zeros(di,dj,sa-1) ;     MUCAPE_all(:) = NaN;



disp('   ')
disp( ' loading post-processed AGL LASSO files for: '    )
disp(runlab)
disp('   ')

for tt = 1 : sa

    disp('   ')
    disp('loading 3d vars for time')
    disp( num2str(tt) )

    metfile = char(metlist(tt));   % ncdisp(metfile)
    cldfile = char(cldlist(tt));
    parcelfile = char(parcellist(tt));

    start3d = [i1 j1 1 1];
    count3d = [di dj Inf 1];
    start2d = [i1 j1 1];
    count2d = [di dj 1];

    W_all = ncread(metfile,'WA',start3d,count3d);
    Temp_all = ncread(metfile,'TEMPERATURE',start3d,count3d);       %K
    QVAP_all = ncread(metfile,'QVAPOR',start3d,count3d);
    P_all = ncread(metfile,'PRESSURE',start3d,count3d);    %mb
    DBZ_all = ncread(cldfile, 'REFL_10CM',start3d,count3d);
    Qcld_all = ncread(cldfile, 'QCLOUD',start3d,count3d);
    MULFC_all = ncread(parcelfile, 'MULFC',start2d,count2d);      % according to bill, this is in ht ASL frame
    Qice_all = ncread(cldfile, 'QICE',start3d,count3d);

    % semi-hardcoded time/date labels for LASSO data:

    metfile = char(metlist(tt)) ;
    YYMMDD_all = metfile(end-17:end-10)
    HHMMSS_all = metfile(end-8:end-3)

    matoutsave = strcat( wrfdir, 'wID_ObjectOutput_NSj_100m5min_qcldice0.0001_dbz10_tt',num2str(tt,'%03.f'),'_',HHMMSS_all,'_',runlab,'.mat' )  ;


    % use this later for low-level updraft ID'ing
    MULFC_agl_all = MULFC_all;   MULFC_agl_all(:) = NaN;   %  MULFC in AGL frame
    MULFC_agl_all = MULFC_all - Terr ;
    clear MULFC_all

    % move on to calculate desireable variables
    CLD_THRESH = 0.0001; DBZ_THRESH = 10;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % define cloudy and clear-air updrafts:
    Wcld_all = W_all; Wcld_all(:) = NaN;
    [xw yw kw]  = size(W_all);
    %[cloudyi cloudyj cloudyk] = find(Qcld_all > CLD_THRESH | DBZ_all > DBZ_THRESH);
    for i = 1:xw
        for j = 1:yw
            for k = 1:kw
                if( ( Qcld_all(i,j,k) + Qice_all(i,j,k) ) > CLD_THRESH | DBZ_all(i,j,k) > DBZ_THRESH)
		    Wcld_all(i,j,k) = W_all(i,j,k);
                end
            end
        end
    end
    
    clear Qcld_all DBZ_all Qice_all
    [Xs Ys Zs] = size(W_all);

    %     %%%%%%%%%%%
    %     %%%%%  vert mass flux calculations:
    %     %%%%%%%%%%%

    P_all = (P_all*100);            %pa
    tmp = (QVAP_all ./ 0.622);
    clear QVAP_all
    
    vappres = (tmp .* P_all) ./ ( 1 + tmp);       %pa
    rhom = (    P_all ./ (287 .* Temp_all)  ) .* (   1 -  (vappres./P_all).* ( 1 - 0.622 )    ) ;
    clear P_all tmp vappres Temp_all

    %     %VMca = Wca_all .* rhom ;
    VM = W_all .* rhom ;
    clear rhom %W_all


    %calculate height ASL field
    Ht_asl= VM(:,:,:,1); Ht_asl(:) = NaN;
    for i = 1:Xs
        for j = 1:Ys
            Ht_asl(i,j,:) = Ht  + Terr(i,j);
        end
    end


    %3d-ize height agl field for later conveneience:
    Ht_agl = VM(:,:,:,1); Ht_agl(:) = NaN;
    for i = 1:Xs
        for j = 1:Ys
            Ht_agl(i,j,:) = Ht;
        end
    end

    %%%%%% define the x,y arrays:
    X = Lon; X(:) = 0; Y = X ;
    xx = dx*[1:Xs]; yy = dx*[1:Ys] ;

    %3D-fy horiz grid for later conveneience:
    X = VM(:,:,:,1); X(:) = NaN; Y = X;
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


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%           wID:  CLOUDY UPDRAFT ID'ing
    %%%%%
    %%%%% catalog 3d traits of updraft masks from the ground to height
    %%%%% (index) = kagltop.
    %%%%%
    %%%%%  see:
    %%%%%  code transplanted from ...500m_v7e.m code
    %%%%%
    %%%%%  this part will connect ID cloudy updraft points in 3D of prescribed tresholds, locate their centroids,
    %%%%%    and compose 2D (vertical) profiles of area, magnitude, massflux. Then I do some cleanup on those fields
    %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %user prescribed input:
    wlowthresh = 2.5;   %weakest updraft to ID
    kagltop = 38;  %first guess of relevant ht index correcponding to:  29 = 4kmAGL;   33 = 5kmAGL;   37 = 9kmAGL
    thresh_updraft_minpts = 3;   % min number of updraft grid points on a horiz surface to calc an area
    maxnum_guess  = 2000;   %take a guess at max number of updrafts per individual LASSO analysis time (not full updraft lifetime or full LASSO period)
    kmaboveLFC = 10000;  %heights above & below LFC (m) in which to ultimately keep updraft info
    kmbelowLFC = 50;
    updraftLFCLAYprofs_zagl = Ht(1:kagltop);  %agl height profile being considered



    %take a first stab at 3d updraft masking at low-levels regardless of the LFC height (will be a first guess for the LFC height-based updraft IDing.
    disp(' ')
    disp(' IDing low-level cloudy updraft objects below kagltop (ignoring LFC constraint at the moment)')
    tic


    mulfcagl = MULFC_agl_all;    clear MULFC_agl_all

    %start by filtering 3d w field on low w threshold
    % only considering the lowest prescribed layer, from k = [1:kagltop]
    W3d = Wcld_all(:,:,1:kagltop);  
    W3d(W3d < wlowthresh) = 0;
    W3d(isnan(W3d)) = 0;
    Wcomp = max(W3d,[],3);



    %checkpoint
    if(  isempty( find(Wcomp ~= 0) )  )  %if all updrafts are filtered out, then skip this time in t-loop;  0 = keep going
        continue  %skip time
    end


    % now finding local w-maxima in 3D field; practially, doing it in a vertically composited 2D area
    % note, this will sometimes ID a few maxima within the exact same 3D updraft shell, so there will
    % be repeats at this point (culled down to 1 for each shell out later)
    wmaxes = imregionalmax(Wcomp,8);
    [wimaxs wjmaxs] = find(wmaxes == 1);


    %checkpoint
    if(  isempty(wimaxs)  )  %if all w-max points filtered out, then skip this time in loop;  0 = keep going
        continue    %skip time
    end


    % now find the height of each i,j max:
    wkmaxs = wimaxs; wkmaxs(:) = NaN;
    for m = 1:length(wimaxs)
        % m = find(imaxs == 104 & jmaxs == 190) ;   %tester
        wprof = W3d(wimaxs(m),wjmaxs(m),:); wprof = permute(wprof,[3 1 2]);
        wkmaxs(m) = find( wprof == max(wprof) );
    end



    %logic-ify the w field to feed through the 3d-pt connection algorithm:
    W3dbin = W3d;
    W3dbin(W3d(:,:,:) >  wlowthresh) = 1;
    W3dbin(W3d(:,:,:) <= wlowthresh) = 0;
    W3dbin = logical(W3dbin);
    [aa bb cc]= size(W3dbin);


    %seed Final output vars
    FupdraftLFCLAYprofs_equivarea  = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_equivarea(:) = NaN;
    FupdraftLFCLAYprofs_numij      = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_numij(:) = NaN;  %numbe of horiz grid pts per height
    FupdraftLFCLAYprofs_maxmag     = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_maxmag(:) = NaN;
    FupdraftLFCLAYprofs_massflux   = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_massflux(:) = NaN;
    FupdraftLFCLAYprofs_centlat    = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_centlat(:) = NaN;
    FupdraftLFCLAYprofs_centlon    = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_centlon(:) = NaN;
    FupdraftLFCLAYprofs_centHTagl  = zeros( maxnum_guess, kagltop);    FupdraftLFCLAYprofs_centHTagl(:) = NaN;
    FupdraftLFCLAYprofs_kbotlayer  = zeros( maxnum_guess,1);             FupdraftLFCLAYprofs_kbotlayer(:) = NaN;
    FupdraftLFCLAYprofs_ktoplayer  = zeros( maxnum_guess,1);             FupdraftLFCLAYprofs_ktoplayer(:) = NaN;
    FupdraftLFCLAYprofs_klfcagl    = zeros( maxnum_guess,1);             FupdraftLFCLAYprofs_klfcagl(:) = NaN;
%     FupdraftLFCLAYprofs_htasl_wtop = zeros( maxnum_guess,1);             FupdraftLFCLAYprofs_htasl_wtop(:) = NaN;  %top and bottom of each updraft obj in mASL
%     FupdraftLFCLAYprofs_htasl_wbot = zeros( maxnum_guess,1);             FupdraftLFCLAYprofs_htasl_wbot(:) = NaN;
    %fullest analysis time wind per day (over which zhe tracks cells)
    %%%%1marHOV
    celltrack_hr = [12:5/60:24];
    %seed Final output vars of wID stats along SDC (y,t) space - for Hovmueller diagrams:
    FupdraftLFCLAYprofs_northj  = zeros( maxnum_guess,kagltop);     FupdraftLFCLAYprofs_northj(:) = NaN;  %will be the northern and southern j extents of mask at each height
    FupdraftLFCLAYprofs_southj  = zeros( maxnum_guess,kagltop);     FupdraftLFCLAYprofs_southj(:) = NaN;



    % %for diagnostics
    % PREFILTupdraftLFCLAYprofs_centlat    = zeros( maxnum_guess, kagltop);    PREFILTupdraftLFCLAYprofs_centlat(:) = NaN;
    % PREFILTupdraftLFCLAYprofs_centlon    = zeros( maxnum_guess, kagltop);    PREFILTupdraftLFCLAYprofs_centlon(:) = NaN;
    %PREFILTupdraftLFCLAYprofs_centHTagl  = zeros( maxnum_guess, kagltop);    PREFILTupdraftLFCLAYprofs_centHTagl(:) = NaN;

    %%%%
    %%%% Now catalog all the updraft profile characteristics in each 3D updraft mask:
    %%%%

    disp(['num wmaxes at t',num2str(tt),' is  ',num2str(numel(wimaxs))])

    %parfor m = 1: maxnum_guess  %


    disp(['starting parallel wID'])
    tic
    parfor m = 1: numel(wimaxs)  % maxnum_guess  %

        %disp( ['start  m',num2str(m)] )

        updraftLFCLAYprofs_equivarea  = zeros(kagltop,1);    updraftLFCLAYprofs_equivarea(:) = NaN;
        updraftLFCLAYprofs_numij      = zeros(kagltop,1);    updraftLFCLAYprofs_numij(:) = NaN;  %numbe of horiz grid pts per height
        updraftLFCLAYprofs_maxmag     = zeros(kagltop,1);    updraftLFCLAYprofs_maxmag(:) = NaN;
        updraftLFCLAYprofs_massflux   = zeros(kagltop,1);    updraftLFCLAYprofs_massflux(:) = NaN;
        updraftLFCLAYprofs_centlat    = zeros(kagltop,1);    updraftLFCLAYprofs_centlat(:) = NaN;
        updraftLFCLAYprofs_centlon    = zeros(kagltop,1);    updraftLFCLAYprofs_centlon(:) = NaN;
        updraftLFCLAYprofs_centHTagl  = zeros(kagltop,1);    updraftLFCLAYprofs_centHTagl(:) = NaN;
        updraftLFCLAYprofs_kbotlayer  = NaN;              
        updraftLFCLAYprofs_ktoplayer  = NaN;           
        updraftLFCLAYprofs_klfcagl    = NaN;       
%         updraftLFCLAYprofs_htasl_wtop = NaN;            
%         updraftLFCLAYprofs_htasl_wbot = NaN;          
        %%%%%%%%1marHOV
        %record the north/south extent of the updraft's mask for Hovmuller purposes 
        updraftLFCLAYprofs_northj  = zeros(kagltop,1);      updraftLFCLAYprofs_northj(:) = NaN;   % N extent as f(z)
        updraftLFCLAYprofs_southj  = zeros(kagltop,1);      updraftLFCLAYprofs_southj(:) = NaN;   % S extent as f(z)


        if( m <= numel(wimaxs) )

            % current 3d binary updraft object
            Wobj3d = bwselect3(W3dbin,wjmaxs(m),wimaxs(m),wkmaxs(m),26) ;

            wob = Wobj3d; %wob(wob==0) = NaN;
            wob_zasl = Ht_asl .* Wobj3d;  wob_zasl(wob_zasl==0) = NaN;

            %check if there is any ht with 3 adjoining updraft points
            %(you found a few with a single point at a height(s))
            wht_robust_flag = 0;
            for k = 1:kagltop
                if( length( find( Wobj3d(:,:,k)==1 ) ) > 2 )
                    wht_robust_flag = 1;
                end
            end

            if( wht_robust_flag == 1 ) % there's at least one ht with an updraft with at least 3 adjoining pts.

                mzi = [];  mzj = [];  mzr = [];
                for k = 1:kagltop
                    % k = 37

                    ik = [];
                    jk = [];
                    [ik jk] = find(  Wobj3d(:,:,k) == 1  );  %all i,j points in the updraft mask at this time & height
                    numpts = length(ik)  ;

                    %calculate the vertical profiles of updraft metrics. Filter the small ones out:
                    if( isempty(ik) == 0   )

                        %%%%%% 1marHOV
                        %%% log the N-S shadow of wID at each height
                        updraftLFCLAYprofs_northj(k) =  max(jk) ;   % N extent as f(z)
                        updraftLFCLAYprofs_southj(k) =  min(jk) ;   % S extent as f(z)
                        %%%%%% 1marHOV

                        %updraft area column stuff
                        if( rem(numpts, 2)==0 & numpts >= thresh_updraft_minpts)  %if even number of numpts (and not a 1d updraft)
                            mza = (dx * 1000) .* ( (numpts/2) - 1) .* (dx * 1000)  ;   %updraft area at this ht [m^2]
                        elseif( rem(numpts, 2)==1 & numpts >= thresh_updraft_minpts )%if odd number of numpts (and not a 1d updraft)
                            mza = (dx * 1000) .* ( floor(numpts/2) - 1) .* (dx * 1000)   +   ((dx * 1000).*(dx * 1000))*0.5    ;   %updraft area at this ht [m^2] -  (rounds down numpts, then adds in a triangle-half-area)
                        else
                            mza = 0.0;
                        end

                        mzr = vertcat( mzr, 3 * ( abs( (mza/3.14159).^0.5 )/(dx*1000) ) ) ;  % 3x radius [index space: 3*(R/dx), not geographical distance] of updraft
%                         updraftLFCLAYprofs_equivarea(m,k) = mza ;   % m^2
%                         updraftLFCLAYprofs_numij(m,k) = numpts;
                        updraftLFCLAYprofs_equivarea(k) = mza ;   % m^2
                        updraftLFCLAYprofs_numij(k) = numpts;

                        %updraft peak magnitude column
                        updraftLFCLAYprofs_maxmag(k) = max( max(  W3d(ik,jk,k)  ) ) ;   % m^2

                        %upward mass flux at each height:
                        updraftLFCLAYprofs_massflux(k) =   mean( mean(  VM(ik,jk,k), 'omitnan'  ) ) .* mza ;  %making to sure to use full VM, not clearair VM

                        % figure; contourf(VM(:,:,10,15),20,'LineColor','none');  hold on; contour(Wcld_all(:,:,10,15),[2.5:1:11.5],'k')

                        %stuff for local LFC detection and updraft placement
                        updraftLFCLAYprofs_centlat(k)   = mean ( mean( Lat(ik,jk), 'omitnan' ) ) ;
                        updraftLFCLAYprofs_centlon(k)   = mean ( mean( Lon(ik,jk), 'omitnan' ) ) ;
                        updraftLFCLAYprofs_centHTagl(k) = mean ( mean( Ht_agl(ik,jk,k), 'omitnan' ) ) ;
                        mzi = vertcat(   mzi, mean(ik, 'omitnan' ) ); % mean i of updraft mask (f(z))
                        mzj = vertcat(   mzj, mean(jk, 'omitnan' ) ); % mean j of updraft mask (f(z))

                    end

                end  %kloop


                % now find updraft's local MULFC_agl within area surronding updraft centroid (using mean i/j of mask - rounded down)
                mki =  floor(mean( mzi,'omitnan') ) ; % z-mean i of updraft mask
                mkj =  floor(mean( mzj,'omitnan') ) ; % z-mean j of updraft mask
                mkr =  floor(max( mzr,[],'omitnan') ) ;

                %indices of averaging box surrounding the updraft
                meanbox_i1 = mki - mkr;
                meanbox_i2 = mki + mkr;
                meanbox_j1 = mkj - mkr;
                meanbox_j2 = mkj + mkr;

                %edge of horiz domain check - if goes outside of domina, just set it to edge of domain
                if( meanbox_i1 < 1)
                    meanbox_i1 = 1;
                end
                if( meanbox_j1 < 1)
                    meanbox_j1 = 1;
                end
                if( meanbox_i2 > Xs )
                    meanbox_i2 = Xs;
                end
                if( meanbox_j2 > Ys)
                    meanbox_j2 = Ys;
                end

                %mean LFC (m agl) surrounding updraft:
                updrafts_mulfcagl = mean( mean( mulfcagl(meanbox_i1:meanbox_i2,meanbox_j1:meanbox_j2), 'omitnan' ),'omitnan' ) ;


                if( isnan(updrafts_mulfcagl) == 0)
                    %agl height index of this LFC (closest k):
                    mk_mulfcagl   = find( abs(Ht -  updrafts_mulfcagl) == min(abs(Ht - updrafts_mulfcagl) )  )  ;
                    mk_mulfcaglm1 = find( abs(Ht - (updrafts_mulfcagl - kmbelowLFC ) ) == min(   abs( Ht - (updrafts_mulfcagl - kmbelowLFC) ) )  ) ;
                    mk_mulfcaglp1 = find( abs(Ht - (updrafts_mulfcagl + kmaboveLFC ) ) == min(   abs( Ht - (updrafts_mulfcagl + kmaboveLFC) ) )  ) ;

                    %edge of vert domain check - if goes outside of domina, just set it to edge of domain
                    if( mk_mulfcaglm1 < 1 )
                        mk_mulfcaglm1 = 1;
                    end
                    if( mk_mulfcaglp1 > Zs )
                        mk_mulfcaglp1 = Zs;
                    end

                    updraftLFCLAYprofs_klfcagl = mk_mulfcagl;
                    updraftLFCLAYprofs_kbotlayer = mk_mulfcaglm1;
                    updraftLFCLAYprofs_ktoplayer = mk_mulfcaglp1;
                end  %updrafts_mulfcagl = nan check

            end % wht_robust_flag:  there's at least one ht with an updraft with at least 3 adjoining pts.

        end %numel(wimaxs check)


	%disp(' done with m loop')

    %fill Final output vars
    FupdraftLFCLAYprofs_equivarea(m,:)    = updraftLFCLAYprofs_equivarea ;
    FupdraftLFCLAYprofs_numij(m,:)        = updraftLFCLAYprofs_numij ;
    FupdraftLFCLAYprofs_maxmag(m,:)       = updraftLFCLAYprofs_maxmag ;
    FupdraftLFCLAYprofs_massflux(m,:)     = updraftLFCLAYprofs_massflux ;
    FupdraftLFCLAYprofs_centlat(m,:)      = updraftLFCLAYprofs_centlat ;
    FupdraftLFCLAYprofs_centlon(m,:)      = updraftLFCLAYprofs_centlon ;
    FupdraftLFCLAYprofs_centHTagl(m,:)    = updraftLFCLAYprofs_centHTagl ;

    FupdraftLFCLAYprofs_kbotlayer(m)    = updraftLFCLAYprofs_kbotlayer ;
    FupdraftLFCLAYprofs_ktoplayer(m)    = updraftLFCLAYprofs_ktoplayer ;
    FupdraftLFCLAYprofs_klfcagl(m)      = updraftLFCLAYprofs_klfcagl ;
%     FupdraftLFCLAYprofs_htasl_wtop   = updraftLFCLAYprofs_htasl_wtop ;
%     FupdraftLFCLAYprofs_htasl_wbot   = updraftLFCLAYprofs_htasl_wbot ;

    %%%%% 1marHOV
    FupdraftLFCLAYprofs_northj(m,:)    =  updraftLFCLAYprofs_northj ;
    FupdraftLFCLAYprofs_southj(m,:)    =  updraftLFCLAYprofs_southj ;
    %%%%% 1marHOV

    %disp( ['tt',num2str(tt),' DONE m',num2str(m)] )


    end  %updraft mask par loop (m)
    % disp('333')

    toc



    disp('done with parallel')

    %rename for contiuned use:
    updraftLFCLAYprofs_equivarea = FupdraftLFCLAYprofs_equivarea;       clear FupdraftLFCLAYprofs_equivarea
    updraftLFCLAYprofs_numij = FupdraftLFCLAYprofs_numij;               clear FupdraftLFCLAYprofs_numij
    updraftLFCLAYprofs_maxmag = FupdraftLFCLAYprofs_maxmag;             clear FupdraftLFCLAYprofs_maxmag
    updraftLFCLAYprofs_massflux = FupdraftLFCLAYprofs_massflux;         clear FupdraftLFCLAYprofs_massflux
    updraftLFCLAYprofs_centlat = FupdraftLFCLAYprofs_centlat;           clear FupdraftLFCLAYprofs_centlat
    updraftLFCLAYprofs_centlon = FupdraftLFCLAYprofs_centlon;           clear FupdraftLFCLAYprofs_centlon
    updraftLFCLAYprofs_centHTagl = FupdraftLFCLAYprofs_centHTagl;       clear FupdraftLFCLAYprofs_centHTagl
    updraftLFCLAYprofs_kbotlayer = FupdraftLFCLAYprofs_kbotlayer;       clear FupdraftLFCLAYprofs_kbotlayer
    updraftLFCLAYprofs_ktoplayer = FupdraftLFCLAYprofs_ktoplayer;       clear FupdraftLFCLAYprofs_ktoplayer
    updraftLFCLAYprofs_klfcagl = FupdraftLFCLAYprofs_klfcagl;           clear FupdraftLFCLAYprofs_klfcagl
    %%%%1marHOV
    updraftLFCLAYprofs_northj   =  FupdraftLFCLAYprofs_northj;           clear FupdraftLFCLAYprofs_northj
    updraftLFCLAYprofs_southj   =  FupdraftLFCLAYprofs_southj;           clear FupdraftLFCLAYprofs_southj
    

    %  find(updraftLFCLAYprofs_kbotlayer > kagltop)
    %  find(updraftLFCLAYprofs_ktoplayer > kagltop)

    % %     % I found a few points in a few days where
    % %     % updraftLFCLAYprofs_ktoplayer(m,t) > kagltop;  Manually reset
    % %  those to kagltop.... nevermind, just raised kagltop
    % %     for m = 1:maxnum_guess
    % %         for t = 1 : sa-1
    % %             updraftLFCLAYprofs_kbotlayer(m,t) > kagltop
    % %             updraftLFCLAYprofs_ktoplayer(m,t) > kagltop
    % %         end
    % %     end




    %%%%%% note: at this point, the raw updraft locations and metric profiles that you'll want are done,
    %%%%%% but you'll want to clean up and filter some out next:


    %%%%%
    %%%%% cleanup updraft metric arrays:
    %%%%%


    
    % %backup before you start manipulating:
    %pres_updraftLFCLAYprofs_centHTagl = updraftLFCLAYprofs_centHTagl;
    prefilt_updraftLFCLAYprofs_centlat = updraftLFCLAYprofs_centlat;
    prefilt_updraftLFCLAYprofs_centlon = updraftLFCLAYprofs_centlon;
    %pres_updraftLFCLAYprofs_equivarea = updraftLFCLAYprofs_equivarea;
    %prefilt_updraftLFCLAYprofs_ID = updraftLFCLAYprofs_ID;
    %pres_updraftLFCLAYprofs_massflux = updraftLFCLAYprofs_massflux;
    %pres_updraftLFCLAYprofs_maxmag = updraftLFCLAYprofs_maxmag;
    prefilt_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;
    %pres_updraftLFCLAYprofs_kbotlayer = updraftLFCLAYprofs_kbotlayer;
    %pres_updraftLFCLAYprofs_ktoplayer = updraftLFCLAYprofs_ktoplayer;
    %
    % % blah1 = pres_updraftLFCLAYprofs_numij(:,:,end);

    %%%%%%1marHOV
    prefilt_updraftLFCLAYprofs_northj  =  updraftLFCLAYprofs_northj  ; 
    prefilt_updraftLFCLAYprofs_southj  =  updraftLFCLAYprofs_southj  ;     

    % revive preservatives:
    %  updraftLFCLAYprofs_centHTagl = pres_updraftLFCLAYprofs_centHTagl;  updraftLFCLAYprofs_centlat = pres_updraftLFCLAYprofs_centlat; updraftLFCLAYprofs_centlon = pres_updraftLFCLAYprofs_centlon; updraftLFCLAYprofs_equivarea = pres_updraftLFCLAYprofs_equivarea; updraftLFCLAYprofs_ID = pres_updraftLFCLAYprofs_ID;  updraftLFCLAYprofs_massflux = pres_updraftLFCLAYprofs_massflux; updraftLFCLAYprofs_maxmag = pres_updraftLFCLAYprofs_maxmag; updraftLFCLAYprofs_numij = pres_updraftLFCLAYprofs_numij; updraftLFCLAYprofs_kbotlayer = pres_updraftLFCLAYprofs_kbotlayer; updraftLFCLAYprofs_ktoplayer = pres_updraftLFCLAYprofs_ktoplayer;


    disp('filtering updraft profiles outside of targeted LFC layer')

    % [0] NaN-out heights above/below the LFC layer - to just keep the prescribed layer
    % this loop doesnt shrink the size of the arrays, just nans unwanted values
    %for t = 1:sa-1
    for m = 1:maxnum_guess

        if( isnan(updraftLFCLAYprofs_kbotlayer(m))==0  &  isnan(updraftLFCLAYprofs_ktoplayer(m))==0 )

            %disp(['m',num2str(m),'___updraftLFCLAYprofs_kbotlayer',num2str(updraftLFCLAYprofs_kbotlayer(m)),'___updraftLFCLAYprofs_ktoplayer',num2str(updraftLFCLAYprofs_ktoplayer(m))])

            updraftLFCLAYprofs_equivarea(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_equivarea(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_numij(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_numij(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_maxmag(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_maxmag(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_massflux(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_massflux(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_centlat(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_centlat(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_centlon(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_centlon(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_centHTagl(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_centHTagl(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            %%%%%%1marHOV
            updraftLFCLAYprofs_northj(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_northj(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;

            updraftLFCLAYprofs_southj(m,1:updraftLFCLAYprofs_kbotlayer(m))  = NaN; % -1111;
            updraftLFCLAYprofs_southj(m,updraftLFCLAYprofs_ktoplayer(m):end)  = NaN; % -1111;
            %%%%%%1marHOV

        end
    end  %mloop
    %end


    filt1_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;




    disp('filtering updrafts on lat/lon box')

    % [1] filter out those outside of lat/lon box that is near SDC:
    [al bl cl] = size(updraftLFCLAYprofs_centlat) ;
    %for t = 1:cl
    % t = 15;
    kill = [];
    for m = 1:al
        % m = 110; m = 150; t = 15;
        blah = [];
        %locate updrafts with column-mean centroids in the prescribed lat/lon box:
        blah =   mean( updraftLFCLAYprofs_centlat(m,:),2,'omitnan') > -33.0    &    mean( updraftLFCLAYprofs_centlat(m,:),2,'omitnan') < -31.5    &  ...
            mean( updraftLFCLAYprofs_centlon(m,:),2,'omitnan') > -65.1     &     mean( updraftLFCLAYprofs_centlon(m,:),2,'omitnan') < -64.6  & ...
            isnan(mean( updraftLFCLAYprofs_centlat(m,:),2,'omitnan'))==0  &  isnan(mean( updraftLFCLAYprofs_centlon(m,:),2,'omitnan'))==0 ;
        if(  blah == 0  )
            kill = vertcat(kill,m);
        end
    end % m
    updraftLFCLAYprofs_equivarea(kill,:)    = NaN; % -999;
    updraftLFCLAYprofs_numij(kill,:)        = NaN; % -999;
    updraftLFCLAYprofs_maxmag(kill,:)       = NaN; % -999;
    updraftLFCLAYprofs_massflux(kill,:)     = NaN; % -999;
    updraftLFCLAYprofs_centlat(kill,:)      = NaN; % -999;
    updraftLFCLAYprofs_centlon(kill,:)      = NaN; % -999;
    updraftLFCLAYprofs_centHTagl(kill,:)    = NaN; % -999;
    updraftLFCLAYprofs_klfcagl(kill)    = NaN;
    %end %time

    %%%%%%1marHOV
    updraftLFCLAYprofs_northj(kill,:)    = NaN;
    updraftLFCLAYprofs_southj(kill,:)    = NaN;
    %end %time

    filt2_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;


    disp('filtering updraft hts that are too narrow')

    %%% [2] kill updrafts that don't have at least one height in their profile with area > 2 pixels at one height in their profile (and other junk profiles):
    [al bl cl] = size(updraftLFCLAYprofs_centlat) ;
    %for t = 1:cl
    killarea = [];
    for m = 1:al
        % m = 110; t = 15;
        ma = [];
        %look for the max num areal points in height:
        ma = max( updraftLFCLAYprofs_numij(m,:),[],2,'omitnan' ) ;
        if( isnan(ma) | isempty(ma) | ma < thresh_updraft_minpts )
            killarea = vertcat(killarea,m) ;
        end
    end
    updraftLFCLAYprofs_equivarea(killarea,:)    = NaN; % -888;
    updraftLFCLAYprofs_numij(killarea,:)        = NaN; % -888;
    updraftLFCLAYprofs_maxmag(killarea,:)       = NaN; % -888;
    updraftLFCLAYprofs_massflux(killarea,:)     = NaN; % -888;
    updraftLFCLAYprofs_centlat(killarea,:)      = NaN; % -888;
    updraftLFCLAYprofs_centlon(killarea,:)      = NaN; % -888;
    updraftLFCLAYprofs_centHTagl(killarea,:)    = NaN; % -888;
    updraftLFCLAYprofs_klfcagl(killarea)        = NaN;
    %end
        %%%%%%1marHOV
    updraftLFCLAYprofs_northj(killarea,:)    = NaN;
    updraftLFCLAYprofs_southj(killarea,:)    = NaN;

    filt3_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;



    %%% [3] filter out repeated 3d w masks that happen because of multiple maxima in exact same 3d updraft shell
    disp('filtering repeated 3d updraft shells')

    [al bl cl] = size(updraftLFCLAYprofs_centlat) ;
    %for t = 1:cl
    %t = 15;
    kill_repeats = [];
    % look for masks with exact same mean properties and only keep the 1st one:
    meanlat =  max( updraftLFCLAYprofs_centlat(:,:) ,[], 2,  'omitnan')  ;
    meanlon =  max( updraftLFCLAYprofs_centlon(:,:),[], 2,  'omitnan')  ;
    meanagl =  max( updraftLFCLAYprofs_centHTagl(:,:) ,[], 2,  'omitnan')  ;
    meanmaxw =  max( updraftLFCLAYprofs_maxmag(:,:),[], 2,  'omitnan')  ;
    meanarea =  max( updraftLFCLAYprofs_equivarea(:,:) ,[], 2,  'omitnan')  ;
    for m = 1 : al - 1
        % m = find(imaxs == 104 & jmaxs == 190) ;   %tester
        rep = find( meanlat(m) == meanlat(m+1:end)  &  meanlon(m) == meanlon(m+1:end)  &  meanagl(m) == meanagl(m+1:end)  &  meanmaxw(m) == meanmaxw(m+1:end)  &   meanarea(m) == meanarea(m+1:end)  ) + m;
        if( isempty(rep)==0 )
            kill_repeats = cat(1,kill_repeats,rep(1:end));
        end
    end
    unique_kills = unique( kill_repeats ) ;  % these should be the repeats, test with 2D composite updraft number maps and i,jmaxes plotted
    updraftLFCLAYprofs_equivarea(unique_kills,:)    = NaN; % -777;
    updraftLFCLAYprofs_numij(unique_kills,:)        = NaN; % -777;
    updraftLFCLAYprofs_maxmag(unique_kills,:)       = NaN; % -777;
    updraftLFCLAYprofs_massflux(unique_kills,:)     = NaN; % -777;
    updraftLFCLAYprofs_centlat(unique_kills,:)      = NaN; % -777;
    updraftLFCLAYprofs_centlon(unique_kills,:)      = NaN; % -777;
    updraftLFCLAYprofs_centHTagl(unique_kills,:)    = NaN; % -777;
    updraftLFCLAYprofs_klfcagl(unique_kills)        = NaN;
    %end

    %%%%%%1marHOV
    updraftLFCLAYprofs_northj(unique_kills,:)    = NaN;
    updraftLFCLAYprofs_southj(unique_kills,:)    = NaN;

    %notes: a) equivarea <= 0 means <4 numij pts, but the updraft is kept because somewhere in the profile there is numij >=4



    filt4_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;




    disp('shrinking updraft prof arrays')
    %%% [4] shrink the array sizes for memory sake (and log the Ws' ID # to correspond with mask numbers):
    %define how many updrafts to keep (with some buffer) to define dimensions of final cleaned array

    howmany = [];
    %for t = 1:cl
    %using numij field as proxy to look for the keeper updrafts (left after cleanup)
    may =  max(updraftLFCLAYprofs_numij(:,:),[],2,'omitnan');
    may(may < 0 ) = []; may(isnan(may)) = [];
    howmany = vertcat(howmany,length(may));
    %end
    newnum = max(howmany)+5;   %largest num of updrafts in all times (+5 for a buffer)

    % initialize shrunken fields
    updraftLFCLAYprofs_ID                = zeros(newnum, 1);         updraftLFCLAYprofs_ID(:) = NaN;  % ID number of updraft on original (unshrunken) array - not super relevant if we're not tracking updrafts, but it could come in handy in other ways (?)
    shrink_updraftLFCLAYprofs_equivarea  = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_equivarea(:) = NaN;
    shrink_updraftLFCLAYprofs_numij      = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_numij(:) = NaN;  %numbe of horiz grid pts per height
    shrink_updraftLFCLAYprofs_maxmag     = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_maxmag(:) = NaN;
    shrink_updraftLFCLAYprofs_massflux   = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_maxmag(:) = NaN;
    shrink_updraftLFCLAYprofs_centlat    = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_centlat(:) = NaN;
    shrink_updraftLFCLAYprofs_centlon    = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_centlon(:) = NaN;
    shrink_updraftLFCLAYprofs_centHTagl  = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_centHTagl(:) = NaN;
    shrink_updraftLFCLAYprofs_klfcagl    = zeros(newnum,1);          shrink_updraftLFCLAYprofs_klfcagl(:) = NaN;
%%%%%%1marHOV
    shrink_updraftLFCLAYprofs_northj  = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_northj(:) = NaN;
    shrink_updraftLFCLAYprofs_southj  = zeros(newnum,kagltop);    shrink_updraftLFCLAYprofs_southj(:) = NaN;    
    
    
    %for t = 1:cl
    mm = 1;  %index counter, leave as = 1 here
    for m = 1:al
        keep = find( isnan(updraftLFCLAYprofs_numij(m,:))==0  & updraftLFCLAYprofs_numij(m,:) > 0) ;   %keep the just the updrafts that actually have info left in them (just going to use numij as the metric)
        if(isempty(keep)==0)
            %mm
            updraftLFCLAYprofs_ID(mm)                  =    m  ;
            shrink_updraftLFCLAYprofs_equivarea(mm,:)  = updraftLFCLAYprofs_equivarea(m,:);
            shrink_updraftLFCLAYprofs_numij(mm,:)      = updraftLFCLAYprofs_numij(m,:);
            shrink_updraftLFCLAYprofs_maxmag(mm,:)     = updraftLFCLAYprofs_maxmag(m,:);
            shrink_updraftLFCLAYprofs_massflux(mm,:)   = updraftLFCLAYprofs_massflux(m,:);
            shrink_updraftLFCLAYprofs_centlat(mm,:)    = updraftLFCLAYprofs_centlat(m,:);
            shrink_updraftLFCLAYprofs_centlon(mm,:)    = updraftLFCLAYprofs_centlon(m,:);
            shrink_updraftLFCLAYprofs_centHTagl(mm,:)  = updraftLFCLAYprofs_centHTagl(m,:);
            shrink_updraftLFCLAYprofs_klfcagl(mm)      = updraftLFCLAYprofs_klfcagl(m);
        %%%%%%1marHOV
            shrink_updraftLFCLAYprofs_northj(mm,:)   =  updraftLFCLAYprofs_northj(m,:);   
            shrink_updraftLFCLAYprofs_southj(mm,:)   =  updraftLFCLAYprofs_southj(m,:);      
        %%%%%%1marHOV	    
            mm = mm + 1;
        end
    end
    %end

    % rename the shrunken fields back to what I would rather name them:
    clear  updraftLFCLAYprofs_equivarea  updraftLFCLAYprofs_numij  updraftLFCLAYprofs_maxmag  updraftLFCLAYprofs_massflux
    clear  updraftLFCLAYprofs_centlat  updraftLFCLAYprofs_centlon  updraftLFCLAYprofs_centHTagl updraftLFCLAYprofs_klfcagl
    updraftLFCLAYprofs_equivarea = shrink_updraftLFCLAYprofs_equivarea ;
    updraftLFCLAYprofs_numij = shrink_updraftLFCLAYprofs_numij ;
    updraftLFCLAYprofs_maxmag = shrink_updraftLFCLAYprofs_maxmag ;
    updraftLFCLAYprofs_massflux = shrink_updraftLFCLAYprofs_massflux ;
    updraftLFCLAYprofs_centlat = shrink_updraftLFCLAYprofs_centlat ;
    updraftLFCLAYprofs_centlon = shrink_updraftLFCLAYprofs_centlon ;
    updraftLFCLAYprofs_centHTagl = shrink_updraftLFCLAYprofs_centHTagl ;
    updraftLFCLAYprofs_klfcagl = shrink_updraftLFCLAYprofs_klfcagl ;
    %%%%%%1marHOV
    updraftLFCLAYprofs_northj   =  shrink_updraftLFCLAYprofs_northj; 
    updraftLFCLAYprofs_southj   =  shrink_updraftLFCLAYprofs_southj;   
    %%%%%%1marHOV

    clear shrink_updraftLFCLAYprofs_equivarea shrink_updraftLFCLAYprofs_numij shrink_updraftLFCLAYprofs_maxmag shrink_updraftLFCLAYprofs_massflux
    clear shrink_updraftLFCLAYprofs_centlat shrink_updraftLFCLAYprofs_centlon shrink_updraftLFCLAYprofs_centHTagl shrink_updraftLFCLAYprofs_klfcagl
    clear shrink_updraftLFCLAYprofs_northj shrink_updraftLFCLAYprofs_southj



    filt5_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;


    disp('calculating a ht_asl for updraft prof layers')
    % calcualte ht asl within each updraft profile, based on htagl/lat/lon @centroid(z) - terrain there.

    [nl kl tl] = size(updraftLFCLAYprofs_centlat) ;
    updraftLFCLAYprofs_centHTasl = updraftLFCLAYprofs_centHTagl;
    %for t = 1:tl
    for m = 1:nl
        for k = 1:kl
            %   t = 1;  k = 21;  m = 1

            %  updraftLFCLAYprofs_centlat(m,k,t)
            %  updraftLFCLAYprofs_centlon(m,k,t)

            if( isnan(updraftLFCLAYprofs_centlat(m,k)) ==0 )   % updraftLFCLAYprofs_centlon(m,k,t)

                latdiff = abs( Lat - updraftLFCLAYprofs_centlat(m,k) );
                londiff = abs( Lon - updraftLFCLAYprofs_centlon(m,k) );
                diff = latdiff + londiff;
                [ilat ilon] = find( diff == min(min(diff)) )  ;
                if( length(ilat)> 0 )
                    ilat = ilat(1);   ilon = ilon(1);
                end
                updraftLFCLAYprofs_centHTasl(m,k) = updraftLFCLAYprofs_centHTagl(m,k) + Terr(ilat,ilon);
            else
                updraftLFCLAYprofs_centHTasl(m,k) = NaN;
            end
        end
    end
    %end




    filt6_updraftLFCLAYprofs_numij = updraftLFCLAYprofs_numij;




%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OPTIONAL: reprocess exact same way, but this time to save these now thinned/QC'ed W 3D masks:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %lowest-level updraft masks:

    disp(' working on saving w shell masks')

    % seed the thinned W mask field:
    [um blah] = size(  updraftLFCLAYprofs_ID ) ;
    [ux uy uz] = size(  Wcld_all ) ;
    Wll_2d_masks = zeros(ux,uy,um); %[x,y,t,updraft#]   %the final thinned answer - 2D (vertical composite, below kagltop) masks
    %W_3d_masks = zeros(ux,uy,uz,ut,um); %[x,y,z,t,updraft#]   %the final thinned answer - 3D masks

    %for t = 1:sa-1

    %start by filtering 3d low-level cloudy w field on low w threshold:
    W3d = Wcld_all(:,:,1:kagltop);
    W3d(W3d < wlowthresh) = 0;
    W3d(isnan(W3d)) = 0;
    Wcomp = max(W3d,[],3);

    %checkpoint
    if(  isempty( find(Wcomp ~= 0) )  )  %if all updrafts are filtered out (if = 1), then skip this time in t-loop
        continue  %skip time
    end

    % then finding local w-maxima in 3D field; practially, doing it in a vertically composited 2d areal plane
    wmaxes = imregionalmax(Wcomp,8) ;
    [wimaxs wjmaxs] = find(wmaxes == 1) ;

    %checkpoint
    if(  isempty(wimaxs)  )  %if all w-max points filtered out (if = 1), then skip this time in loop
        continue    %skip time
    end

    % now find the height of each i,j max:
    wkmaxs = wimaxs; wkmaxs(:) = NaN;
    for m = 1:length(wimaxs)
        % m = find(imaxs == 104 & jmaxs == 190) ;   %tester
        wprof = W3d(wimaxs(m),wjmaxs(m),:); wprof = permute(wprof,[3 1 2]);
        wkmaxs(m) = find( wprof == max(wprof) );
    end

    %logic-ify the w field:
    W3dbin = W3d;
    W3dbin(W3d(:,:,:) >  wlowthresh) = 1;
    W3dbin(W3d(:,:,:) <= wlowthresh) = 0;
    W3dbin = logical(W3dbin);
    [aa bb cc]= size(W3dbin);  %i,j,k

    %loop thru just the kept (after shrinking arrays) updraft IDs
    kept_w_masks = updraftLFCLAYprofs_ID(:);  kept_w_masks(isnan(kept_w_masks)) = [];
    for m = 1:um
        if( m <= length(kept_w_masks) )
            % flood fill 3d binary updraft object
            Wobj3d = bwselect3(W3dbin,wjmaxs(kept_w_masks(m)),wimaxs(kept_w_masks(m)),wkmaxs(kept_w_masks(m)),26) ;
            Wobj3di = double( Wobj3d);   Wobj3di(Wobj3di==1) = kept_w_masks(m) ;

            %FINAL ANSWER - it is the mask number in 2D vertical composite mask.
            Wll_2d_masks(:,:,m) = max(Wobj3di,[],3);

            %put a 3d mask in here somewhere if you eventually need that... gulp
            %matlab did no like at ~300x300x40x30x100 sized array, so if you want this, you'll have to trim things more!

        end
    end

    %end
    toc
    % done with 3D wID
%}



    %RENAME:
    wIDprofs_equivarea = updraftLFCLAYprofs_equivarea;
    wIDprofs_numij = updraftLFCLAYprofs_numij;
    wIDprofs_maxmag = updraftLFCLAYprofs_maxmag;
    wIDprofs_massflux = updraftLFCLAYprofs_massflux;
    wIDprofs_centlat = updraftLFCLAYprofs_centlat;
    wIDprofs_centlon = updraftLFCLAYprofs_centlon;
    wIDprofs_centHTagl = updraftLFCLAYprofs_centHTagl;
    wIDprofs_centHTasl = updraftLFCLAYprofs_centHTasl;
    wIDprofs_ID = updraftLFCLAYprofs_ID;
    wIDprofs_zagl = updraftLFCLAYprofs_zagl;
    % wIDprofs_2d_masks = Wll_2d_masks;
    wIDprofs_klfcagl = updraftLFCLAYprofs_klfcagl;
    %%%%1marHOV
    wIDprofs_northj = updraftLFCLAYprofs_northj;
    wIDprofs_southj = updraftLFCLAYprofs_southj;


    clear updraftLFCLAYprofs_equivarea updraftLFCLAYprofs_numij updraftLFCLAYprofs_maxmag updraftLFCLAYprofs_massflux updraftLFCLAYprofs_centlat ...
        updraftLFCLAYprofs_centlon updraftLFCLAYprofs_centHTagl updraftLFCLAYprofs_ID updraftLFCLAYprofs_zagl Wll_2d_masks updraftLFCLAYprofs_klfcagl ...
        updraftLFCLAYprofs_northj updraftLFCLAYprofs_southj

    disp('     ')
    disp(' saving wID output into mat file: ' )
    disp(matoutsave)
    save(matoutsave,'wrfdir','Lon','Lat','W_all','Wcld_all','Terr','dx','HHMMSS_all','YYMMDD_all',...
        'wIDprofs_equivarea','wIDprofs_numij','wIDprofs_maxmag','wIDprofs_massflux','wIDprofs_centlat',...
        'wIDprofs_centlon','wIDprofs_centHTagl','wIDprofs_centHTasl','wIDprofs_ID','wIDprofs_zagl',...
        'wIDprofs_klfcagl','wIDprofs_northj','wIDprofs_southj',...
	'filt6_updraftLFCLAYprofs_numij','filt5_updraftLFCLAYprofs_numij','filt4_updraftLFCLAYprofs_numij','filt3_updraftLFCLAYprofs_numij','filt2_updraftLFCLAYprofs_numij','filt1_updraftLFCLAYprofs_numij',...
	'prefilt_updraftLFCLAYprofs_centlat','prefilt_updraftLFCLAYprofs_centlon','prefilt_updraftLFCLAYprofs_numij','-v7.3')
    toc
    % 'wIDprofs_2d_masks',

    disp('    ')
    disp( ['DONEski with WID time: ',num2str(tt) ])
    


end   %serial time loop





%%

















