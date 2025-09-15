

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%% load cell/amf stats file
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

netcdf_out = 1;
savemats = 0;

% location of'raw' interpsonde data:
%rootdir = '/global/cfs/projectdirs/m1657/jmarquis/EnvironmentParameters/matlab/ZheLASSO_2prof_moistLFC/';


inputs = ls('/pscratch/sd/j/jmarquis/cacti/lasso/d4_100m_5min/22Jan100m/20190122/gefs18/base/les/subset_d4/f15min/wID_3Denv_*mat')
%inputs = ls('/Users/marq789/Documents/PROJECTS/ICLASS/moist_lfc_code_ZHelasso2loc/input/20181129/*/base/d2/stats_env1d_2location_20181129.1200_20181130.0000.nc');
inputs = split(inputs);
inputs = cell(inputs);
inputs(end,:) = []

[inlen blah] = size(inputs)


outdirprefix = '/pscratch/sd/j/jmarquis/testoutwidenv/';



%inlen = 5;

for in = 1 : inlen

    filname = string(char( inputs(in,:) ) )   

    splidirs = split(filname,"/") ;   splidirs(1) = [];
    outdir = strcat( outdirprefix, splidirs(9), '/', splidirs(10), '/', splidirs(11), '/', splidirs(12),'/' )  ;
    outdir = string(outdir);
    disp(outdir)
    mkdir(outdir)

    %disp('subfile')
    subfile = splidirs(15);  subfile = split(subfile,"_");
    widtime = subfile(4);

    %rootdir = '/Users/marq789/Documents/PROJECTS/ICLASS/moist_lfc_code_ZHelasso2loc/';
    %filname = [rootdir,'stats_env1d_2location_20181129.1200_20181130.0000.nc']

    %where you want the output:
    %outdir = [rootdir,'output/']
%     disp(outdir)
%     mkdir(outdir)


    load(filname,'WID_3Dtemper','WID_3Dqvapor','WID_3Dpress','WID_3Du','WID_3Dv','WID_3Dhtasl','WID_3Dhtagl','WID_3Drh')

    %size(WID_3Dpress)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % load in data arrays:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tdry_IN  = WID_3Dtemper;   %temperature (K)
    mr_IN    = WID_3Dqvapor;   %  vapor mix ratio [kg/kg]
    p_IN     = WID_3Dpress;    % pressure [mb]
    p_IN     = p_IN * 100;     % pressure [Pa]
    u_IN     = WID_3Du;        % [m/s]
    v_IN     = WID_3Dv;        % [m/s]
    htasl_IN = WID_3Dhtasl;    %model ht asl
    rh_IN    = WID_3Drh;       %rh 0-100%

    %tdry_IN = ncread(filname,'temperature');   %temperature (K)
    %mr_IN = ncread(filname,'qv'); %  vapor mix ratio [kg/kg]
    %p_IN = ncread(filname,'pressure'); % pressure [mb]
    %p_IN = p_IN * 100;   % pressure [Pa]
    %u_IN = ncread(filname,'u');  % [m/s]
    %v_IN = ncread(filname,'v'); % [m/s]
    %htasl_IN = ncread(filname,'height');  %model ht asl
    %rh_IN = ncread(filname,'rh');  %r 0-100%




    %permute to fit in right shaped hole:
    tdry_IN = permute(tdry_IN,[2 4 3 1]);
    mr_IN   = permute(mr_IN,[2 4 3 1]);
    p_IN    = permute(p_IN,[2 4 3 1]);
    u_IN    = permute(u_IN,[2 4 3 1]);
    v_IN    = permute(v_IN,[2 4 3 1]);
    rh_IN    = permute(rh_IN,[2 4 3 1]);
    htasl_IN    = permute(htasl_IN,[2 4 3 1]);

    % NL = locations, NT = times, MK = vert levels, NC = num cell tracks
    %[Nlocations Nklevs Ntimes Ncells] = size(tdry_IN)

    [Nlocations Ntimes Nklevs Ncells] = size(tdry_IN)


    % %convert gph to geometric height asl
    %Rea = 6.37e6;  %m  mean earth rad
    %grav = 9.8;   %constant m2/s2
    %%breaking up location 1 & 2 because Adam's current input ahs one GP and one GPH  [FIX THIS WHEN HE FIXES HIS INPUT]:
    %htasl_IN = gph; htasl_IN(:) = NaN;
    %htasl_IN(1,:,:,:) = (gph(1,:,:,:) .* Rea)./(Rea  - gph(1,:,:,:)) ;   %asl height
    %htasl_IN(2,:,:,:) = (gph(2,:,:,:) .* Rea)./(grav*Rea  - gph(2,:,:,:)) ;   %asl height
    %%htagl = htasl - htasl(1);  %% ht is now m AGL  (nvm, you do this later).





    %calc some other things:
    pres = p_IN/100;
    T = tdry_IN ;

    %convert mr -> rh,TdC
    E = 0.622;
    esat = 6.112*exp( (17.67*(tdry_IN-273.15))./(243.5+(tdry_IN-273.15)) )*100; % pa
    e = ((mr_IN.*p_IN)/E)./(1 + (mr_IN/E) );   %Pa
    rh_IN = 100*(e./esat); %  0-100percent
    blorg = log(e/611.2) ;
    TdC_IN =  (blorg .* 243.5) ./ (17.67 - blorg) ;  %degC

    clear blorg esat e T  gph Rea

       dummy = zeros(Nlocations, Ntimes, 1, Ncells);   dummy(:) = NaN;
      % [ds1 ds2 ds3 ds4]  = size(dummy)
       

       CAPE_mu = dummy;   %
        %    CIN_IB_mu = dummy;    %
        LCL_height_mu = dummy;
        LFC_height_mu= dummy;
        EL_height_mu= dummy;
        CIN_NA_fract_mu = dummy;
        CIN_mu = dummy;
        CAPEacbl_mu = dummy;
        %    CAPElcl_IB_mu = dummy;
        initial_ht_parcel_mu = dummy;
        pCAPE_mu = dummy;
        tallenough_mu  = dummy;
        EL_temp_mu  = dummy;
        LCL_temp_mu   = dummy;
        LFC_temp_mu   = dummy;
        pCIN_mu     = dummy;
        pCAPEacbl_mu   = dummy;
        pLFC_height_mu = dummy;
        pEL_height_mu  = dummy;
        pLFC_temp_mu   = dummy;
        pEL_temp_mu    = dummy;
        ptallenough_mu = dummy;
        CAPE_sfc = dummy;   %
        %    CIN_IB_sfc = dummy;    %
        LCL_height_sfc = dummy;
        LFC_height_sfc= dummy;
        EL_height_sfc= dummy;
        %     LCL_pres_sfc = dummy;
        %     LFC_pres_sfc= dummy;
        %     EL_pres_sfc= dummy;
        %    CIN_NA_fract_sfc= dummy;
        CIN_sfc = dummy;
        CAPEacbl_sfc = dummy;
        %    CAPElcl_IB_sfc = dummy;
        %    tallenough_sfc = dummy;
        pCAPE_sfc = dummy;
        %    pCIN_NA_sfc = dummy;
        tallenough_sfc  = dummy;
        EL_temp_sfc     = dummy;
        LCL_temp_sfc    = dummy;
        LFC_temp_sfc    = dummy;
        pCIN_sfc     = dummy;
        pCAPEacbl_sfc   = dummy;
        pLFC_height_sfc = dummy;
        pEL_height_sfc  = dummy;
        pLCL_temp_sfc   = dummy;
        pLFC_temp_sfc   = dummy;
        pEL_temp_sfc    = dummy;
        ptallenough_sfc = dummy;
       U_mu = dummy;
        V_mu = dummy;
        theta_mu = dummy;
        thetav_mu = dummy;
        CAPE_ml= dummy;   %
        %    CIN_IB_ml = dummy;    %
        LCL_height_ml = dummy;
        LFC_height_ml = dummy;
        EL_height_ml = dummy;
        %     LCL_pres_ml = dummy;
        %     LFC_pres_ml = dummy;
        %     EL_pres_ml = dummy;
        %    CIN_NA_fract_ml = dummy;
        CIN_ml = dummy;
        CAPEacbl_ml = dummy;
        %    CAPElcl_IB_ml = dummy;
        %    tallenough_ml = dummy;
        pCAPE_ml = dummy;
        %    pCIN_NA_ml = dummy;
        tallenough_ml  = dummy;
        EL_temp_ml     = dummy;
        LCL_temp_ml    = dummy;
        LFC_temp_ml    = dummy;
        pCIN_ml     = dummy;
        pCAPEacbl_ml   = dummy;
        pLFC_height_ml = dummy;
        pEL_height_ml  = dummy;
        pLCL_temp_ml   = dummy;
        pLFC_temp_ml   = dummy;
        pEL_temp_ml    = dummy;
        ptallenough_ml = dummy;
        PW  = dummy;
        freezing_ht = dummy;
        thetae_mu = dummy;
        thetae_mean_subcloud_mu = dummy;
        thetae_mean_ACBL_mu = dummy;
        %LCL_temp_mu = dummy;
        %LFC_temp_mu = dummy;
        %EL_temp_mu = dummy;
        RH_mean_ACBL = dummy;
        rvap_mu = dummy;
        rvap_850mb = dummy;
        rvap_700mb = dummy;
        rvap_500mb = dummy;
        rh_850mb = dummy;
        rh_700mb = dummy;
        rh_500mb = dummy;
        Ucrel_mu = dummy;
        Vcrel_mu = dummy;
        shear_mag_bulk_ACBL_mu = dummy;
        shear_mag_bulk_FT_mu = dummy;
        shear_dir_bulk_FT_mu = dummy;
        shear_mag_bulk_0to1km = dummy;
        shear_dir_bulk_0to1km = dummy;
        shear_mag_bulk_0to3km = dummy;
        shear_dir_bulk_0to3km = dummy;
        shear_mag_bulk_0to6km = dummy;
        shear_dir_bulk_0to6km = dummy;
        NANdiag = dummy;
        %     Ucrel_mean_ACBL_ml = dummy;
        %     Vcrel_mean_ACBL_ml = dummy;
        %     Ucrel_mean_0toACBL_ml = dummy;
        %     Vcrel_mean_0toACBL_ml = dummy;
        Ucrel_mean_ACBL_mu = dummy;
        Vcrel_mean_ACBL_mu = dummy;
        %    Ucrel_mean_0toACBL_mu = dummy;
        %    Vcrel_mean_0toACBL_mu = dummy;
        %     Ucrel_mean_ACBL_sfc = dummy;
        %     Vcrel_mean_ACBL_sfc = dummy;
        %     Ucrel_mean_0toACBL_sfc = dummy;
        %     Vcrel_mean_0toACBL_sfc = dummy;
        %     Ucrel_mean_EIL10 = dummy;
        %     Vcrel_mean_EIL10 = dummy;
        %     Ucrel_mean_EIL25 = dummy;
        %     Vcrel_mean_EIL25 = dummy;
        %     Ucrel_mean_EIL50 = dummy;
        %     Vcrel_mean_EIL50 = dummy;
        shear_mag_bulk_0to9km = dummy;
        shear_dir_bulk_0to9km = dummy;
        rvap_925mb = dummy;
        rh_925mb = dummy;
        rvap_600mb = dummy;
        rh_600mb = dummy;
        %Enoch vars:
        rh_at_muLFCplus1500m = dummy;
        rh_at_muLCLplus1500m = dummy;
        %adam vars
        rh_at_pmuLFCplus1500m = dummy;
        RH_mean_ACBL_pmu = dummy;
        shear_mag_bulk_ACBL_pmu = dummy;
        thetae_mean_ACBL_pmu = dummy;
        Ucrel_mean_ACBL_pmu = dummy;
        Vcrel_mean_ACBL_pmu = dummy;
        rvap_600mb = dummy;
        rvap_400mb = dummy;
        rvap_300mb = dummy;
        rh_600mb = dummy;
        rh_400mb = dummy;
        rh_300mb = dummy;
        rvap_3000masl = dummy;
        rvap_4000masl = dummy;
        rvap_5000masl = dummy;
        rvap_6000masl = dummy;
        rh_3000masl = dummy;
        rh_4000masl = dummy;
        rh_5000masl = dummy;
        rh_6000masl = dummy;


	upslope_max  = dummy;
	upslope_mean = dummy;
        upslope_fract = dummy;
	upslope_depth = dummy;


    todolist = [1:Ncells];  % [1:Ns];

    %todolist = [1:3];

    %parfor cccc = 1:length(todolist)      % cell track loop
    for cccc = 1:length(todolist)      % cell track loop
        c = todolist(cccc)
        % c = 1;



        labc = strcat('start cell ',num2str(c) );
        disp( labc )

        %     % calculate unique cell motion vectors:
        %     [cell_motion_x, cell_motion_y] = CELLindivmotion1(cellstats,40.0);
        %     %just current cell's motion:
        %     cell_motion_x = single( cell_motion_x(:,c) );
        %     cell_motion_y = single( cell_motion_y(:,c) );

        %%%%%%%%%%%%  WARNING: this is just temporary cell motion = 0 hack until you figure out advection file stuff (commented above)
        cell_motion_x = 0;  %zeros(Ntimes,1);
        cell_motion_y = 0;  %zeros(Ntimes,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %diagnostic used later when nco'ing files back together
        cellnum = c;


        %check = 1;

        %wall clock start time of the cell
        start_calculate_time  = datestr(now) ;

        %     %hack to look at a subset of data from 1-x2, 1-y2:
        %     x2 = Xs;
        %     y2 = Ys;

        % % % just load data for the current cell:
        %wrf_times = ncread(filname,'times');
        % % wrf_basetime = ncread(filname,'base_time');
        wrf_hasl = htasl_IN(:,:,:,c); %masl
        wrf_pres = p_IN(:,:,:,c)/100; %mb  (converts to pa where needed later)
        wrf_T = tdry_IN(:,:,:,c);     %K  (converts to C later)
        wrf_mr = mr_IN(:,:,:,c);      %kg/kg
        wrf_u = u_IN(:,:,:,c);        %m/s
        wrf_v = v_IN(:,:,:,c);        %m/s
        wrf_rh = rh_IN(:,:,:,c);      %0-100percent


        %ID label for cell output
        CELLID = pad(num2str(c),5,'left');  CELLID(find(CELLID==' '))='0';


        %intialize variables:
        dummy = single( zeros(Nlocations, Ntimes, 1) );    %note: 1 is dummy dimension to adapt to array uinput with fewer dimentions that we want the capability to do
        dummy(:) = -911;


%{
        CAPE_mu = dummy;   %
        %    CIN_IB_mu = dummy;    %
        LCL_height_mu = dummy;
        LFC_height_mu= dummy;
        EL_height_mu= dummy;
        CIN_NA_fract_mu = dummy;
        CIN_mu = dummy;
        CAPEacbl_mu = dummy;
        %    CAPElcl_IB_mu = dummy;
        initial_ht_parcel_mu = dummy;
        pCAPE_mu = dummy;
        tallenough_mu  = dummy;
        EL_temp_mu  = dummy;
        LCL_temp_mu   = dummy;
        LFC_temp_mu   = dummy;
        pCIN_mu     = dummy;
        pCAPEacbl_mu   = dummy;
        pLFC_height_mu = dummy;
        pEL_height_mu  = dummy;
        pLFC_temp_mu   = dummy;
        pEL_temp_mu    = dummy;
        ptallenough_mu = dummy;
        CAPE_sfc = dummy;   %
        %    CIN_IB_sfc = dummy;    %
        LCL_height_sfc = dummy;
        LFC_height_sfc= dummy;
        EL_height_sfc= dummy;
        %     LCL_pres_sfc = dummy;
        %     LFC_pres_sfc= dummy;
        %     EL_pres_sfc= dummy;
        %    CIN_NA_fract_sfc= dummy;
        CIN_sfc = dummy;
        CAPEacbl_sfc = dummy;
        %    CAPElcl_IB_sfc = dummy;
        %    tallenough_sfc = dummy;
        pCAPE_sfc = dummy;
        %    pCIN_NA_sfc = dummy;
        tallenough_sfc  = dummy;
        EL_temp_sfc     = dummy;
        LCL_temp_sfc    = dummy;
        LFC_temp_sfc    = dummy;
        pCIN_sfc     = dummy;
        pCAPEacbl_sfc   = dummy;
        pLFC_height_sfc = dummy;
        pEL_height_sfc  = dummy;
        pLCL_temp_sfc   = dummy;
        pLFC_temp_sfc   = dummy;
        pEL_temp_sfc    = dummy;
        ptallenough_sfc = dummy;
        U_mu = dummy;
        V_mu = dummy;
        theta_mu = dummy;
        thetav_mu = dummy;
        CAPE_ml= dummy;   %
        %    CIN_IB_ml = dummy;    %
        LCL_height_ml = dummy;
        LFC_height_ml = dummy;
        EL_height_ml = dummy;
        %     LCL_pres_ml = dummy;
        %     LFC_pres_ml = dummy;
        %     EL_pres_ml = dummy;
        %    CIN_NA_fract_ml = dummy;
        CIN_ml = dummy;
        CAPEacbl_ml = dummy;
        %    CAPElcl_IB_ml = dummy;
        %    tallenough_ml = dummy;
        pCAPE_ml = dummy;
        %    pCIN_NA_ml = dummy;
        tallenough_ml  = dummy;
        EL_temp_ml     = dummy;
        LCL_temp_ml    = dummy;
        LFC_temp_ml    = dummy;
        pCIN_ml     = dummy;
        pCAPEacbl_ml   = dummy;
        pLFC_height_ml = dummy;
        pEL_height_ml  = dummy;
        pLCL_temp_ml   = dummy;
        pLFC_temp_ml   = dummy;
        pEL_temp_ml    = dummy;
        ptallenough_ml = dummy;






        %     EIL10_top_height = dummy;
        %     EIL10_bot_height = dummy;
        %     EIL25_top_height = dummy;
        %     EIL25_bot_height = dummy;
        %     EIL50_top_height = dummy;
        %     EIL50_bot_height = dummy;

        %     rvap_mean_EIL50 = dummy;
        %     rvap_mean_EIL25 = dummy;
        %     rvap_mean_EIL10 = dummy;

        PW  = dummy;
        freezing_ht = dummy;

        thetae_mu = dummy;
        thetae_mean_subcloud_mu = dummy;
        thetae_mean_ACBL_mu = dummy;
        %LCL_temp_mu = dummy;
        %LFC_temp_mu = dummy;
        %EL_temp_mu = dummy;

        RH_mean_ACBL = dummy;
        rvap_mu = dummy;
        rvap_850mb = dummy;
        rvap_700mb = dummy;
        rvap_500mb = dummy;
        rh_850mb = dummy;
        rh_700mb = dummy;
        rh_500mb = dummy;

        Ucrel_mu = dummy;
        Vcrel_mu = dummy;

        shear_mag_bulk_ACBL_mu = dummy;
        shear_mag_bulk_FT_mu = dummy;
        shear_dir_bulk_FT_mu = dummy;

        shear_mag_bulk_0to1km = dummy;
        shear_dir_bulk_0to1km = dummy;
        shear_mag_bulk_0to3km = dummy;
        shear_dir_bulk_0to3km = dummy;
        shear_mag_bulk_0to6km = dummy;
        shear_dir_bulk_0to6km = dummy;

        NANdiag = dummy;

        %     Ucrel_mean_ACBL_ml = dummy;
        %     Vcrel_mean_ACBL_ml = dummy;
        %     Ucrel_mean_0toACBL_ml = dummy;
        %     Vcrel_mean_0toACBL_ml = dummy;

        Ucrel_mean_ACBL_mu = dummy;
        Vcrel_mean_ACBL_mu = dummy;
        %    Ucrel_mean_0toACBL_mu = dummy;
        %    Vcrel_mean_0toACBL_mu = dummy;

        %     Ucrel_mean_ACBL_sfc = dummy;
        %     Vcrel_mean_ACBL_sfc = dummy;
        %     Ucrel_mean_0toACBL_sfc = dummy;
        %     Vcrel_mean_0toACBL_sfc = dummy;

        %     Ucrel_mean_EIL10 = dummy;
        %     Vcrel_mean_EIL10 = dummy;
        %     Ucrel_mean_EIL25 = dummy;
        %     Vcrel_mean_EIL25 = dummy;
        %     Ucrel_mean_EIL50 = dummy;
        %     Vcrel_mean_EIL50 = dummy;

        shear_mag_bulk_0to9km = dummy;
        shear_dir_bulk_0to9km = dummy;

        %     U_4000m = dummy;
        %     V_4000m = dummy;

        rvap_925mb = dummy;
        rh_925mb = dummy;
        rvap_600mb = dummy;
        rh_600mb = dummy;

        %Enoch vars:
        rh_at_muLFCplus1500m = dummy;
        rh_at_muLCLplus1500m = dummy;

        %adam vars
        rh_at_pmuLFCplus1500m = dummy;
        RH_mean_ACBL_pmu = dummy;
        shear_mag_bulk_ACBL_pmu = dummy;
        thetae_mean_ACBL_pmu = dummy;
        Ucrel_mean_ACBL_pmu = dummy;
        Vcrel_mean_ACBL_pmu = dummy;

        rvap_600mb = dummy;
        rvap_400mb = dummy;
        rvap_300mb = dummy;
        rh_600mb = dummy;
        rh_400mb = dummy;
        rh_300mb = dummy;
        rvap_3000masl = dummy;
        rvap_4000masl = dummy;
        rvap_5000masl = dummy;
        rvap_6000masl = dummy;
        rh_3000masl = dummy;
        rh_4000masl = dummy;
        rh_5000masl = dummy;
        rh_6000masl = dummy;
%}









        %   check = 3;

        for x = 1:Nlocations
            %    x
            for y = 1:Ntimes
                for t = 1:1  %dummy var

                    %disp('****');
                    %disp('x')
                    %disp(num2str(x))
                    %disp('y')
                    %disp(num2str(y))
                    %disp('t')
                    %disp(num2str(t))

                    %   x = 2;   y = 3;   t = 1;

                    %YYMMDD = filname(107:114)
                    %wrf_base_time = ncread(filname,'base_time',[t c],[t c]); %[seconds since 1/1/1970 of 0000 utc - day of (i.e., single value per day)]

                    %wrf_base_time = 0.

                    %define the sounding from the 3D wrf chunk:
                    height = 0; pres = 0; tdry = 0; mr = 0; U = 0; V = 0; rh = 0;

                    height = wrf_hasl(x,y,:,t) - wrf_hasl(x,y,1,t);         %m AGL
                    pres = wrf_pres(x,y,:,t);           %mb
                    tdry = wrf_T(x,y,:,t) - 273.15;      %C
                    mr = wrf_mr(x,y,:,t) ;              %kg/kg
                    U = wrf_u(x,y,:,t);                 %m/s
                    V = wrf_v(x,y,:,t);                 %m/s
                    rh = wrf_rh(x,y,:,t);               %  0-100%

                    htasl = wrf_hasl(x,y,:,t);   %masl ht profile

                    %fix index order via permutation (so  vars are f(k) only):
                    U = permute(U,[3 2 1]);
                    V = permute(V,[3 2 1]);
                    ht = permute(height,[3 2 1]);
                    pres = permute(pres,[3 2 1]);
                    tdry = permute(tdry,[3 2 1]);
                    mr = permute(mr,[3 2 1]);
                    rh = permute(rh,[3 2 1]);

                    htasl = permute(htasl,[3 2 1]);

                    [sHt sTi] = size(tdry);


                    % figure; plot(ht,'or')

                    %%%%% chop out stray mr < 0 quirk, replace with interpolation:
                    deadpixels = find (  mr <= 0 );
                    mr(deadpixels) = NaN;
                    for nm = 1:length(deadpixels)
                        mr(deadpixels(nm)) = mean ( mean(  mr(deadpixels(nm)-1:deadpixels(nm)+1)  ,'omitnan') ) ;
                    end
                    %%%%%


                    if( isnan(height)==0 & isnan(pres)==0 & isnan(tdry)==0 & isnan(mr)==0 & isnan(U)==0 & isnan(V)==0 )

                        NANdiag(x,y,t,c) = 0;


                        %calc some other things:
                        p = pres*100;     %pres -> Pa
                        E = 0.622;
                        T = tdry + 273.15;       % temperature ->  K
                        f = mr / 0.622;
                        e = (p.*f)./(1+f);   %pa
                        sphu  = mr./(1+mr);  %specific humidity

                        % calculate potential temp:
                        theta = T.*(100000./p).^0.286; %K

                        %emanual dew point:
                        eee = e/100; %vap press [mb]
                        TdC =273.15+  243.5./( (17.67./log(eee./6.112)) -1) - 273.15;    %  deg C      K.Emanuel


                        %%% define some constants:
                        cpd = 1006; % J/kg K
                        g = 9.8;    % m/s2
                        cpv = 1870; %J/kg K
                        E = 0.622;
                        lv0 = 2501000;
                        Rv = 461.5;
                        Rd = 287.04;
                        cw = 4190;             % heat capacity of water
                        cc = 2320;
                        cl = 4200;
                        cvv = 1410;
                        cvd = 719;

                        alv = lv0 - cc*(T-273.15); %J/kg

                        Tv    = T.*(1+0.608.*mr);  %K
                        thetav = Tv.*(100000./p).^0.286; %K

                        % teA = T.*(100000./p).^(Rd./(cpd + cl.*mr));
                        % teB = exp( (lv0.*mr)./((cpd + mr*cl).*T));
                        % thetae = teA .* teB;
                        teA = T.*(100000./p).^(Rd./(cpd + cl.*mr));
                        teB = exp( (lv0.*mr)./((cpd + mr*cl).*T));
                        teC = (rh/100).^(  (-mr*Rv) ./ (cpd + cl.*mr)   );
                        thetae = teA .* teB .* teC;



                        %height indices of some desired geometric altitudes:
                        h100m  = 100 + ht(1); k100m =  find( abs(ht - h100m) == min(abs(ht - h100m)  ) ) ; k100m = k100m(1);
                        h1km   = 1000 + ht(1); k1km  =  find( abs(ht - h1km)  == min(abs(ht - h1km)  ) ) ; k1km = k1km(1);
                        h1500m = 1500 + ht(1); k1500m  =  find( abs(ht - h1500m)  == min(abs(ht - h1500m)  ) ) ; k1500m = k1500m(1);
                        h2km   = 2000 + ht(1); k2km  =  find( abs(ht - h2km)  == min(abs(ht - h2km)  ) ) ; k2km = k2km(1);
                        h3km   = 3000 + ht(1); k3km  =  find( abs(ht - h3km)  == min(abs(ht - h3km)  ) ) ; k3km = k3km(1);
                        h4km   = 4000 + ht(1); k4km  =  find( abs(ht - h4km)  == min(abs(ht - h4km)  ) ) ; k4km = k4km(1);
                        h5km   = 5000 + ht(1); k5km  =  find( abs(ht - h5km)  == min(abs(ht - h5km)  ) ) ; k5km = k5km(1);
                        h6km   = 6000 + ht(1); k6km  =  find( abs(ht - h6km)  == min(abs(ht - h6km)  ) ) ; k6km = k6km(1);
                        h9km   = 9000 + ht(1); k9km  =  find( abs(ht - h9km)  == min(abs(ht - h9km)  ) ) ; k9km = k9km(1);
                        h10km  = 10000 + ht(1); k10km  =  find( abs(ht - h10km)  == min(abs(ht - h10km)  ) ) ; k10km = k10km(1);
                        h12km  = 12000 + ht(1); k12km  =  find( abs(ht - h12km)  == min(abs(ht - h12km)  ) ) ; k12km = k12km(1);


                        %height indices of some desired geometric altitudes ASL:
                        h3kmasl   = 3000; k3kmasl  =  find( abs(htasl - h3kmasl)  == min(abs(htasl - h3kmasl)  ) ) ; k3kmasl = k3kmasl(1);
                        h4kmasl   = 4000; k4kmasl  =  find( abs(htasl - h4kmasl)  == min(abs(htasl - h4kmasl)  ) ) ; k4kmasl = k4kmasl(1);
                        h5kmasl   = 5000; k5kmasl  =  find( abs(htasl - h5kmasl)  == min(abs(htasl - h5kmasl)  ) ) ; k5kmasl = k5kmasl(1);
                        h6kmasl   = 6000; k6kmasl  =  find( abs(htasl - h6kmasl)  == min(abs(htasl - h6kmasl)  ) ) ; k6kmasl = k6kmasl(1);


                        %height indices of some desired pressure levels:
                        k925mb  =  find( abs(pres - 925.0) == min(abs(pres - 925.0) ) );  k925mb = k925mb(1);
                        k850mb  =  find( abs(pres - 850.0) == min(abs(pres - 850.0) ) );  k850mb = k850mb(1);
                        k700mb  =  find( abs(pres - 700.0) == min(abs(pres - 700.0) ) );  k700mb = k700mb(1);
                        k600mb  =  find( abs(pres - 600.0) == min(abs(pres - 600.0) ) );  k600mb = k600mb(1);
                        k500mb  =  find( abs(pres - 500.0) == min(abs(pres - 500.0) ) );  k500mb = k500mb(1);
                        k400mb  =  find( abs(pres - 400.0) == min(abs(pres - 400.0) ) );  k400mb = k400mb(1);
                        k300mb  =  find( abs(pres - 300.0) == min(abs(pres - 300.0) ) );  k300mb = k300mb(1);

                        hterr = 2600 - htasl(1); % height diff between AMF and ~ peak of terrain (ie ht agl of peak from amf's pov)
                        kterr =  find( abs(ht - hterr) == min(abs(ht - hterr)  ) ) ;    kterr = kterr(1);  %k index of terrain (agl)
                        %
                        %                 hterrp1 = 1000 + 2600 - htasl(1); % height diff between AMF and 1 km above peak of terrain (ie ht agl of peak from amf's pov)
                        %                 kterrp1 =  find( abs(ht - hterrp1) == min(abs(ht - hterrp1)  ) ) ;    kterrp1 = kterrp1(1);  %k index of terrain (agl) + 1 km
                        %
                        %                 hterrp1 = 3550;
                        %                 kterrp1 =  find( abs(ht - hterrp1) == min(abs(ht - hterrp1)  ) ) ;




                        %                     %%%%%%%%% find the adv_x,y corresponding to the current
                        %                     %%%%%%%%% cell basetime, using cell advection scene
                        %                     %%%%%%%%% (not unique cell motions per time)
                        %
                        %                     %note, wrfbasetimes fot cell t = 1-5 seem to not
                        %                     %necessarily be on discrete 15min intervals, so need to
                        %                     %account for that
                        %
                        %                     wbt = wrf_base_time(t,c) ;  %current cell basetime
                        %                     %wbt = wrf_base_time ;  %current cell basetime
                        %
                        %                     dtthresh = 449.99;  %seconds  - tolerance of wrf_basetime relative to the 15-min intervals
                        %
                        %                     %if cell's time is before 12utc:
                        %                     if( wbt < adv_basetime(1) - (dtthresh+1) )
                        %                         adv_x =  adv_cx(1);  %mean(U(1:k10km));
                        %                         adv_y =  adv_cy(1);  %mean(V(1:k10km));
                        %                     end
                        %
                        %                     %if cell's time is within dt intervals of 1200-2345 UTC
                        %                     for r = 1 : length(adv_basetime)
                        %                         if( wbt > adv_basetime(r) - dtthresh &  wbt < adv_basetime(r) + dtthresh )
                        %                             adv_x = adv_cx(r);
                        %                             adv_y = adv_cy(r);
                        %                         end
                        %                     end
                        %
                        %                     %now calc the cloud-rel wind:
                        %                     u_crel = U - adv_x ;
                        %                     v_crel = V - adv_y ;



                        %%%% calculate a unique cell-relative wind vector per time.
                        if( isnan(cell_motion_x)==0 &  isnan(cell_motion_y)==0  )
                            u_crel = U - cell_motion_x ;
                            v_crel = V - cell_motion_y ;
                        end



                        muCAPE = -999;
                        %muCIN = -999;
                        muLCL = -999;
                        muLFC = -999;
                        muEL = -999;
                        %                     muLCLp = -999;
                        %                     muLFCp = -999;
                        %                     muELp = -999;
                        %                    mu_negCINfract = -999;
                        mu_CIN = -999;
                        mu_CAPEACBL = -999;
                        %mu_CAPELCL = -999;
                        %                     eil50 = -999;
                        %                     eil25 = -999;
                        %                     eil10 = -999;
                        hMU = -999;
                        %teMU = -999;


                        %		    check = 5;

                        % [sfcCAPE, sfcCIN, sfcLCL, sfcLFC, sfcEL, sfcLCLp, sfcLFCp, sfcELp, sfc_negCINfract, sfc_negCINsum, sfc_CAPEACBL, sfc_CAPELCL, muCAPE, muCIN, muLCL, muLFC, muEL, muLCLp, muLFCp, muELp, mu_negCINfract, mu_negCINsum, mu_CAPEACBL, mu_CAPELCL, mlCAPE, mlCIN, mlLCL, mlLFC, mlEL, mlLCLp, mlLFCp, mlELp, ml_negCINfract, ml_negCINsum, ml_CAPEACBL, ml_CAPELCL, eil50, eil25, eil10, hMU, teSFC, teMU, teML] = SOUNDING_SFC_MU_ML_i_Adam_EIL(pres, tdry, TdC, ht);


            		    %disp('pre-SOUNDING ')
            		    %disp('cell = ')
            		    %disp(num2str(c))
            		    %disp(num2str(x))
            		    %disp(num2str(y))
            		    %disp(num2str(t))
            		    %disp('----')

            		    [sfcCAPE, sfcCIN, sfcLCL, sfcLFC, sfcEL, sfc_CAPEACBL, sfcCAPEpse, sfcCINpse, sfcLFCpse, sfcELpse, sfc_CAPEACBLpse,...
                            muCAPE,  muCIN,  muLCL, muLFC,  muEL,  mu_CAPEACBL,  muCAPEpse,  muCINpse,  muLFCpse,  muELpse,  mu_CAPEACBLpse, ...
                            mlCAPE,  mlCIN,  mlLCL, mlLFC,  mlEL,  ml_CAPEACBL,  mlCAPEpse,  mlCINpse,  mlLFCpse,  mlELpse,  ml_CAPEACBLpse,...
                            hMU, teSFC, teMU, teML, teSFCpse, teMUpse, teMLpse] = SOUNDING_SFC_MU_ML_WIDlasso( pres, tdry, TdC, ht);



            		    %disp('post-SOUNDING ')
                        %disp('cell = ')
            		    %disp(num2str(c))
                        %disp(num2str(x))
                        %disp(num2str(y))
                        %disp(num2str(t))
                        %disp('----')

                        pCAPE_sfc(x,y,t,c)   = sfcCAPEpse;
                        pCAPE_ml(x,y,t,c)    = mlCAPEpse;
                        pCAPE_mu(x,y,t,c)    = muCAPEpse;


                        CAPE_mu(x,y,t,c) = muCAPE;
                        %CIN_IB_mu(x,y,t) = muCIN;
                        LCL_height_mu(x,y,t,c) = muLCL;
                        LFC_height_mu(x,y,t,c) = muLFC;
                        EL_height_mu(x,y,t,c) = muEL;
                        %LCL_pres_mu(x,y,t) = muLCLp;
                        %LFC_pres_mu(x,y,t) = muLFCp;
                        %EL_pres_mu(x,y,t) = muELp;
                        %CIN_NA_fract_mu(x,y,t) = mu_negCINfract;
                        CIN_mu(x,y,t,c) = muCIN;
                        CAPEacbl_mu(x,y,t,c) = mu_CAPEACBL;
                        %CAPElcl_IB_mu(x,y,t) = mu_CAPELCL;
                        initial_ht_parcel_mu(x,y,t,c) = hMU;
                        %tallenough_mu(x,y,t) = teMU;


                        CAPE_sfc(x,y,t,c) = sfcCAPE;
                        %CIN_IB_sfc(x,y,t) = sfcCIN;
                        LCL_height_sfc(x,y,t,c) = sfcLCL;
                        LFC_height_sfc(x,y,t,c) = sfcLFC;
                        EL_height_sfc(x,y,t,c) = sfcEL;
                        %LCL_pres_sfc(x,y,t) = sfcLCLp;
                        %LFC_pres_sfc(x,y,t) = sfcLFCp;
                        %EL_pres_sfc(x,y,t) = sfcELp;
                        %CIN_NA_fract_sfc(x,y,t) = sfc_negCINfract;
                        CIN_sfc(x,y,t,c) = sfcCIN;
                        CAPEacbl_sfc(x,y,t,c) = sfc_CAPEACBL;
                        %CAPElcl_IB_sfc(x,y,t) = sfc_CAPELCL;
                        %tallenough_sfc(x,y,t) = teSFC;

                        CAPE_ml(x,y,t,c) = mlCAPE;
                        %CIN_IB_ml(x,y,t) = mlCIN;
                        LCL_height_ml(x,y,t,c) = mlLCL;
                        LFC_height_ml(x,y,t,c) = mlLFC;
                        EL_height_ml(x,y,t,c) = mlEL;
                        %LCL_pres_ml(x,y,t) = mlLCLp;
                        %LFC_pres_ml(x,y,t) = mlLFCp;
                        %EL_pres_ml(x,y,t) = mlELp;
                        %CIN_NA_fract_ml(x,y,t) = ml_negCINfract;
                        CIN_ml(x,y,t,c) = mlCIN;
                        CAPEacbl_ml(x,y,t,c) = ml_CAPEACBL;
                        %CAPElcl_IB_ml(x,y,t) = ml_CAPELCL;
                        %tallenough_ml(x,y,t) = teML;


                        tallenough_mu(x,y,t,c)    = teMU;

                        pCIN_mu(x,y,t,c)       = muCINpse;
                        pCAPEacbl_mu(x,y,t,c)     = mu_CAPEACBLpse;
                        pLFC_height_mu(x,y,t,c)   = muLFCpse;
                        pEL_height_mu(x,y,t,c)    = muELpse;
                        ptallenough_mu(x,y,t,c)   = teMUpse;

                        tallenough_sfc(x,y,t,c)    = teSFC;

                        pCIN_sfc(x,y,t,c)       = sfcCINpse;
                        pCAPEacbl_sfc(x,y,t,c)     = sfc_CAPEACBLpse;
                        pLFC_height_sfc(x,y,t,c)   = sfcLFCpse;
                        pEL_height_sfc(x,y,t,c)    = sfcELpse;
                        ptallenough_sfc(x,y,t,c)   = teSFCpse;

                        tallenough_ml(x,y,t,c)    = teML;

                        pCIN_ml(x,y,t,c)       = mlCINpse;
                        pCAPEacbl_ml(x,y,t,c)     = ml_CAPEACBLpse;
                        pLFC_height_ml(x,y,t,c)   = mlLFCpse;
                        pEL_height_ml(x,y,t,c)    = mlELpse;
                        ptallenough_ml(x,y,t,c)   = teMLpse;

                        %EIL10_top_height(x,y,t) = eil10(2);
                        %EIL10_bot_height(x,y,t) = eil10(1);
                        %EIL25_top_height(x,y,t) = eil25(2);
                        %EIL25_bot_height(x,y,t) = eil25(1);
                        %EIL50_top_height(x,y,t) = eil50(2);
                        %EIL50_bot_height(x,y,t) = eil50(1);



        	            % record values of var-height temperatures:


            		    % LCL temperatures
                        if(   isempty(sfcLCL)==0 & isnan(sfcLCL)==0   )
                            k_vari = find( abs(ht - (sfcLCL)) == min(abs(ht - (sfcLCL))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            LCL_temp_sfc(x,y,t,c) = T(k_vari);
                        else
                            LCL_temp_sfc(x,y,t,c) = NaN;
                        end

                        if(   isempty(muLCL)==0 & isnan(muLCL)==0   )
                            k_vari = find( abs(ht - (muLCL)) == min(abs(ht - (muLCL))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            LCL_temp_mu(x,y,t,c) = T(k_vari);
                        else
                            LCL_temp_mu(x,y,t,c) = NaN;
                        end


                        if(   isempty(mlLCL)==0 & isnan(mlLCL)==0   )
                            k_vari = find( abs(ht - (mlLCL)) == min(abs(ht - (mlLCL))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            LCL_temp_ml(x,y,t,c) = T(k_vari);
                        else
                            LCL_temp_ml(x,y,t,c) = NaN;
                        end

            		    %revEL temperatures

                        if(   isempty(sfcEL)==0 & isnan(sfcEL)==0   )
                            k_vari = find( abs(ht - (sfcEL)) == min(abs(ht - (sfcEL))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            EL_temp_sfc(x,y,t,c) = T(k_vari);
                        else
                            EL_temp_sfc(x,y,t,c) = NaN;
                        end


                        if(   isempty(muEL)==0 & isnan(muEL)==0   )
                            k_vari = find( abs(ht - (muEL)) == min(abs(ht - (muEL))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            EL_temp_mu(x,y,t,c) = T(k_vari);
                        else
                            EL_temp_mu(x,y,t,c) = NaN;
                        end


                        if(   isempty(mlEL)==0 & isnan(mlEL)==0   )
                            k_vari = find( abs(ht - (mlEL)) == min(abs(ht - (mlEL))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            EL_temp_ml(x,y,t,c) = T(k_vari);
                        else
                            EL_temp_ml(x,y,t,c) = NaN;
                        end


                        %pseudoadiab EL temperature

                        if(   isempty(sfcELpse)==0 & isnan(sfcELpse)==0   )
                            k_vari = find( abs(ht - (sfcELpse)) == min(abs(ht - (sfcELpse))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            pEL_temp_sfc(x,y,t,c) = T(k_vari);
                        else
                            pEL_temp_sfc(x,y,t,c) = NaN;
                        end


                        if(   isempty(muELpse)==0 & isnan(muELpse)==0   )
                            k_vari = find( abs(ht - (muELpse)) == min(abs(ht - (muELpse))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            pEL_temp_mu(x,y,t,c) = T(k_vari);
                        else
                            pEL_temp_mu(x,y,t,c) = NaN;
                        end


                        if(   isempty(mlELpse)==0 & isnan(mlELpse)==0   )
                            k_vari = find( abs(ht - (mlELpse)) == min(abs(ht - (mlELpse))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            pEL_temp_ml(x,y,t,c) = T(k_vari);
                        else
                            pEL_temp_ml(x,y,t,c) = NaN;
                        end


            		    %rev LFC temperature
                        if(   isempty(sfcLFC)==0 & isnan(sfcLFC)==0   )
                            k_vari = find( abs(ht - (sfcLFC)) == min(abs(ht - (sfcLFC))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            LFC_temp_sfc(x,y,t,c) = T(k_vari);
                        else
                            LFC_temp_sfc(x,y,t,c) = NaN;
                        end

                        if(   isempty(muLFC)==0 & isnan(muLFC)==0   )
                            k_vari = find( abs(ht - (muLFC)) == min(abs(ht - (muLFC))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            LFC_temp_mu(x,y,t,c) = T(k_vari);
                        else
                            LFC_temp_mu(x,y,t,c) = NaN;
                        end


                        if(   isempty(mlLFC)==0 & isnan(mlLFC)==0   )
                            k_vari = find( abs(ht - (mlLFC)) == min(abs(ht - (mlLFC))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            LFC_temp_ml(x,y,t,c) = T(k_vari);
                        else
                            LFC_temp_ml(x,y,t,c) = NaN;
                        end


                        %pseudoadiab LFC temperature
                        if(   isempty(sfcLFCpse)==0 & isnan(sfcLFCpse)==0   )
                            k_vari = find( abs(ht - (sfcLFCpse)) == min(abs(ht - (sfcLFCpse))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            pLFC_temp_sfc(x,y,t,c) = T(k_vari);
                        else
                            pLFC_temp_sfc(x,y,t,c) = NaN;
                        end

                        if(   isempty(muLFCpse)==0 & isnan(muLFCpse)==0   )
                            k_vari = find( abs(ht - (muLFCpse)) == min(abs(ht - (muLFCpse))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            pLFC_temp_mu(x,y,t,c) = T(k_vari);
                        else
                            pLFC_temp_mu(x,y,t,c) = NaN;
                        end


                        if(   isempty(mlLFCpse)==0 & isnan(mlLFCpse)==0   )
                            k_vari = find( abs(ht - (mlLFCpse)) == min(abs(ht - (mlLFCpse))  ) ) ;  k_vari = k_vari(1);   %ht index vari
                            pLFC_temp_ml(x,y,t,c) = T(k_vari);
                        else
                            pLFC_temp_ml(x,y,t,c) = NaN;
                        end




                        %                   check = 6;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % variables that dont care about lifting parcels:
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                        %%%%%%%%%%%%%%%%%%%%%%%
                        % precipitable water
                        %%%%%%%%%%%%%%%%%%%%%%%


                        chunks = mr; chunks(:) = 0;
                        pp = pres.*100 ;
                        for k = 1:length(mr)-1
                            chunks(k) = 0.5*( mr(k+1) + mr(k) ) .* ( pp(k+1)-pp(k) ) ;
                        end
                        %   PW_interpsonde = vertcat( PW_interpsonde, -(1/9.8).*sum(chunks,'omitnan') ) ; %m - spec hum
                        PW(x,y,t,c) = -(100/9.8).*(1/997).*sum(chunks,'omitnan') ; %the 100 is to make it cm



                        %%%%%%%%%%%%%%%%%%%%%%%
                        %freezing level:
                        %%%%%%%%%%%%%%%%%%%%%%%

                        frz = find( abs( tdry(:)) == min( abs(tdry(:) ) ) )  ; frz = frz(1);
                        freezing_ht(x,y,t,c) = ht( frz ) ;         % m AGL


                        %%%%%%%%%%%%%%%%%%%%%%%
                        %upslope variables between the ground and peak terrain altitude:
                        %%%%%%%%%%%%%%%%%%%%%%%
                       
                        
                        NegU = find( U(1:kterr) < 0 );   %data levels below peak terrain with upslope flow component
                        
                       
                        if(isempty(NegU) == 0) %there are easterly low-level winds
                            upslope_max(x,y,t,c) = min(U(NegU)); %max upsloping wind below terrain top
                        elseif(isempty(NegU) == 1) %there are not easterly low-level winds
                            upslope_max(x,y,t,c) = NaN ; %max upsloping wind below terrain top
                        end
                       
                        upslope_mean(x,y,t,c) = LAYERMEAN(U(:),ht,1,kterr); %mean upsloping wind b/w AMF sfc and terrain top (neg = upslope)
                        upslope_fract(x,y,t,c) = length(NegU)/kterr ;  %fraction of vertical data points between AMF surface and the terrain peak with upslope flow
                      
                        %attempt to find depth of upslope flow by locating first westerly U
                        %above the ground. This will probably not work out if there are micro-layers of
                        %non-upslope flow:
                      
                        PosU = find( U(1:k5km) > 0 );         %jnm_new_sep8  - allowing the upslope detection to go deeper
                        %PosU = find( u(1:kterr,t) > 0 );       %jnm_new_sep8
                        if (isempty(PosU)==0)   %upslope stops before terrain ht
                           upslope_depth(x,y,t,c) = ht(PosU(1));   %first height between AMF sfc and terrain peak at which U no longer upslopes
                        elseif (isempty(PosU)==1)  % all upslope
                           upslope_depth(x,y,t,c) = hterr;
                        end



                        %         u_mean_0to1 = vertcat(u_mean_0to1,LAYERMEAN(u(:,t),ht,1,k1km));
                        %         v_mean_0to1 = vertcat(v_mean_0to1,LAYERMEAN(v(:,t),ht,1,k1km));
                        %         u_mean_5to6 = vertcat(u_mean_5to6,LAYERMEAN(u(:,t),ht,k5km,k6km ));
                        %         v_mean_5to6 = vertcat(v_mean_5to6,LAYERMEAN(v(:,t),ht,k5km,k6km ));

                        %         V_mean_terr = vertcat(V_mean_terr,LAYERMEAN(v(:,t),ht,kterr,kterrp1));
                        %         U_mean_terr = vertcat(U_mean_terr,LAYERMEAN(u(:,t),ht,kterr,kterrp1));

                        %         Ucrel_mean_0to4 = vertcat(Ucrel_mean_0to4, LAYERMEAN(u_crel(:),ht,1,k4km ) );
                        %         Vcrel_mean_0to4 = vertcat(Vcrel_mean_0to4, LAYERMEAN(v_crel(:),ht,1,k4km ) );


                        %%%%%%%%%%%%%%%%%%%%%%%
                        % boundary layer depth
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         % calculate lapse rate of Theta
                        %         LRth = theta(:,t); LRth(:) = NaN;
                        %         for k = 2:length(LRth)-1
                        %             LRth(k) = (theta(k+1,t) - theta(k-1,t))./(ht(k+1) - ht(k-1))*1000  ;  %K/km
                        %         end
                        %         lapse_rate_theta = horzcat(lapse_rate_theta,LRth);   %full LR array
                        %         LRth = movmean(LRth,9);  %smoothed with 9-pt pass
                        %         smoothed_lapse_rate_theta = horzcat(smoothed_lapse_rate_theta,LRth);   %full LR array


                        %         % call BL (Liu-Liang) function here:
                        %         clear BLtype BLtop_LL
                        %         [BLtype, BLtop_LL] = BoundaryLayerDepth( pres(:,t),theta(:,t),u(:,t),v(:,t),ht ) ;
                        %         SBLorCBL = vertcat( SBLorCBL, BLtype );
                        %         height_BLtop_LiuLiang = vertcat(height_BLtop_LiuLiang,BLtop_LL) ;

                        %         %find k of BLtop on native vertical grid
                        %         if(isnan(BLtop_LL) == 0)
                        %             kBLtop  =  find( abs(ht - BLtop_LL) == min(abs(ht - BLtop_LL)  ) ) ;
                        %             kBLtop = kBLtop(1);
                        %             presBLtop = vertcat(presBLtop,pres(kBLtop));
                        %         elseif(isnan(BLtop_LL))
                        %             presBLtop = vertcat(presBLtop,NaN);
                        %         end



                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%  humidity metrics
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         if(isnan(BLtop_LL) == 0)
                        %             BLmean_RH = vertcat( BLmean_RH, LAYERMEAN(rh(:,t),ht,1,kBLtop) );
                        %             BLmean_DewptDepr = vertcat( BLmean_DewptDepr, LAYERMEAN(TdDepress(:,t),ht,1,kBLtop) );  %deg C
                        %             BLmean_vapmixrat = vertcat( BLmean_vapmixrat, LAYERMEAN(mr(:,t)*1000,ht,1,kBLtop) );   %g/kg
                        %         elseif(isnan(BLtop_LL))
                        %             BLmean_RH = vertcat(BLmean_RH, NaN);
                        %             BLmean_DewptDepr = vertcat( BLmean_DewptDepr,NaN);
                        %             BLmean_vapmixrat = vertcat( BLmean_vapmixrat,NaN);
                        %         end


                        %         % calculate lapse rate of vapor mixing ratio
                        %         MRgrad = theta(:,t); MRgrad(:) = NaN;
                        %         for k = 2:length(MRgrad)-1
                        %             MRgrad(k) = (mr(k+1,t) - mr(k-1,t))./(ht(k+1) - ht(k-1))*1000*1000  ;  %g/kg /km
                        %         end
                        %         MRgrad = movmean(MRgrad,9);  %smoothed with 9-pt pass
                        %
                        %         if(isnan(BLtop_LL) == 0)
                        %             BLmean_vapmixratgrad = vertcat(BLmean_vapmixratgrad, LAYERMEAN( MRgrad(:),ht,1,kBLtop) ) ;
                        %         elseif(isnan(BLtop_LL))
                        %             BLmean_vapmixratgrad = vertcat(BLmean_vapmixratgrad, NaN ) ;
                        %         end


                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%  lapse rate metrics
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         % calculate lapse rate of Tv
                        %         LRtv = Tv(:,t); LRtv(:) = NaN;
                        %         for k = 2:length(LRtv)-1
                        %             LRtv(k) = (Tv(k+1,t) - Tv(k-1,t))./(ht(k+1) - ht(k-1))*1000  ;  %K/km
                        %         end
                        %         LRtv = movmean(LRtv,9);  %smoothed with 9-pt pass
                        %         lapse_rate_tv = horzcat(lapse_rate_tv,LRtv);  %full-sized LRtv array


                        %         if(isnan(BLtop_LL) == 0)
                        %             BLmean_LapseRate = vertcat( BLmean_LapseRate, LAYERMEAN( LRtv(:),ht,k100m,kBLtop) );
                        %
                        %             %height = PBL_top + XX km
                        %             kBLtopp05  =  find( abs(ht - (BLtop_LL + 500)) == min(abs(ht - (BLtop_LL + 500)  ))) ;
                        %             kBLtopp1  =  find( abs(ht - (BLtop_LL + 1000)) == min(abs(ht - (BLtop_LL + 1000)  ))) ;
                        %             kBLtopp2  =  find( abs(ht - (BLtop_LL + 2000)) == min(abs(ht - (BLtop_LL + 2000)  ) ) );
                        %
                        %             LapseRate_Tv_mean_BLinver = vertcat( LapseRate_Tv_mean_BLinver, LAYERMEAN( LRtv(:),ht,kBLtop,kBLtopp05(1)) ) ;
                        %
                        %             %diagnostic:
                        %             LRdiag = horzcat(sonday,t,LAYERMEAN( LRtv(:),ht,kBLtop,kBLtopp05(1)));
                        %
                        %             LapseRate_Tv_max_BLinver = vertcat( LapseRate_Tv_max_BLinver, max( LRtv(kBLtop:kBLtopp1)) ) ;
                        %
                        %             LRTVthresh = -4.0;
                        %             inv = find(LRtv(kBLtop:kBLtopp2) >= LRTVthresh );
                        %             if( isempty(inv) == 0  )
                        %                 LapseRate_Tv_depth_BLinver = vertcat( LapseRate_Tv_depth_BLinver, ht(inv(end)) - ht(inv(1)) );
                        %             else
                        %                 LapseRate_Tv_depth_BLinver = vertcat( LapseRate_Tv_depth_BLinver, 0.0 );
                        %             end
                        %
                        %         elseif(isnan(BLtop_LL))
                        %
                        %             BLmean_LapseRate = vertcat( BLmean_LapseRate, NaN );
                        %             LapseRate_Tv_mean_BLinver = vertcat( LapseRate_Tv_mean_BLinver, NaN );
                        %             LapseRate_Tv_max_BLinver = vertcat( LapseRate_Tv_max_BLinver, NaN );
                        %             LapseRate_Tv_depth_BLinver = vertcat( LapseRate_Tv_depth_BLinver, NaN );
                        %         end




                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%  surface metrics
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         clear u10 v10 rh10
                        u10 = interp1(ht(1:5),U(1:5),10);
                        v10 = interp1(ht(1:5),V(1:5),10);
                        %         rh10 = interp1(ht(1:5),rh(1:5,t),10);
                        %         U_10m = vertcat(U_10m, u10 );
                        %         V_10m = vertcat(V_10m, v10 );
                        %         rh_10m = vertcat(rh_10m, rh10 );





                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%  shear metrics
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         % calculate lapse rate of Tv
                        %         ShearU = u(:,t); ShearU(:) = NaN;
                        %         ShearV = u(:,t); ShearV(:) = NaN;
                        %         for k = 2:length(ShearU)-1
                        %             ShearU(k) = (u(k+1,t) - u(k-1,t))./(ht(k+1) - ht(k-1))*1000  ;  %m/s /km
                        %             ShearV(k) = (v(k+1,t) - v(k-1,t))./(ht(k+1) - ht(k-1))*1000  ;  %m/s /km
                        %         end
                        %         ShearU = movmean(ShearU,9);  %smoothed with 9-pt pass
                        %         ShearV = movmean(ShearV,9);  %smoothed with 9-pt pass
                        %         magshear = ( ShearU.*ShearU + ShearV.*ShearV ).^0.5;
                        %         Mag_Shear = horzcat(Mag_Shear,magshear); %full-sized array of shear magnitude


                        %         ksfcLCL  =  find( abs(ht - sfcLCL) == min(abs(ht - sfcLCL)  ) ) ;
                        %         shearmag_mean_subcloud_sfc = vertcat( shearmag_mean_subcloud_sfc, LAYERMEAN( magshear(:),ht,1,ksfcLCL ) );



                        du = U(k1km) - u10; dv = V(k1km) - v10;
                        shear_mag_bulk_0to1km(x,y,t,c) = ( du*du + dv*dv ).^0.5 ;

                        sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        if( sdir > 360.0 )
                            sdir = sdir - 360 ;
                        end
                        shear_dir_bulk_0to1km(x,y,t,c) =  sdir;




                        du = U(k3km) - u10; dv = V(k3km) - v10;
                        shear_mag_bulk_0to3km(x,y,t,c) =  ( du*du + dv*dv ).^0.5 ;

                        sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        if( sdir > 360.0 )
                            sdir = sdir - 360 ;
                        end
                        shear_dir_bulk_0to3km(x,y,t,c) =  sdir;




                        du = U(k6km) - u10; dv = V(k6km) - v10;
                        shear_mag_bulk_0to6km(x,y,t,c) =  ( du*du + dv*dv ).^0.5 ;

                        sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        if( sdir > 360.0 )
                            sdir = sdir - 360 ;
                        end
                        shear_dir_bulk_0to6km(x,y,t,c) = sdir;



                        du = U(k9km) - u10; dv = V(k9km) - v10;
                        shear_mag_bulk_0to9km(x,y,t,c) =  ( du*du + dv*dv ).^0.5 ;

                        sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        if( sdir > 360.0 )
                            sdir = sdir - 360 ;
                        end
                        shear_dir_bulk_0to9km(x,y,t,c) = sdir;



                        %         clear du dv sdir
                        %         du = u(k12km,t) - u10; dv = v(k12km,t) - v10;
                        %         shear_mag_bulk_0to1km2 = vertcat( shear_mag_bulk_0to1km2, ( du*du + dv*dv ).^0.5 );
                        %
                        %         sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        %         if( sdir > 360.0 )
                        %             sdir = sdir - 360 ;
                        %         end
                        %         shear_dir_bulk_0to12km = vertcat(shear_dir_bulk_0to12km, sdir);
                        %         %end new 22 Nov 2021
                        %
                        %
                        %         %new 13 Dec
                        %         clear du dv sdir
                        %         du = u(kterrp1,t) - u(kterr,t); dv = v(kterrp1,t) - v(kterr,t);
                        %         shearmag_bulk_terr = vertcat( shearmag_bulk_terr, ( du*du + dv*dv ).^0.5 );
                        %
                        %         sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        %         if( sdir > 360.0 )
                        %             sdir = sdir - 360 ;
                        %         end
                        %         shear_dir_bulk_terr = vertcat(shear_dir_bulk_terr, sdir);
                        %         %end new 13 Dec 2021




                        %         if(isnan(BLtop_LL) == 0)
                        %             clear du dv sdir
                        %             du = u(kBLtop,t) - u10; dv = v(kBLtop,t) - v10;
                        %             shearmag_bulk_BL = vertcat( shearmag_bulk_BL, ( du*du + dv*dv ).^0.5 );
                        %             shearmag_mean_BL  = vertcat( shearmag_mean_BL, LAYERMEAN( magshear(:),ht,1,kBLtop ) );
                        %
                        %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                        %             if( sdir > 360.0 )
                        %                 sdir = sdir - 360 ;
                        %             end
                        %             shear_dir_bulk_BL = vertcat(shear_dir_bulk_BL, sdir);
                        %
                        %         elseif(isnan(BLtop_LL))
                        %             shearmag_bulk_BL = vertcat( shearmag_bulk_BL, NaN );
                        %             shearmag_mean_BL  = vertcat( shearmag_mean_BL, NaN );
                        %             shear_dir_bulk_BL = vertcat(shear_dir_bulk_BL, NaN);
                        %         end
                        %
                        %
                        %         shearmag_max  = vertcat( shearmag_max, max( magshear) );  %max shear throughout entire sounding
                        %         shearmag_max_height  = vertcat( shearmag_max_height,  ht(  find(  magshear == max( magshear) ) ) ); %height of max shear throughout entire sounding




                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%  etc. metrics
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         if(isnan(BLtop_LL) == 0)
                        %             BL_mean_theta = vertcat(BL_mean_theta, LAYERMEAN( theta(:,t),ht,1,kBLtop) ) ;
                        %             thetae_mean_BL = vertcat(thetae_mean_BL, LAYERMEAN( thetae(:,t),ht,1,kBLtop) ) ;  %jnm_new_23mar2022
                        %         elseif(isnan(BLtop_LL))
                        %             BL_mean_theta = vertcat(BL_mean_theta, NaN ) ;
                        %             thetae_mean_BL = vertcat(thetae_mean_BL, NaN) ; %jnm_new_23mar2022
                        %         end


                        %%%%%%%%%%%%%%%%%%%%%%%%
                        %%% some standard height or pressure level metrics
                        %%%%%%%%%%%%%%%%%%%%%%%%

                        %         U_1500m = vertcat(U_1500m,u(k1500m,t));
                        %         V_1500m = vertcat(V_1500m,v(k1500m,t));
                        %         U_3000m = vertcat(U_3000m,u(k3km,t));
                        %         V_3000m = vertcat(V_3000m,v(k3km,t));
                        %                              U_4000m(x,y,t) = U(k4km);
                        %                              V_4000m(x,y,t) = V(k4km);
                        %         U_6000m = vertcat(U_6000m,u(k6km,t));
                        %         V_6000m = vertcat(V_6000m,v(k6km,t));

                        %         rvap_1500m = vertcat(rvap_1500m,mr(k1500m,t)*1000);
                        %         rvap_3000m = vertcat(rvap_3000m,mr(k3km,t)*1000);
                        %         rvap_6000m = vertcat(rvap_6000m,mr(k6km,t)*1000);
                        %         rh_1500m = vertcat(rh_1500m,rh(k1500m,t));
                        %         rh_3000m = vertcat(rh_3000m,rh(k3km,t));
                        %         rh_6000m = vertcat(rh_6000m,rh(k6km,t));

                        %         U_850mb = vertcat(U_850mb,u(k850mb,t));
                        %         V_850mb = vertcat(V_850mb,v(k850mb,t));
                        %         U_700mb = vertcat(U_700mb,u(k700mb,t));
                        %         V_700mb = vertcat(V_700mb,v(k700mb,t));
                        %         U_500mb = vertcat(U_500mb,u(k500mb,t));
                        %         V_500mb = vertcat(V_500mb,v(k500mb,t));
                        %         U_300mb = vertcat(U_300mb,u(k300mb,t));
                        %         V_300mb = vertcat(V_300mb,v(k300mb,t));


                        rvap_925mb(x,y,t,c) = mr(k925mb)*1000;
                        rvap_850mb(x,y,t,c) = mr(k850mb)*1000;
                        rvap_700mb(x,y,t,c) = mr(k700mb)*1000;
                        rvap_600mb(x,y,t,c) = mr(k600mb)*1000;
                        rvap_500mb(x,y,t,c) = mr(k500mb)*1000;
                        rvap_400mb(x,y,t,c) = mr(k400mb)*1000;
                        rvap_300mb(x,y,t,c) = mr(k300mb)*1000;

                        rh_925mb(x,y,t,c) = rh(k925mb);
                        rh_850mb(x,y,t,c) = rh(k850mb);
                        rh_600mb(x,y,t,c) = rh(k600mb);
                        rh_700mb(x,y,t,c) = rh(k700mb);
                        rh_500mb(x,y,t,c) = rh(k500mb);
                        rh_400mb(x,y,t,c) = rh(k400mb);
                        rh_300mb(x,y,t,c) = rh(k300mb);

                        rvap_3000masl(x,y,t,c) = mr(k3kmasl)*1000;
                        rvap_4000masl(x,y,t,c) = mr(k4kmasl)*1000;
                        rvap_5000masl(x,y,t,c) = mr(k5kmasl)*1000;
                        rvap_6000masl(x,y,t,c) = mr(k6kmasl)*1000;

                        rh_3000masl(x,y,t,c) = rh(k3kmasl);
                        rh_4000masl(x,y,t,c) = rh(k4kmasl);
                        rh_5000masl(x,y,t,c) = rh(k5kmasl);
                        rh_6000masl(x,y,t,c) = rh(k6kmasl);



            		    if( isempty(muLCL)==0 & isnan(muLCL)==0 )
                            kLCLmu1500 = find( abs(ht - (muLCL+1500)) == min(abs(ht - (muLCL+1500))  ) ) ;  kLCLmu1500 = kLCLmu1500(1);   %ht index of z = 1.5km + LCL
                            rh_at_muLCLplus1500m(x,y,t,c) = rh(kLCLmu1500);
            		    end


                        %adam pseudoadiab vars
                        if( isnan(muLFCpse) == 0 )

                            kLFCmupse   =  find( abs(ht - muLFCpse) == min(abs(ht - muLFCpse)  ) ) ;                  kLFCmupse = kLFCmupse(1);
                            kACBLmupse  =  find( abs(ht - (muLFCpse+1500)) == min(abs(ht - (muLFCpse+1500))  ) ) ;   kACBLmupse = kACBLmupse(1);

                            rh_at_pmuLFCplus1500m(x,y,t,c) = rh(kACBLmupse);

                            RH_mean_ACBL_pmu(x,y,t,c) =       LAYERMEAN( rh(:),    ht,kLFCmupse,kACBLmupse )  ;
                            thetae_mean_ACBL_pmu(x,y,t,c) =   LAYERMEAN( thetae(:),ht,kLFCmupse,kACBLmupse )  ;
                            Ucrel_mean_ACBL_pmu(x,y,t,c) =    LAYERMEAN( u_crel(:),ht,kLFCmupse,kACBLmupse )  ;
                            Vcrel_mean_ACBL_pmu(x,y,t,c) =    LAYERMEAN( v_crel(:),ht,kLFCmupse,kACBLmupse )  ;

                            du = U(kACBLmupse,t) - U(kLFCmupse,t); dv = V(kACBLmupse,t) - V(kLFCmupse,t);
                            shear_mag_bulk_ACBL_pmu(x,y,t,c) = ( du*du + dv*dv ).^0.5 ;

                        else

                            rh_at_pmuLFCplus1500m(x,y,t,c) = NaN;
                            RH_mean_ACBL_pmu(x,y,t,c) = NaN;
                            shear_mag_bulk_ACBL_pmu(x,y,t,c) = NaN;
                            thetae_mean_ACBL_pmu(x,y,t,c) = NaN;
                            Ucrel_mean_ACBL_pmu(x,y,t,c) = NaN;
                            Vcrel_mean_ACBL_pmu(x,y,t,c) = NaN;

                        end

                        %%%%%%%%%%%%%%%%%%%%%%%
                        %%%  EIL stuff
                        %%%%%%%%%%%%%%%%%%%%%%%

                        %         if(isnan(eil250(1)) == 0 & isnan(eil250(2)) == 0)
                        %             kBotEil250  =  find( abs( ht - eil250(1) ) == min(abs( ht - eil250(1) )  ) ) ; kBotEil250 = kBotEil250(1);
                        %             kTopEil250  =  find( abs( ht - eil250(2) ) == min(abs( ht - eil250(2) )  ) ) ; kTopEil250 = kTopEil250(1);
                        %             EIL250_mean_U = vertcat(EIL250_mean_U,  LAYERMEAN(u(:,t),ht,kBotEil250,kTopEil250) );
                        %             EIL250_mean_V = vertcat(EIL250_mean_V,  LAYERMEAN(v(:,t),ht,kBotEil250,kTopEil250) );
                        %             EIL250_mean_VapMixRat = vertcat(EIL250_mean_VapMixRat,  LAYERMEAN(mr(:,t)*1000,ht,kBotEil250,kTopEil250) );
                        %             EIL250_mean_RH = vertcat(EIL250_mean_RH, LAYERMEAN(rh(:,t),ht,kBotEil250,kTopEil250 ) );
                        %         else
                        %             EIL250_mean_U = vertcat(EIL250_mean_U,  NaN );
                        %             EIL250_mean_V = vertcat(EIL250_mean_V,  NaN );
                        %             EIL250_mean_VapMixRat = vertcat(EIL250_mean_VapMixRat,  NaN );
                        %             EIL250_mean_RH = vertcat(EIL250_mean_RH, NaN );
                        %         end
                        %
                        %         if(isempty(adv_ind) == 0 & noadv == 0 & isnan(eil25(1)) == 0 & isnan(eil250(2)) == 0)
                        %             Ucrel_mean_EIL250 = vertcat(Ucrel_mean_EIL250, LAYERMEAN( u_crel(:),ht,kBotEil250,kTopEil250 ) ) ;
                        %             Vcrel_mean_EIL250 = vertcat(Vcrel_mean_EIL250, LAYERMEAN( v_crel(:),ht,kBotEil250,kTopEil250 ) ) ;
                        %         else
                        %             Ucrel_mean_EIL250 = vertcat(Ucrel_mean_EIL250, NaN ) ;
                        %             Vcrel_mean_EIL250 = vertcat(Vcrel_mean_EIL250, NaN ) ;
                        %         end
                        %
                        %
                        %
                        %
                        %         if(isnan(eil100(1)) == 0 & isnan(eil100(2)) == 0)
                        %             kBotEil100  =  find( abs( ht - eil100(1) ) == min(abs( ht - eil100(1) )  ) ) ; kBotEil100 = kBotEil100(1);
                        %             kTopEil100  =  find( abs( ht - eil100(2) ) == min(abs( ht - eil100(2) )  ) ) ; kTopEil100 = kTopEil100(1);
                        %             EIL100_mean_U = vertcat(EIL100_mean_U,  LAYERMEAN(u(:,t),ht,kBotEil100,kTopEil100) );
                        %             EIL100_mean_V = vertcat(EIL100_mean_V,  LAYERMEAN(v(:,t),ht,kBotEil100,kTopEil100) );
                        %             EIL100_mean_VapMixRat = vertcat(EIL100_mean_VapMixRat,  LAYERMEAN(mr(:,t)*1000,ht,kBotEil100,kTopEil100) );
                        %             EIL100_mean_RH = vertcat(EIL100_mean_RH, LAYERMEAN(rh(:,t)*1000,ht,kBotEil100,kTopEil100) );
                        %         else
                        %             EIL100_mean_U = vertcat(EIL100_mean_U,  NaN );
                        %             EIL100_mean_V = vertcat(EIL100_mean_V,  NaN );
                        %             EIL100_mean_VapMixRat = vertcat(EIL100_mean_VapMixRat,  NaN );
                        %             EIL100_mean_RH = vertcat(EIL100_mean_RH, NaN );
                        %         end
                        %
                        %         if(isempty(adv_ind) == 0 & noadv == 0 & isnan(eil10(1)) == 0 & isnan(eil10(2)) == 0)
                        %             Ucrel_mean_EIL100 = vertcat(Ucrel_mean_EIL100, LAYERMEAN( u_crel(:),ht,kBotEil100,kTopEil100) ) ;
                        %             Vcrel_mean_EIL100 = vertcat(Vcrel_mean_EIL100, LAYERMEAN( v_crel(:),ht,kBotEil100,kTopEil100) ) ;
                        %         else
                        %             Ucrel_mean_EIL100 = vertcat(Ucrel_mean_EIL100, NaN ) ;
                        %             Vcrel_mean_EIL100 = vertcat(Vcrel_mean_EIL100, NaN ) ;
                        %         end
                        %




                        %                     if(isnan(eil50(1)) == 0 & isnan(eil50(2)) == 0)
                        %                         kBotEil50  =  find( abs( ht - eil50(1) ) == min(abs( ht - eil50(1) )  ) ) ; kBotEil50 = kBotEil50(1);
                        %                         kTopEil50  =  find( abs( ht - eil50(2) ) == min(abs( ht - eil50(2) )  ) ) ; kTopEil50 = kTopEil50(1);
                        %                         %EIL50_mean_U = vertcat(EIL50_mean_U,  LAYERMEAN(u(:,t),ht,kBotEil50,kTopEil50) );
                        %                         %EIL50_mean_V = vertcat(EIL50_mean_V,  LAYERMEAN(v(:,t),ht,kBotEil50,kTopEil50) );
                        %                         rvap_mean_EIL50(x,y,t) =  LAYERMEAN(mr(:),ht,kBotEil50,kTopEil50) ;
                        %                         %EIL50_mean_RH = vertcat(EIL50_mean_RH, LAYERMEAN(rh(:,t),ht,kBotEil50,kTopEil50) );
                        %                         Ucrel_mean_EIL50(x,y,t) =  LAYERMEAN(u_crel(:),ht,kBotEil50,kTopEil50) ;
                        %                         Vcrel_mean_EIL50(x,y,t) =  LAYERMEAN(v_crel(:),ht,kBotEil50,kTopEil50) ;
                        %                     else
                        %                         %EIL50_mean_U = vertcat(EIL50_mean_U,  NaN );
                        %                         %EIL50_mean_V = vertcat(EIL50_mean_V,  NaN );
                        %                         rvap_mean_EIL50(x,y,t) =  NaN ;
                        %                         %EIL50_mean_RH = vertcat(EIL50_mean_RH, NaN );
                        %                         Ucrel_mean_EIL50(x,y,t) =  NaN;
                        %                         Vcrel_mean_EIL50(x,y,t) =  NaN;
                        %                     end

                        %if(isempty(adv_ind) == 0 & noadv == 0 & isnan(eil10(1)) == 0 & isnan(eil10(2)) == 0)
                        %                     Ucrel_mean_EIL50(x,y,t) =  LAYERMEAN(u_crel(:),ht,kBotEil50,kTopEil50) );
                        %                     Vcrel_mean_EIL50(x,y,t) =  LAYERMEAN(v_crel(:),ht,kBotEil50,kTopEil50) );
                        % end

                        %                              else
                        %                                 Ucrel_mean_EIL50 = vertcat(Ucrel_mean_EIL50, NaN ) ;
                        %                                 Vcrel_mean_EIL50 = vertcat(Vcrel_mean_EIL50, NaN ) ;
                        %                              end




                        %                     if(isnan(eil25(1)) == 0 & isnan(eil25(2)) == 0)
                        %                         kBotEil25  =  find( abs( ht - eil25(1) ) == min(abs( ht - eil25(1) )  ) ) ; kBotEil25 = kBotEil25(1);
                        %                         kTopEil25  =  find( abs( ht - eil25(2) ) == min(abs( ht - eil25(2) )  ) ) ; kTopEil25 = kTopEil25(1);
                        %                         %EIL25_mean_U = vertcat(EIL25_mean_U,  LAYERMEAN(u(:,t),ht,kBotEil25,kTopEil25) );
                        %                         %EIL25_mean_V = vertcat(EIL25_mean_V,  LAYERMEAN(v(:,t),ht,kBotEil25,kTopEil25) );
                        %                         rvap_mean_EIL25(x,y,t) =  LAYERMEAN(mr(:),ht,kBotEil25,kTopEil25) ;
                        %                         %EIL25_mean_RH = vertcat(EIL25_mean_RH, LAYERMEAN(rh(:,t),ht,kBotEil25,kTopEil25) );
                        %                         Ucrel_mean_EIL25(x,y,t) =  LAYERMEAN(u_crel(:),ht,kBotEil25,kTopEil25) ;
                        %                         Vcrel_mean_EIL25(x,y,t) =  LAYERMEAN(v_crel(:),ht,kBotEil25,kTopEil25) ;
                        %                     else
                        %                         %EIL25_mean_U = vertcat(EIL25_mean_U,  NaN );
                        %                         %EIL25_mean_V = vertcat(EIL25_mean_V,  NaN );
                        %                         rvap_mean_EIL25(x,y,t) = NaN ;
                        %                         %EIL25_mean_RH = vertcat(EIL25_mean_RH, NaN );
                        %                         Ucrel_mean_EIL25(x,y,t) =  NaN;
                        %                         Vcrel_mean_EIL25(x,y,t) =  NaN;
                        %                     end


                        %
                        %          if(isempty(adv_ind) == 0 & noadv == 0 & isnan(eil10(1)) == 0 & isnan(eil10(2)) == 0)
                        %             Ucrel_mean_EIL25 = vertcat(Ucrel_mean_EIL25, LAYERMEAN(u_crel(:),ht,kBotEil25,kTopEil25) );
                        %             Vcrel_mean_EIL25 = vertcat(Vcrel_mean_EIL25, LAYERMEAN(v_crel(:),ht,kBotEil25,kTopEil25) );
                        %          else
                        %             Ucrel_mean_EIL25 = vertcat(Ucrel_mean_EIL25, NaN ) ;
                        %             Vcrel_mean_EIL25 = vertcat(Vcrel_mean_EIL25, NaN ) ;
                        %          end
                        %



                        %                     if(isnan(eil10(1)) == 0 & isnan(eil10(2)) == 0)
                        %                         kBotEil10  =  find( abs( ht - eil10(1) ) == min(abs( ht - eil10(1) )  ) ) ; kBotEil10 = kBotEil10(1);
                        %                         kTopEil10  =  find( abs( ht - eil10(2) ) == min(abs( ht - eil10(2) )  ) ) ; kTopEil10 = kTopEil10(1);
                        %                         %EIL10_mean_U = vertcat(EIL10_mean_U,  LAYERMEAN(u(:,t),ht,kBotEil10,kTopEil10) );
                        %                         %EIL10_mean_V = vertcat(EIL10_mean_V,  LAYERMEAN(v(:,t),ht,kBotEil10,kTopEil10) );
                        %                         rvap_mean_EIL10(x,y,t) = LAYERMEAN(mr(:),ht,kBotEil10,kTopEil10) ;
                        %                         %EIL10_mean_RH = vertcat(EIL10_mean_RH,  LAYERMEAN(rh(:,t),ht,kBotEil10,kTopEil10) );
                        %                         Ucrel_mean_EIL10(x,y,t) =  LAYERMEAN(u_crel(:),ht,kBotEil10,kTopEil10) ;
                        %                         Vcrel_mean_EIL10(x,y,t) =  LAYERMEAN(v_crel(:),ht,kBotEil10,kTopEil10) ;
                        %                     else
                        %                         %EIL10_mean_U = vertcat(EIL10_mean_U,  NaN );
                        %                         %EIL10_mean_V = vertcat(EIL10_mean_V,  NaN );
                        %                         rvap_mean_EIL10(x,y,t) = NaN ;
                        %                         %EIL10_mean_RH = vertcat(EIL10_mean_RH, NaN );
                        %                         Ucrel_mean_EIL10(x,y,t) =  NaN;
                        %                         Vcrel_mean_EIL10(x,y,t) =  NaN;
                        %                     end
                        %
                        %          if(isempty(adv_ind) == 0 & noadv == 0 & isnan(eil10(1)) == 0 & isnan(eil10(2)) == 0)
                        %                 Ucrel_mean_EIL10 = vertcat(Ucrel_mean_EIL10, LAYERMEAN(u_crel(:),ht,kBotEil10,kTopEil10) );
                        %                 Vcrel_mean_EIL10 = vertcat(Vcrel_mean_EIL10, LAYERMEAN(v_crel(:),ht,kBotEil10,kTopEil10) );
                        %          else
                        %                 Ucrel_mean_EIL10 = vertcat(Ucrel_mean_EIL10, NaN ) ;
                        %                 Vcrel_mean_EIL10 = vertcat(Vcrel_mean_EIL10, NaN ) ;
                        %          end



                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% surface-based LFC/CAPE/CIN/EL-dependent variables





                        if(isnan(sfcLFC) == 1 )  % lfc doesnt exist
                            %
                            %             VAPmix_diff_grndtoACBLsfc = vertcat(VAPmix_diff_grndtoACBLsfc, NaN );
                            %             ACBLmax_sfc_LapseRate = vertcat( ACBLmax_sfc_LapseRate, NaN );                     %steepest lapse rate in ACBL
                            %             ACBLmean_sfc_LapseRate = vertcat( ACBLmean_sfc_LapseRate, NaN );        %mean lapse rate in ACBL
                            %             shearmag_mean_acbl_sfc  = vertcat( shearmag_mean_acbl_sfc, NaN );
                            %             shearmag_bulk_ACBLsfc = vertcat( shearmag_bulk_ACBLsfc, NaN );
                            %             V_mean_ACBLsfc = vertcat(V_mean_ACBLsfc,NaN);
                            %             U_mean_ACBLsfc = vertcat(U_mean_ACBLsfc,NaN);
                            %             pres_ACBL_sfc = vertcat( pres_ACBL_sfc, NaN );
                            %
                            %             LCL_temp_sfc = vertcat(LCL_temp_sfc,NaN);
                            %             LFC_temp_sfc = vertcat(LFC_temp_sfc,NaN);
                            %             EL_temp_sfc = vertcat(EL_temp_sfc,NaN);
                            %
                            %                                 Ucrel_mean_ACBL_sfc(x,y,t)    = NaN;
                            %                                 Vcrel_mean_ACBL_sfc(x,y,t)    = NaN;
                            %                                 Ucrel_mean_0toACBL_sfc(x,y,t) = NaN;
                            %                                 Vcrel_mean_0toACBL_sfc(x,y,t) = NaN;

                            %             Ucrel_mean_anvil_sfc = vertcat(Ucrel_mean_anvil_sfc ,NaN);
                            %             Vcrel_mean_anvil_sfc = vertcat(Vcrel_mean_anvil_sfc ,NaN);
                            %
                            %             shear_dir_bulk_acbl_sfc = vertcat(shear_dir_bulk_acbl_sfc, NaN);
                            %
                            %             shear_mag_bulk_FT_sfc = vertcat(shear_mag_bulk_FT_sfc, NaN); %jnm_new_sep8
                            %             shear_dir_bulk_FT_sfc = vertcat(shear_dir_bulk_FT_sfc, NaN); %jnm_new_sep8
                            %
                            %
                            %             shearmag_bulk_LFCp1_sfc = vertcat( shearmag_bulk_LFCp1_sfc, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp1_sfc = vertcat( shear_dir_bulk_LFCp1_sfc, NaN ); %jnm new_22oct
                            %
                            %             shearmag_bulk_LFCp2_sfc = vertcat( shearmag_bulk_LFCp2_sfc, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp2_sfc = vertcat( shear_dir_bulk_LFCp2_sfc, NaN ); %jnm new_22oct
                            %
                            %             shearmag_bulk_LFCp3_sfc = vertcat( shearmag_bulk_LFCp3_sfc, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp3_sfc = vertcat( shear_dir_bulk_LFCp3_sfc, NaN ); %jnm new_22oct
                            %
                            %             shearmag_bulk_LFCp4_sfc = vertcat( shearmag_bulk_LFCp4_sfc, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp4_sfc = vertcat( shear_dir_bulk_LFCp4_sfc, NaN ); %jnm new_22oct
                            %
                            %             RH_mean_LFCp1_sfc = vertcat(RH_mean_LFCp1_sfc, NaN ); %jnm new_7dec
                            %             RH_mean_LFCp2_sfc = vertcat(RH_mean_LFCp2_sfc, NaN ); %jnm new_7dec
                            %             RH_mean_LFCp3_sfc = vertcat(RH_mean_LFCp3_sfc, NaN ); %jnm new_7dec
                            %             RH_mean_LFCp4_sfc = vertcat(RH_mean_LFCp4_sfc, NaN ); %jnm new_7dec
                            %
                            %
                        elseif(isnan(sfcLFC) == 0 ) % lfc does exist
                            %
                            %             kELsfc =  find( abs(ht - sfcEL) == min(abs(ht - sfcEL)  ) ) ;  kELsfc = kELsfc(1);   %ht index of EL
                            %             kELsfcm1 = find( abs(ht - (sfcEL - 1000)) == min(abs(ht - (sfcEL - 1000))  ) ) ;  kELsfcm1 = kELsfcm1(1);   %ht index at 1 km below EL
                            kLFCsfc =  find( abs(ht - sfcLFC) == min(abs(ht - sfcLFC)  ) ) ;  kLFCsfc = kLFCsfc(1);     %ht index of lfc
                            kACBLsfc =  find( abs(ht - (sfcLFC+1500)) == min(abs(ht - (sfcLFC+1500))  ) ) ;  kACBLsfc = kACBLsfc(1);  %find TOP of active cloud bearing layer (LFC + 1.5 km)
                            %             kLCLsfc = find( abs(ht - sfcLCL) == min(abs(ht - sfcLCL)  ) ) ;  kLCLsfc = kLCLsfc(1);   %ht index of LCL
                            %
                            %
                            %             %new_22 oct 2021:
                            %             kLFCp1sfc =  find( abs(ht - (sfcLFC+1000)) == min(abs(ht - (sfcLFC+1000))  ) ) ;  kLFCp1sfc = kLFCp1sfc(1);  %find k index of LFC + 1 km
                            %             kLFCp2sfc =  find( abs(ht - (sfcLFC+2000)) == min(abs(ht - (sfcLFC+2000))  ) ) ;  kLFCp2sfc = kLFCp2sfc(1);  %find k index of LFC + 2 km
                            %             kLFCp3sfc =  find( abs(ht - (sfcLFC+3000)) == min(abs(ht - (sfcLFC+3000))  ) ) ;  kLFCp3sfc = kLFCp3sfc(1);  %find k index of LFC + 3 km
                            %             kLFCp4sfc =  find( abs(ht - (sfcLFC+4000)) == min(abs(ht - (sfcLFC+4000))  ) ) ;  kLFCp4sfc = kLFCp4sfc(1);  %find k index of LFC + 4 km
                            %
                            %
                            %             LCL_temp_sfc = vertcat(LCL_temp_sfc,tdry(kLCLsfc,t)+273.15);
                            %             LFC_temp_sfc = vertcat(LFC_temp_sfc,tdry(kLFCsfc,t)+273.15);
                            %             EL_temp_sfc  = vertcat(EL_temp_sfc, tdry(kELsfc,t)+273.15);
                            %
                            %             pres_ACBL_sfc = vertcat( pres_ACBL_sfc, pres(kACBLsfc) );
                            %
                            %
                            %
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %             %%%  cloud-rel wind metrics
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %             if(isempty(adv_ind) == 0 & noadv == 0)
                            %                                     Ucrel_mean_ACBL_sfc(x,y,t) = LAYERMEAN( u_crel(:),ht,kLFCsfc,kACBLsfc );
                            %                                     Vcrel_mean_ACBL_sfc(x,y,t) = LAYERMEAN( v_crel(:),ht,kLFCsfc,kACBLsfc );
                            %
                            %                                     Ucrel_mean_0toACBL_sfc(x,y,t) = LAYERMEAN( u_crel(:),ht,1,kACBLsfc );
                            %                                     Vcrel_mean_0toACBL_sfc(x,y,t) = LAYERMEAN( v_crel(:),ht,1,kACBLsfc );
                            %                     %
                            %                 Ucrel_mean_anvil_sfc = vertcat(Ucrel_mean_anvil_sfc , LAYERMEAN( u_crel(:),ht,kELsfcm1,kELsfc ));
                            %                 Vcrel_mean_anvil_sfc = vertcat(Vcrel_mean_anvil_sfc , LAYERMEAN( v_crel(:),ht,kELsfcm1,kELsfc ));
                            %             elseif( isempty(adv_ind) == 1 | noadv == 1 )
                            %                 Ucrel_mean_ACBL_sfc = vertcat(Ucrel_mean_ACBL_sfc , NaN);
                            %                 Vcrel_mean_ACBL_sfc = vertcat(Vcrel_mean_ACBL_sfc , NaN);
                            %                 Ucrel_mean_0toACBL_sfc = vertcat(Ucrel_mean_0toACBL_sfc , NaN );
                            %                 Vcrel_mean_0toACBL_sfc = vertcat(Vcrel_mean_0toACBL_sfc , NaN );
                            %                 Ucrel_mean_anvil_sfc = vertcat(Ucrel_mean_anvil_sfc , NaN);
                            %                 Vcrel_mean_anvil_sfc = vertcat(Vcrel_mean_anvil_sfc , NaN);
                            %             end
                            %
                            %
                            %
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %             %%%  humidity metrics
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %             VAPmix_diff_grndtoACBLsfc = vertcat(VAPmix_diff_grndtoACBLsfc, 1000*(mr(3,t) - mr(kACBLsfc,t)) );   %g/kg
                            %
                            %             RH_mean_LFCp1_sfc = vertcat( RH_mean_LFCp1_sfc, LAYERMEAN(rh(:,t),ht,kLFCsfc,kLFCp1sfc) );  % jnm new_7dec
                            %             RH_mean_LFCp2_sfc = vertcat( RH_mean_LFCp2_sfc, LAYERMEAN(rh(:,t),ht,kLFCp1sfc,kLFCp2sfc) );  % jnm new_7dec
                            %             RH_mean_LFCp3_sfc = vertcat( RH_mean_LFCp3_sfc, LAYERMEAN(rh(:,t),ht,kLFCp2sfc,kLFCp3sfc) );  % jnm new_7dec
                            %             RH_mean_LFCp4_sfc = vertcat( RH_mean_LFCp4_sfc, LAYERMEAN(rh(:,t),ht,kLFCp3sfc,kLFCp4sfc) );  % jnm new_7dec
                            %
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %             %%%  lapse rate metrics
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %             ACBLmax_sfc_LapseRate = vertcat( ACBLmax_sfc_LapseRate, min( LRtv(kLFCsfc:kACBLsfc)) );                     %steepest lapse rate in ACBL
                            %             ACBLmean_sfc_LapseRate = vertcat( ACBLmean_sfc_LapseRate, LAYERMEAN( LRtv(:),ht,kLFCsfc,kACBLsfc ));        %mean lapse rate in ACBL
                            %
                            %
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %             %%%  shear metrics
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %             V_mean_ACBLsfc = vertcat(V_mean_ACBLsfc, LAYERMEAN( v(:,t),ht,kLFCsfc,kACBLsfc) );
                            %             U_mean_ACBLsfc = vertcat(U_mean_ACBLsfc, LAYERMEAN( u(:,t),ht,kLFCsfc,kACBLsfc) );
                            %
                            %             shearmag_mean_acbl_sfc = vertcat( shearmag_mean_acbl_sfc, LAYERMEAN( magshear(:),ht,kLFCsfc,kACBLsfc ) );
                            %
                            %             du = u(kACBLsfc,t) - u(kLFCsfc,t); dv = v(kACBLsfc,t) - v(kLFCsfc,t);
                            %             shearmag_bulk_ACBLsfc = vertcat( shearmag_bulk_ACBLsfc, ( du*du + dv*dv ).^0.5 );
                            %
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_acbl_sfc = vertcat(shear_dir_bulk_acbl_sfc, sdir);
                            %
                            %
                            %             du = u(kELsfc,t) - u(kLFCsfc,t); dv = v(kELsfc,t) - v(kLFCsfc,t);                   %jnm_new_sep8
                            %             shear_mag_bulk_FT_sfc = vertcat(shear_mag_bulk_FT_sfc,  ( du*du + dv*dv ).^0.5 );   %jnm_new_sep8
                            %
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;                                             %jnm_new_sep8
                            %             if( sdir > 360.0 )                                                                  %jnm_new_sep8
                            %                 sdir = sdir - 360 ;                                                             %jnm_new_sep8
                            %             end                                                                                 %jnm_new_sep8
                            %             shear_dir_bulk_FT_sfc = vertcat(shear_dir_bulk_FT_sfc, sdir);                       %jnm_new_sep8
                            %
                            %
                            %
                            %
                            %
                            %             % % new_22 oct 2021:
                            %             du = u(kLFCp1sfc,t) - u(kLFCsfc ,t); dv = v(kLFCp1sfc,t) - v(kLFCsfc ,t);
                            %             shearmag_bulk_LFCp1_sfc = vertcat( shearmag_bulk_LFCp1_sfc, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp1_sfc = vertcat(shear_dir_bulk_LFCp1_sfc, sdir);
                            %
                            %
                            %             du = u(kLFCp2sfc,t) - u(kLFCp1sfc ,t); dv = v(kLFCp2sfc,t) - v(kLFCp1sfc ,t);
                            %             shearmag_bulk_LFCp2_sfc = vertcat( shearmag_bulk_LFCp2_sfc, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp2_sfc = vertcat(shear_dir_bulk_LFCp2_sfc, sdir);
                            %
                            %
                            %             du = u(kLFCp3sfc,t) - u(kLFCp2sfc ,t); dv = v(kLFCp3sfc,t) - v(kLFCp2sfc ,t);
                            %             shearmag_bulk_LFCp3_sfc = vertcat( shearmag_bulk_LFCp3_sfc, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp3_sfc = vertcat(shear_dir_bulk_LFCp3_sfc, sdir);
                            %
                            %
                            %             du = u(kLFCp4sfc,t) - u(kLFCp3sfc ,t); dv = v(kLFCp4sfc,t) - v(kLFCp3sfc ,t);
                            %             shearmag_bulk_LFCp4_sfc = vertcat( shearmag_bulk_LFCp4_sfc, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp4_sfc = vertcat(shear_dir_bulk_LFCp4_sfc, sdir);
                            %             %end new_22 Oct
                            %
                        end


                        %ok to here






                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% ML-based LFC/CAPE/CIN/EL-dependent variables




                        if(isnan(mlLFC)== 1 )
                            %
                            %             ACBLmax_ml_LapseRate  = vertcat( ACBLmax_ml_LapseRate,  NaN );                    %steepest lapse rate in ACBL
                            %             ACBLmean_ml_LapseRate  = vertcat( ACBLmean_ml_LapseRate,  NaN );      %mean lapse rate in ACBL
                            %
                            %             pres_ACBL_ml = vertcat( pres_ACBL_ml, NaN );
                            %
                            %             V_mean_ACBLml = vertcat(V_mean_ACBLml,NaN);
                            %             U_mean_ACBLml = vertcat(U_mean_ACBLml,NaN);
                            %
                            %                                 Ucrel_mean_ACBL_ml(x,y,t) = NaN;
                            %                                 Vcrel_mean_ACBL_ml(x,y,t) = NaN;
                            %
                            %                                 Ucrel_mean_0toACBL_ml(x,y,t) =  NaN;
                            %                                 Vcrel_mean_0toACBL_ml(x,y,t) =  NaN;
                            %
                            %             Ucrel_mean_anvil_ml = vertcat(Ucrel_mean_anvil_ml , NaN);
                            %             Vcrel_mean_anvil_ml = vertcat(Vcrel_mean_anvil_ml , NaN);
                            %
                            %             LCL_temp_ml = vertcat(LCL_temp_ml,NaN);
                            %             LFC_temp_ml = vertcat(LFC_temp_ml,NaN);
                            %             EL_temp_ml  = vertcat(EL_temp_ml,NaN);
                            %
                        elseif(isnan(mlLFC) == 0 )
                            %             kELml  =  find( abs(ht - mlEL) == min(abs(ht - mlEL)  ) ) ;   kELml = kELml(1);
                            %             kELmlm1 = find( abs(ht - (mlEL - 1000)) == min(abs(ht - (mlEL - 1000))  ) ) ;  kELmlm1 = kELmlm1(1);   %ht index at 1 km below EL
                            kLFCml  =  find( abs(ht - mlLFC) == min(abs(ht - mlLFC)  ) ) ;  kLFCml = kLFCml(1);
                            kACBLml  =  find( abs(ht - (mlLFC+1500)) == min(abs(ht - (mlLFC+1500))  ) ) ; kACBLml = kACBLml(1);
                            %             kLCLml = find( abs(ht - mlLCL) == min(abs(ht - mlLCL)  ) ) ;  kLCLml = kLCLml(1);   %ht index of LCL
                            %
                            %
                            %             LCL_temp_ml = vertcat(LCL_temp_ml,tdry(kLCLml,t)+273.15);
                            %             LFC_temp_ml = vertcat(LFC_temp_ml,tdry(kLFCml,t)+273.15);
                            %             EL_temp_ml  = vertcat(EL_temp_ml, tdry(kELml,t)+273.15);
                            %
                            %             pres_ACBL_ml = vertcat( pres_ACBL_ml, pres(kACBLml) );
                            %
                            %
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %             %%%  lapse rate metrics
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %             ACBLmax_ml_LapseRate  = vertcat( ACBLmax_ml_LapseRate,  min( LRtv(kLFCml:kACBLml)) );                    %steepest lapse rate in ACBL
                            %             ACBLmean_ml_LapseRate  = vertcat( ACBLmean_ml_LapseRate,  LAYERMEAN( LRtv(:),ht,kLFCml,kACBLml) );      %mean lapse rate in ACBL
                            %
                            %             %%%%%
                            %             % other
                            %             %%%%%
                            %
                            %             V_mean_ACBLml = vertcat(V_mean_ACBLml, LAYERMEAN( v(:,t),ht,kLFCml,kACBLml));
                            %             U_mean_ACBLml = vertcat(U_mean_ACBLml, LAYERMEAN( u(:,t),ht,kLFCml,kACBLml));
                            %
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %             %%%  cloud-rel wind metrics
                            %             %%%%%%%%%%%%%%%%%%%%%%%
                            %
                            %             if( isempty(adv_ind) == 0 & noadv == 0 )
                            %                                     Ucrel_mean_ACBL_ml(x,y,t) = LAYERMEAN( u_crel(:),ht,kLFCml,kACBLml );
                            %                                     Vcrel_mean_ACBL_ml(x,y,t) = LAYERMEAN( v_crel(:),ht,kLFCml,kACBLml );
                            %
                            %                                     Ucrel_mean_0toACBL_ml(x,y,t) = LAYERMEAN( u_crel(:),ht,1,kACBLml );
                            %                                     Vcrel_mean_0toACBL_ml(x,y,t) = LAYERMEAN( v_crel(:),ht,1,kACBLml );
                            %
                            %                 Ucrel_mean_anvil_ml = vertcat(Ucrel_mean_anvil_ml , LAYERMEAN( u_crel(:),ht,kELmlm1,kELml));
                            %                 Vcrel_mean_anvil_ml = vertcat(Vcrel_mean_anvil_ml , LAYERMEAN( v_crel(:),ht,kELmlm1,kELml));
                            %             elseif( isempty(adv_ind) == 1 | noadv == 1 )
                            %                 Ucrel_mean_ACBL_ml = vertcat(Ucrel_mean_ACBL_ml ,NaN);
                            %                 Vcrel_mean_ACBL_ml = vertcat(Vcrel_mean_ACBL_ml ,NaN);
                            %
                            %                 Ucrel_mean_0toACBL_ml = vertcat(Ucrel_mean_0toACBL_ml ,NaN);
                            %                 Vcrel_mean_0toACBL_ml = vertcat(Vcrel_mean_0toACBL_ml ,NaN);
                            %
                            %                 Ucrel_mean_anvil_ml = vertcat(Ucrel_mean_anvil_ml , NaN);
                            %                 Vcrel_mean_anvil_ml = vertcat(Vcrel_mean_anvil_ml , NaN);
                            %             end
                            %
                            %
                        end






                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables
                        %%%%%%%%%%%%%%% MU-based LFC/CAPE/CIN/EL-dependent variables




                        if(isnan(muLFC)== 1 )

                            %             FTmean_RH = vertcat( FTmean_RH, NaN );
                            %             FTmean_DewptDepr = vertcat( FTmean_DewptDepr, NaN );
                            %             FTmax_DewptDepr  = vertcat( FTmax_DewptDepr, NaN );  %the driest dewpt depress
                            %             thetae_mean_FT   = vertcat( thetae_mean_FT, NaN ); %new_jnm_23mar2022
                            RH_mean_ACBL(x,y,t,c) =  NaN ;
                            %             ACBLmin_RH         = vertcat( ACBLmin_RH, NaN );   %driest RH
                            %             ACBLmean_DewptDepr = vertcat( ACBLmean_DewptDepr, NaN );
                            %             ACBLmax_DewptDepr  = vertcat( ACBLmax_DewptDepr,  NaN );  %the driest dewpt depress
                            %             VAPmix_diff_grndtoACBLmu = vertcat(VAPmix_diff_grndtoACBLmu, NaN );
                            %             dZstar_mu = vertcat(dZstar_mu, NaN  );  %from Lock and Houston 2014 (diff between initial parcel ht and lfc)
                            %             ACBLmax_mu_LapseRate  = vertcat( ACBLmax_mu_LapseRate, NaN );     %steepest lapse rate in ACBL
                            %             ACBLmean_mu_LapseRate  = vertcat( ACBLmean_mu_LapseRate,  NaN);      %mean lapse rate in ACBL
                            %             shearmag_mean_subcloud_mu  = vertcat( shearmag_mean_subcloud_mu, NaN);
                            %             shearmag_mean_acbl_mu   = vertcat( shearmag_mean_acbl_mu, NaN);
                            %
                            shear_mag_bulk_ACBL_mu(x,y,t,c) = NaN;
                            %             shear_dir_bulk_acbl_mu = vertcat( shear_dir_bulk_acbl_mu, NaN);
                            %
                            %             pres_ACBL_mu = vertcat( pres_ACBL_mu, NaN );
                            %
                            %             V_mean_ACBLmu = vertcat(V_mean_ACBLmu,NaN);
                            %             U_mean_ACBLmu = vertcat(U_mean_ACBLmu,NaN);

            		        rh_at_muLFCplus1500m(x,y,t,c) = NaN;
                            Ucrel_mu(x,y,t,c) =  NaN ;
                            Vcrel_mu(x,y,t,c) =  NaN ;

                            U_mu(x,y,t,c) =  NaN ;
                            V_mu(x,y,t,c) =  NaN ;
                			theta_mu(x,y,t,c) = NaN;
                			thetav_mu(x,y,t,c) = NaN;

                            Ucrel_mean_ACBL_mu(x,y,t,c)    =  NaN;
                            Vcrel_mean_ACBL_mu(x,y,t,c)    = NaN;
                            %Ucrel_mean_0toACBL_mu(x,y,t) =  NaN;
                            %Vcrel_mean_0toACBL_mu(x,y,t) = NaN;

                            %             Ucrel_mean_anvil_mu = vertcat(Ucrel_mean_anvil_mu ,NaN);
                            %             Vcrel_mean_anvil_mu = vertcat(Vcrel_mean_anvil_mu ,NaN);

                            %             Ucrel_mean_FT = vertcat(Ucrel_mean_FT, NaN);
                            %             Vcrel_mean_FT = vertcat(Vcrel_mean_FT, NaN);
                            %             Ucrel_mean_subcloud_mu  = vertcat( Ucrel_mean_subcloud_mu, NaN);
                            %             Vcrel_mean_subcloud_mu  = vertcat( Vcrel_mean_subcloud_mu, NaN);

                            %             U_mu = vertcat(U_mu, NaN);
                            %             V_mu = vertcat(V_mu, NaN);
                            %             theta_mu = vertcat(theta_mu,NaN);
                            %             thetav_mu = vertcat(thetav_mu,NaN);

                            thetae_mu(x,y,t,c) = NaN; %jnm_new_23mar2022
                            thetae_mean_subcloud_mu(x,y,t,c) = NaN; %jnm_new_23mar2022
                            rvap_mu(x,y,t,c) =  NaN;

                            thetae_mean_ACBL_mu(x,y,t,c) =  NaN; %jnm_new_23mar2022

                            %                         LCL_temp_mu(x,y,t) = NaN;
                            %                         LFC_temp_mu(x,y,t) = NaN;
                            %                         EL_temp_mu(x,y,t)  = NaN;

                            shear_mag_bulk_FT_mu(x,y,t,c) =  NaN; %jnm_new_sep8
                            shear_dir_bulk_FT_mu(x,y,t,c) =  NaN; %jnm_new_sep8

                            %             shearmag_bulk_LFCp1_mu = vertcat( shearmag_bulk_LFCp1_mu, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp1_mu = vertcat( shear_dir_bulk_LFCp1_mu, NaN ); %jnm new_22oct
                            %
                            %             shearmag_bulk_LFCp2_mu = vertcat( shearmag_bulk_LFCp2_mu, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp2_mu = vertcat( shear_dir_bulk_LFCp2_mu, NaN ); %jnm new_22oct
                            %
                            %             shearmag_bulk_LFCp3_mu = vertcat( shearmag_bulk_LFCp3_mu, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp3_mu = vertcat( shear_dir_bulk_LFCp3_mu, NaN ); %jnm new_22oct
                            %
                            %             shearmag_bulk_LFCp4_mu = vertcat( shearmag_bulk_LFCp4_mu, NaN ); %jnm new_22oct
                            %             shear_dir_bulk_LFCp4_mu = vertcat( shear_dir_bulk_LFCp4_mu, NaN ); %jnm new_22oct


                            %             RH_mean_LFCp1_mu = vertcat(RH_mean_LFCp1_mu, NaN ); %jnm new_7dec
                            %             RH_mean_LFCp2_mu = vertcat(RH_mean_LFCp2_mu, NaN ); %jnm new_7dec
                            %             RH_mean_LFCp3_mu = vertcat(RH_mean_LFCp3_mu, NaN ); %jnm new_7dec
                            %             RH_mean_LFCp4_mu = vertcat(RH_mean_LFCp4_mu, NaN ); %jnm new_7dec

                        elseif(isnan(muLFC) == 0 )  %there is an mu lfc

                            kMU    = find( abs(ht - hMU) == min(abs(ht - hMU)  ) ) ;                        kMU = kMU(1);
                            kELmu  =  find( abs(ht - muEL) == min(abs(ht - muEL)  ) ) ;                     kELmu = kELmu(1);
                            kELmum1 = find( abs(ht - (muEL - 1000)) == min(abs(ht - (muEL - 1000))) ) ;       kELmum1 = kELmum1(1);   %ht index at 1 km below EL
                            kLFCmu  =  find( abs(ht - muLFC) == min(abs(ht - muLFC)  ) ) ;                  kLFCmu = kLFCmu(1);
                            kACBLmu  =  find( abs(ht - (muLFC+1500)) == min(abs(ht - (muLFC+1500))  ) ) ;   kACBLmu = kACBLmu(1);
                            kmuLCL  =  find( abs(ht - muLCL) == min(abs(ht - muLCL)  ) ) ;                  kmuLCL = kmuLCL(1);

                            %new_22 oct 2021:
                            %             kLFCp1mu =  find( abs(ht - (muLFC+1000)) == min(abs(ht - (muLFC+1000))  ) ) ;  kLFCp1mu = kLFCp1mu(1);  %find k index of LFC + 1 km
                            %             kLFCp2mu =  find( abs(ht - (muLFC+2000)) == min(abs(ht - (muLFC+2000))  ) ) ;  kLFCp2mu = kLFCp2mu(1);  %find k index of LFC + 2 km
                            %             kLFCp3mu =  find( abs(ht - (muLFC+3000)) == min(abs(ht - (muLFC+3000))  ) ) ;  kLFCp3mu = kLFCp3mu(1);  %find k index of LFC + 3 km
                            %             kLFCp4mu =  find( abs(ht - (muLFC+4000)) == min(abs(ht - (muLFC+4000))  ) ) ;  kLFCp4mu = kLFCp4mu(1);  %find k index of LFC + 4 km


                            %            pres_ACBL_mu = vertcat( pres_ACBL_mu, pres(kACBLmu) );

                            %                         LCL_temp_mu(x,y,t) = tdry(kmuLCL)+273.15;
                            %                         LFC_temp_mu(x,y,t) = tdry(kLFCmu)+273.15;
                            %                         EL_temp_mu(x,y,t)  = tdry(kELmu)+273.15;

                            %%%%%%%%%%%%%%%%%%%%%%%
                            %%%  cloud-rel wind metrics
                            %%%%%%%%%%%%%%%%%%%%%%%

                            if( isnan(cell_motion_x)==0 & isnan(cell_motion_y)==0 )

                                Ucrel_mean_ACBL_mu(x,y,t,c) = LAYERMEAN( u_crel(:),ht,kLFCmu,kACBLmu);
                                Vcrel_mean_ACBL_mu(x,y,t,c) = LAYERMEAN( v_crel(:),ht,kLFCmu,kACBLmu);

                                %Ucrel_mean_0toACBL_mu(x,y,t) = LAYERMEAN( u_crel(:),ht,kMU,kACBLmu);
                                %Vcrel_mean_0toACBL_mu(x,y,t) = LAYERMEAN( v_crel(:),ht,kMU,kACBLmu);

                                Ucrel_mu(x,y,t,c) =  mean( u_crel(kMU),'omitnan');
                                Vcrel_mu(x,y,t,c) =  mean( v_crel(kMU),'omitnan');

                                %                 Ucrel_mean_anvil_mu = vertcat(Ucrel_mean_anvil_mu , LAYERMEAN( u_crel(:),ht,kELmum1,kELmu));
                                %                 Vcrel_mean_anvil_mu = vertcat(Vcrel_mean_anvil_mu , LAYERMEAN( v_crel(:),ht,kELmum1,kELmu));

                            elseif( isnan(cell_motion_x) | isnan(cell_motion_y) )

                                Ucrel_mean_ACBL_mu(x,y,t,c) = NaN;
                                Vcrel_mean_ACBL_mu(x,y,t,c) = NaN;

                                %Ucrel_mean_0toACBL_mu(x,y,t) = NaN;
                                %Vcrel_mean_0toACBL_mu(x,y,t) = NaN;

                                Ucrel_mu(x,y,t,c) =  NaN;
                                Vcrel_mu(x,y,t,c) =  NaN;

                                %                 Ucrel_mean_anvil_mu = vertcat(Ucrel_mean_anvil_mu , NaN );
                                %                 Vcrel_mean_anvil_mu = vertcat(Vcrel_mean_anvil_mu , NaN );

                            end



                            %%%%%%%%%%%%%%%%%%%%%%%
                            %%%  humidity metrics
                            %%%%%%%%%%%%%%%%%%%%%%%

                            %             RH_mean_LFCp1_mu = vertcat( RH_mean_LFCp1_mu, LAYERMEAN(rh(:,t),ht,kLFCmu,kLFCp1mu) );  % jnm new_7dec
                            %             RH_mean_LFCp2_mu = vertcat( RH_mean_LFCp2_mu, LAYERMEAN(rh(:,t),ht,kLFCp1mu,kLFCp2mu) );  % jnm new_7dec
                            %             RH_mean_LFCp3_mu = vertcat( RH_mean_LFCp3_mu, LAYERMEAN(rh(:,t),ht,kLFCp2mu,kLFCp3mu) );  % jnm new_7dec
                            %             RH_mean_LFCp4_mu = vertcat( RH_mean_LFCp4_mu, LAYERMEAN(rh(:,t),ht,kLFCp3mu,kLFCp4mu) );  % jnm new_7dec

                            %             if(isnan(BLtop_LL) == 0)
                            %                 FTmean_RH = vertcat( FTmean_RH, mean(rh(kBLtop+1:kELmu,t),'omitnan') );
                            %                 FTmean_DewptDepr = vertcat( FTmean_DewptDepr, mean(TdDepress(kBLtop+1:kELmu,t),'omitnan') );
                            %                 FTmax_DewptDepr = vertcat( FTmax_DewptDepr, max(TdDepress(kBLtop+1:kELmu,t)) );  %the driest dewpt depress
                            %                 xDD = max(TdDepress(kBLtop+1:kELmu,t))
                            %                 mDD = mean(TdDepress(kBLtop+1:kELmu,t),'omitnan')
                            %             elseif(isnan(BLtop_LL))
                            %                 FTmean_RH = vertcat( FTmean_RH, NaN );
                            %                 FTmean_DewptDepr = vertcat( FTmean_DewptDepr, NaN );
                            %                 FTmax_DewptDepr = vertcat( FTmax_DewptDepr, NaN );  %the driest dewpt depress
                            %                 xDD = NaN
                            %                 mDD = NaN
                            %             end

                            %            if(isnan(BLtop_LL) == 0)
                            %                if( (kBLtop+1) < kELmu) % the occassion with EL is really low and/or BL is freakishly high, which breaks some calculations
                            %                     FTmean_RH = vertcat( FTmean_RH, LAYERMEAN(rh(:,t),ht,kBLtop+1,kELmu) );
                            %                     FTmean_DewptDepr = vertcat( FTmean_DewptDepr, LAYERMEAN(TdDepress(:,t),ht,kBLtop+1,kELmu) );
                            %                     FTmax_DewptDepr = vertcat( FTmax_DewptDepr, max(TdDepress(kBLtop+1:kELmu,t)) );  %the driest dewpt depress
                            %                     thetae_mean_FT   = vertcat( thetae_mean_FT, LAYERMEAN(thetae(:,t),ht,kBLtop+1,kELmu) ); %new_jnm_23mar2022
                            %
                            %                     Ucrel_mean_FT = vertcat(Ucrel_mean_FT, LAYERMEAN(u_crel(:),ht,kBLtop+1,kELmu ) );
                            %                     Vcrel_mean_FT = vertcat(Vcrel_mean_FT, LAYERMEAN(v_crel(:),ht,kBLtop+1,kELmu ) );
                            %
                            %
                            %                 else
                            %                     FTmean_RH = vertcat( FTmean_RH, NaN );
                            %                     FTmean_DewptDepr = vertcat( FTmean_DewptDepr, NaN );
                            %                     FTmax_DewptDepr = vertcat( FTmax_DewptDepr, NaN );  %the driest dewpt depress
                            %                     thetae_mean_FT   = vertcat( thetae_mean_FT, NaN); %new_jnm_23mar2022
                            %
                            %                     Ucrel_mean_FT = vertcat(Ucrel_mean_FT, NaN );
                            %                     Vcrel_mean_FT = vertcat(Vcrel_mean_FT, NaN );
                            %
                            %                 end
                            %             else
                            %                 FTmean_RH = vertcat( FTmean_RH, NaN );
                            %                 FTmean_DewptDepr = vertcat( FTmean_DewptDepr, NaN );
                            %                 FTmax_DewptDepr = vertcat( FTmax_DewptDepr, NaN );  %the driest dewpt depress
                            %                 thetae_mean_FT   = vertcat( thetae_mean_FT, NaN); %new_jnm_23mar2022
                            %
                            %                 Ucrel_mean_FT = vertcat(Ucrel_mean_FT, NaN );
                            %                 Vcrel_mean_FT = vertcat(Vcrel_mean_FT, NaN );
                            %
                            %             end



                            %ACBL humidity:
                			rh_at_muLFCplus1500m(x,y,t,c) = rh(kACBLmu);
                            RH_mean_ACBL(x,y,t,c)        =  LAYERMEAN(rh(:), ht, kLFCmu, kACBLmu);
                            %             ACBLmin_RH         = vertcat( ACBLmin_RH,         min(rh(kLFCmu:kACBLmu,t)) );   %driest RH
                            %             ACBLmean_DewptDepr = vertcat( ACBLmean_DewptDepr, LAYERMEAN(TdDepress(:,t),ht,kLFCmu,kACBLmu) );
                            %             ACBLmax_DewptDepr  = vertcat( ACBLmax_DewptDepr,  max(TdDepress(kLFCmu:kACBLmu,t)) );  %the driest dewpt depress

                            %             if(isnan(hMU) == 0)
                            %                 VAPmix_diff_grndtoACBLmu = vertcat(VAPmix_diff_grndtoACBLmu, 1000*(mr(kMU,t) - mr(kACBLmu,t)) );  %g/kg
                            %             elseif(isnan(hMU) == 1)
                            %                 VAPmix_diff_grndtoACBLmu = vertcat(VAPmix_diff_grndtoACBLmu, NaN );  %g/kg
                            %             end

                            %%%%%%%%%%%%%%%%%%%%%%%
                            %%%  etc. metrics
                            %%%%%%%%%%%%%%%%%%%%%%%

                            if(isnan(hMU) == 0)
                                %dZstar_mu = vertcat(dZstar_mu, muLFC - hMU  );  %from Lock and Houston 2014 (diff between initial parcel ht and lfc)
                                thetae_mean_ACBL_mu(x,y,t,c) =  LAYERMEAN(thetae(:), ht, kLFCmu, kACBLmu) ; %jnm_new_23mar2022
                                thetae_mean_subcloud_mu(x,y,t,c)  =  LAYERMEAN(thetae(:),ht,kMU,kmuLCL ) ;%jnm_new_23mar2022
                            elseif(isnan(hMU) == 1)
                                %dZstar_mu = vertcat(dZstar_mu, NaN  );
                                thetae_mean_ACBL_mu(x,y,t,c) =  NaN; %jnm_new_23mar2022
                                thetae_mean_subcloud_mu(x,y,t,c)  =  NaN ; %jnm_new_23mar2022
                            end



                            %%%%%%%%%%%%%%%%%%%%%%%
                            %%%  lapse rate metrics
                            %%%%%%%%%%%%%%%%%%%%%%%

                            %             ACBLmax_mu_LapseRate  = vertcat( ACBLmax_mu_LapseRate,  min( LRtv(kLFCmu:kACBLmu)) );     %steepest lapse rate in ACBL
                            %             ACBLmean_mu_LapseRate  = vertcat( ACBLmean_mu_LapseRate,  LAYERMEAN( LRtv(:),ht,kLFCmu,kACBLmu) );      %mean lapse rate in ACBL


                            %%%%%%%%%%%%%%%%%%%%%%%
                            %%%  shear metrics
                            %%%%%%%%%%%%%%%%%%%%%%%
                            %             V_mean_ACBLmu = vertcat(V_mean_ACBLmu, LAYERMEAN( v(:,t),ht,kLFCmu,kACBLmu));
                            %             U_mean_ACBLmu = vertcat(U_mean_ACBLmu, LAYERMEAN( u(:,t),ht,kLFCmu,kACBLmu));
                            %
                            %             if(isnan(hMU) == 0)
                            %                 shearmag_mean_subcloud_mu  = vertcat( shearmag_mean_subcloud_mu,  LAYERMEAN( magshear(:),ht,kMU,kmuLCL ) );
                            %
                            %                 Ucrel_mean_subcloud_mu  = vertcat( Ucrel_mean_subcloud_mu, LAYERMEAN(u_crel(:),ht,kMU,kmuLCL ));
                            %                 Vcrel_mean_subcloud_mu  = vertcat( Vcrel_mean_subcloud_mu, LAYERMEAN(u_crel(:),ht,kMU,kmuLCL ));
                            %
                            %             elseif(isnan(hMU) == 1)
                            %                 shearmag_mean_subcloud_mu  = vertcat( shearmag_mean_subcloud_mu,  NaN );
                            %
                            %                 Ucrel_mean_subcloud_mu  = vertcat( Ucrel_mean_subcloud_mu, NaN);
                            %                 Vcrel_mean_subcloud_mu  = vertcat( Vcrel_mean_subcloud_mu, NaN);
                            %             end
                            %
                            %             shearmag_mean_acbl_mu   = vertcat( shearmag_mean_acbl_mu,  LAYERMEAN( magshear(:),ht,kLFCmu,kACBLmu) );



                            du = U(kACBLmu) - U(kLFCmu); dv = V(kACBLmu) - V(kLFCmu);
                            shear_mag_bulk_ACBL_mu(x,y,t,c) =  ( du*du + dv*dv ).^0.5 ;
                            %
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_acbl_mu = vertcat(shear_dir_bulk_acbl_mu, sdir);
                            %


                            du = U(kELmu) - U(kLFCmu); dv = V(kELmu) - V(kLFCmu);                       %jnm_new_sep8
                            shear_mag_bulk_FT_mu(x,y,t,c) =   ( du*du + dv*dv ).^0.5 ;                       %jnm_new_sep8

                            sdir = 270 - rad2deg( atan2(dv,du) )  ;                                             %jnm_new_sep8
                            if( sdir > 360.0 )                                                                  %jnm_new_sep8
                                sdir = sdir - 360 ;                                                             %jnm_new_sep8
                            end                                                                                 %jnm_new_sep8
                            shear_dir_bulk_FT_mu(x,y,t,c) = sdir;                         %jnm_new_sep8




                            %             % % new_22 oct 2021:
                            %             du = u(kLFCp1mu,t) - u(kLFCmu ,t); dv = v(kLFCp1mu,t) - v(kLFCmu ,t);
                            %             shearmag_bulk_LFCp1_mu = vertcat( shearmag_bulk_LFCp1_mu, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp1_mu = vertcat(shear_dir_bulk_LFCp1_mu, sdir);
                            %
                            %
                            %             du = u(kLFCp2mu,t) - u(kLFCp1mu ,t); dv = v(kLFCp2mu,t) - v(kLFCp1mu ,t);
                            %             shearmag_bulk_LFCp2_mu = vertcat( shearmag_bulk_LFCp2_mu, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp2_mu = vertcat(shear_dir_bulk_LFCp2_mu, sdir);
                            %
                            %
                            %             du = u(kLFCp3mu,t) - u(kLFCp2mu ,t); dv = v(kLFCp3mu,t) - v(kLFCp2mu ,t);
                            %             shearmag_bulk_LFCp3_mu = vertcat( shearmag_bulk_LFCp3_mu, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp3_mu = vertcat(shear_dir_bulk_LFCp3_mu, sdir);
                            %
                            %
                            %             du = u(kLFCp4mu,t) - u(kLFCp3mu ,t); dv = v(kLFCp4mu,t) - v(kLFCp3mu ,t);
                            %             shearmag_bulk_LFCp4_mu = vertcat( shearmag_bulk_LFCp4_mu, ( du*du + dv*dv ).^0.5 );
                            %             sdir = 270 - rad2deg( atan2(dv,du) )  ;
                            %             if( sdir > 360.0 )
                            %                 sdir = sdir - 360 ;
                            %             end
                            %             shear_dir_bulk_LFCp4_mu = vertcat(shear_dir_bulk_LFCp4_mu, sdir);
                            %             %end new_22 Oct


                            if(isnan(hMU) == 0)
                                %                 U_mu = vertcat(U_mu,u(kMU,t));
                                %                 V_mu = vertcat(V_mu,v(kMU,t));
                                %                 theta_mu = vertcat(theta_mu,theta(kMU,t));
                                %                 thetav_mu = vertcat(thetav_mu,thetav(kMU,t));
                                U_mu(x,y,t,c) =  U(kMU) ;
                                V_mu(x,y,t,c) =  V(kMU) ;
                                theta_mu(x,y,t,c) = theta(kMU);
                                thetav_mu(x,y,t,c) = thetav(kMU);
                			    thetae_mu(x,y,t,c) = thetae(kMU); %jnm_new_23mar2022
                                rvap_mu(x,y,t,c) = mr(kMU)*1000; %g/kg
                            elseif(isnan(hMU) == 1)
                                %                 U_mu = vertcat(U_mu, NaN);
                                %                 V_mu = vertcat(V_mu, NaN);
                                %                 theta_mu = vertcat(theta_mu,NaN);
                                %                 thetav_mu = vertcat(thetav_mu,NaN);
                			    U_mu(x,y,t,c) =  NaN ;
                                V_mu(x,y,t,c) =  NaN ;
                                theta_mu(x,y,t,c) = NaN ;
                                thetav_mu(x,y,t,c) = NaN ;
                                thetae_mu(x,y,t,c) = NaN; %jnm_new_23mar2022
                                rvap_mu(x,y,t,c) =  NaN;


                            end

                        end   %end of MU's


                        %ok to here




                        %        % yymmdd(sonday,qwq) = YYMMDD;
                        %        % ymd = str2num(yymmdd);
                        %
                        %
                        %
                        %
                        %
                        %
                        %
                        %        % hhh1 = HH(t);
                        %        % mmm1 = MM(t);
                        %        % sss1 = SS(t);
                        %        % hhh(sonday,qwq) = sprintf( '%02d', hhh1 ) ;
                        %        % mmm(sonday,qwq) = sprintf( '%02d', mmm1 ) ;
                        %        % sss(sonday,qwq) = sprintf( '%02d', sss1 ) ;
                        %
                        %        hhh(sonday,qwq) = str2num(sprintf( '%02d', HH(t) )) ;
                        %        mmm(sonday,qwq) = str2num(sprintf( '%02d', MM(t) )) ;
                        %        sss(sonday,qwq) = str2num(sprintf( '%02d', SS(t) )) ;
                        %        YYY(sonday,qwq) = str2num(sprintf( '%02d', sonde_YY)) ;
                        %        MMM(sonday,qwq) = str2num(sprintf( '%02d', sonde_MON)) ;
                        %        DDD(sonday,qwq) = str2num(sprintf( '%02d', sonde_DD)) ;


                        %       % hrhrminminsecsec(sonday,qwq) = horzcat(hhh1,mmm1,sss1);

                        %   end  % time loop


                        %matsav = horzcat('Adam_sondes_i',YYMMDD,'.mat');
                        %outmatdir = '/Users/marq789/Downloads/testskewT/';
                        %save(horzcat(outmatdir,matsav))


                    else %   if the profiles come in as NaNs

                        %		    check = 7;

                        NANdiag(x,y,t,c) = 1;


                        pCAPE_sfc(x,y,t,c)   = NaN;
                        %pCIN_NA_sfc(x,y,t) = NaN;
                        pCAPE_ml(x,y,t,c)    = NaN;
                        %pCIN_NA_ml(x,y,t)  = NaN;
                        pCAPE_mu(x,y,t,c)    = NaN;
                        %pCIN_NA_mu(x,y,t)  = NaN;


                        CAPE_mu(x,y,t,c) =  NaN;
                        %                     CIN_IB_mu(x,y,t) =  NaN;
                        LCL_height_mu(x,y,t,c) =  NaN;
                        LFC_height_mu(x,y,t,c) =  NaN;
                        EL_height_mu(x,y,t,c) =  NaN;
                        %                     LCL_pres_mu(x,y,t) =  NaN;
                        %                     LFC_pres_mu(x,y,t) =  NaN;
                        %                     EL_pres_mu(x,y,t) =  NaN;
                        %                    CIN_NA_fract_mu(x,y,t) =  NaN;
                        CIN_mu(x,y,t,c) =  NaN;
                        CAPEacbl_mu(x,y,t,c) =  NaN;
                        %                     CAPElcl_IB_mu(x,y,t) =  NaN;
                        initial_ht_parcel_mu(x,y,t,c) =  NaN;
                        %                     tallenough_mu(x,y,t) =  NaN;


                        CAPE_sfc(x,y,t,c) =  NaN;
                        %                     CIN_IB_sfc(x,y,t) =  NaN;
                        LCL_height_sfc(x,y,t,c) =  NaN;
                        LFC_height_sfc(x,y,t,c) =  NaN;
                        EL_height_sfc(x,y,t,c) =  NaN;
                        %                     LCL_pres_sfc(x,y,t) =  NaN;
                        %                     LFC_pres_sfc(x,y,t) =  NaN;
                        %                     EL_pres_sfc(x,y,t) =  NaN;
                        %                    CIN_NA_fract_sfc(x,y,t) =  NaN;
                        CIN_sfc(x,y,t,c) =  NaN;
                        CAPEacbl_sfc(x,y,t,c) =  NaN;
                        %                     CAPElcl_IB_sfc(x,y,t) =  NaN;
                        %                     tallenough_sfc(x,y,t) =  NaN;

                        CAPE_ml(x,y,t,c) =  NaN;
                        %                     CIN_IB_ml(x,y,t) =  NaN;
                        LCL_height_ml(x,y,t,c) =  NaN;
                        LFC_height_ml(x,y,t,c) =  NaN;
                        EL_height_ml(x,y,t,c) =  NaN;
                        %                     LCL_pres_ml(x,y,t) =  NaN;
                        %                     LFC_pres_ml(x,y,t) =  NaN;
                        %                     EL_pres_ml(x,y,t) =  NaN;
                        %                    CIN_NA_fract_ml(x,y,t) =  NaN;
                        CIN_ml(x,y,t,c) =  NaN;
                        CAPEacbl_ml(x,y,t,c) =  NaN;
                        %                     CAPElcl_IB_ml(x,y,t) =  NaN;
                        %                     tallenough_ml(x,y,t) =  NaN;


                        %                     EIL10_top_height(x,y,t) =  NaN;
                        %                     EIL10_bot_height(x,y,t) =  NaN;
                        %                     EIL25_top_height(x,y,t) =  NaN;
                        %                     EIL25_bot_height(x,y,t) =  NaN;
                        %                     EIL50_top_height(x,y,t) =  NaN;
                        %                     EIL50_bot_height(x,y,t) =  NaN;

                        RH_mean_ACBL(x,y,t,c) =  NaN;
                        shear_mag_bulk_ACBL_mu(x,y,t,c) =  NaN;
                        Ucrel_mu(x,y,t,c) =  NaN;
                        Vcrel_mu(x,y,t,c) =  NaN;
                        thetae_mu(x,y,t,c) =  NaN;
                        thetae_mean_subcloud_mu(x,y,t,c) =  NaN;
                        rvap_mu(x,y,t,c) =  NaN;
                        thetae_mean_ACBL_mu(x,y,t,c) =  NaN;
                        %                     LCL_temp_mu(x,y,t) =  NaN;
                        %                     LFC_temp_mu(x,y,t) =  NaN;
                        %                     EL_temp_mu(x,y,t) =  NaN;
                        shear_mag_bulk_FT_mu(x,y,t,c) =  NaN;
                        shear_dir_bulk_FT_mu(x,y,t,c) =  NaN;
                        shear_mag_bulk_0to1km(x,y,t,c) =  NaN;
                        shear_dir_bulk_0to1km(x,y,t,c) =  NaN;
                        shear_mag_bulk_0to3km(x,y,t,c) =  NaN;
                        shear_dir_bulk_0to3km(x,y,t,c) =  NaN;
                        shear_mag_bulk_0to6km(x,y,t,c) =  NaN;
                        shear_dir_bulk_0to6km(x,y,t,c) =  NaN;
                        rvap_850mb(x,y,t,c) =  NaN;
                        rvap_700mb(x,y,t,c) =  NaN;
                        rvap_500mb(x,y,t,c) =  NaN;
                        rh_850mb(x,y,t,c) =  NaN;
                        rh_700mb(x,y,t,c) =  NaN;
                        rh_500mb(x,y,t,c) =  NaN;

                        %enoch vars
                        rh_at_muLCLplus1500m(x,y,t,c) =  NaN;
                        rh_at_muLFCplus1500m(x,y,t,c) =  NaN;

                        %adam vars
                        rh_at_pmuLFCplus1500m(x,y,t,c) =  NaN;
                        RH_mean_ACBL_pmu(x,y,t,c) =  NaN;
                        shear_mag_bulk_ACBL_pmu(x,y,t,c) =  NaN;
                        thetae_mean_ACBL_pmu(x,y,t,c) =  NaN;
                        Ucrel_mean_ACBL_pmu(x,y,t,c) =  NaN;
                        Vcrel_mean_ACBL_pmu(x,y,t,c) =  NaN;

                        rvap_600mb(x,y,t,c) =  NaN;
                        rvap_400mb(x,y,t,c) =  NaN;
                        rvap_300mb(x,y,t,c) =  NaN;
                        rh_600mb(x,y,t,c) =  NaN;
                        rh_400mb(x,y,t,c) =  NaN;
                        rh_300mb(x,y,t,c) =  NaN;
                        rvap_3000masl(x,y,t,c) =  NaN;
                        rvap_4000masl(x,y,t,c) =  NaN;
                        rvap_5000masl(x,y,t,c) =  NaN;
                        rvap_6000masl(x,y,t,c) =  NaN;
                        rh_3000masl(x,y,t,c) =  NaN;
                        rh_4000masl(x,y,t,c) =  NaN;
                        rh_5000masl(x,y,t,c) =  NaN;
                        rh_6000masl(x,y,t,c) =  NaN;



                        %                     rvap_mean_EIL50(x,y,t) =  NaN;
                        %                     rvap_mean_EIL25(x,y,t) =  NaN;
                        %                     rvap_mean_EIL10(x,y,t) =  NaN;

                        %                     Ucrel_mean_ACBL_ml(x,y,t) = NaN;
                        %                     Vcrel_mean_ACBL_ml(x,y,t) = NaN;
                        %                     Ucrel_mean_0toACBL_ml(x,y,t) =  NaN;
                        %                     Vcrel_mean_0toACBL_ml(x,y,t) =  NaN;

                        Ucrel_mean_ACBL_mu(x,y,t,c) = NaN;
                        Vcrel_mean_ACBL_mu(x,y,t,c) = NaN;
                        %Ucrel_mean_0toACBL_mu(x,y,t) =  NaN;
                        %Vcrel_mean_0toACBL_mu(x,y,t) =  NaN;

                        %                     Ucrel_mean_ACBL_sfc(x,y,t) = NaN;
                        %                     Vcrel_mean_ACBL_sfc(x,y,t) = NaN;
                        %                     Ucrel_mean_0toACBL_sfc(x,y,t) =  NaN;
                        %                     Vcrel_mean_0toACBL_sfc(x,y,t) =  NaN;

                        %                     Ucrel_mean_EIL10(x,y,t) =  NaN;
                        %                     Vcrel_mean_EIL10(x,y,t) =  NaN;
                        %                     Ucrel_mean_EIL25(x,y,t) =  NaN;
                        %                     Vcrel_mean_EIL25(x,y,t) =  NaN;
                        %                     Ucrel_mean_EIL50(x,y,t) =  NaN;
                        %                     Vcrel_mean_EIL50(x,y,t) =  NaN;

                        shear_mag_bulk_0to9km(x,y,t,c) = NaN;
                        shear_dir_bulk_0to9km(x,y,t,c) = NaN;

                        %                     U_4000m(x,y,t) = NaN;
                        %                     V_4000m(x,y,t) = NaN;

                        rvap_925mb(x,y,t,c) = NaN;
                        rh_925mb(x,y,t,c) = NaN;
                        rvap_600mb(x,y,t,c) = NaN;
                        rh_600mb(x,y,t,c) = NaN;


             		    tallenough_mu(x,y,t,c)    =  NaN;
             		    pCIN_mu(x,y,t,c)       =  NaN;
             		    pCAPEacbl_mu(x,y,t,c)     =  NaN;
             		    pLFC_height_mu(x,y,t,c)   =  NaN;
             		    pEL_height_mu(x,y,t,c)    =  NaN;
             		    ptallenough_mu(x,y,t,c)   =  NaN;


                        tallenough_sfc(x,y,t,c)    =  NaN;
                        pCIN_sfc(x,y,t,c)       =  NaN;
                        pCAPEacbl_sfc(x,y,t,c)     =  NaN;
                        pLFC_height_sfc(x,y,t,c)   =  NaN;
                        pEL_height_sfc(x,y,t,c)    =  NaN;
                        ptallenough_sfc(x,y,t,c)   =  NaN;


                        tallenough_ml(x,y,t,c)    =  NaN;
                        pCIN_ml(x,y,t,c)       =  NaN;
                        pCAPEacbl_ml(x,y,t,c)     =  NaN;
                        pLFC_height_ml(x,y,t,c)   =  NaN;
                        pEL_height_ml(x,y,t,c)    =  NaN;
                        ptallenough_ml(x,y,t,c)   =  NaN;

                        U_mu(x,y,t,c) = NaN;
                        V_mu(x,y,t,c) = NaN;
                        theta_mu(x,y,t,c) = NaN;
             		    thetav_mu(x,y,t,c) = NaN;

             		    EL_temp_mu(x,y,t,c) = NaN;
             		    LCL_temp_mu(x,y,t,c) = NaN;
             		    LFC_temp_mu(x,y,t,c) = NaN;
             	 	    pLFC_temp_mu(x,y,t,c) = NaN;
              		    pEL_temp_mu(x,y,t,c) = NaN;

             		    EL_temp_ml(x,y,t,c) = NaN;
              		    LCL_temp_ml(x,y,t,c) = NaN;
             		    LFC_temp_ml(x,y,t,c) = NaN;
             		    pLFC_temp_ml(x,y,t,c) = NaN;
             		    pEL_temp_ml(x,y,t,c) = NaN;

             		    EL_temp_sfc(x,y,t,c) = NaN;
             		    LCL_temp_sfc(x,y,t,c) = NaN;
                        LFC_temp_sfc(x,y,t,c) = NaN;
             		    pLFC_temp_sfc(x,y,t,c) = NaN;
             		    pEL_temp_sfc(x,y,t,c) = NaN;


			upslope_max(x,y,t,c) = NaN;
                        upslope_mean(x,y,t,c) = NaN;
                        upslope_fract(x,y,t,c) = NaN;			 
                        upslope_depth(x,y,t,c) = NaN;

                    end

                    % ok to here


                end % x
            end % y
        end % time


    end  %end of cell loop



        %%%%%%%%%% write out netcdf fields:
        if(netcdf_out == 1)

            start_cdfwrite_time = strcat(CELLID,' -- ',datestr(now))

            %%%to make individual cell output files
            ncfile = strcat(outdir,'WIDlasso_envparms_',splidirs(9), '_',widtime,'_', splidirs(10), '_', splidirs(11), '_', splidirs(12),'.nc');
            delete(ncfile);

            parcdfer_cellstogether(ncfile, NANdiag, CAPE_mu, LCL_height_mu ,LFC_height_mu, EL_height_mu, ...
                CIN_mu, CAPEacbl_mu, initial_ht_parcel_mu, ...
                CAPE_sfc, LCL_height_sfc, LFC_height_sfc, EL_height_sfc,...
                CIN_sfc, CAPEacbl_sfc, ...
                CAPE_ml, LCL_height_ml, LFC_height_ml, EL_height_ml, CIN_ml, CAPEacbl_ml, ...
                RH_mean_ACBL, shear_mag_bulk_ACBL_mu, Ucrel_mu, Vcrel_mu, thetae_mu, thetae_mean_subcloud_mu, ...
                rvap_mu, thetae_mean_ACBL_mu, shear_mag_bulk_FT_mu, shear_dir_bulk_FT_mu, ...
                shear_mag_bulk_0to1km, shear_dir_bulk_0to1km, ...
                shear_mag_bulk_0to3km, shear_dir_bulk_0to3km, shear_mag_bulk_0to6km, shear_dir_bulk_0to6km, ...
                rvap_850mb, rvap_700mb, rvap_500mb, rh_850mb, rh_700mb, rh_500mb, freezing_ht, PW,...
                Ucrel_mean_ACBL_mu, Vcrel_mean_ACBL_mu,...
                shear_mag_bulk_0to9km, shear_dir_bulk_0to9km,...
                rvap_925mb, rh_925mb, rvap_600mb, rh_600mb,...
                rh_at_muLCLplus1500m,rh_at_muLFCplus1500m,...
                pCAPE_sfc, pCAPE_ml, pCAPE_mu, ...
                tallenough_mu, pCIN_mu, pCAPEacbl_mu, pLFC_height_mu, pEL_height_mu, ptallenough_mu,...
                tallenough_sfc, pCIN_sfc, pCAPEacbl_sfc, pLFC_height_sfc, pEL_height_sfc, ptallenough_sfc,...
                tallenough_ml, pCIN_ml, pCAPEacbl_ml, pLFC_height_ml, pEL_height_ml, ptallenough_ml,...
                U_mu, V_mu, theta_mu, thetav_mu,...
                EL_temp_sfc, LCL_temp_sfc, LFC_temp_sfc, pLFC_temp_sfc, pEL_temp_sfc,...
                EL_temp_mu, LCL_temp_mu, LFC_temp_mu, pLFC_temp_mu, pEL_temp_mu,...
                EL_temp_ml, LCL_temp_ml, LFC_temp_ml, pLFC_temp_ml, pEL_temp_ml,...
                rh_at_pmuLFCplus1500m, RH_mean_ACBL_pmu, shear_mag_bulk_ACBL_pmu, ...
                thetae_mean_ACBL_pmu, Ucrel_mean_ACBL_pmu, Vcrel_mean_ACBL_pmu,...
                rvap_400mb, rvap_300mb, rh_400mb, rh_300mb,...
                rvap_3000masl, rvap_4000masl, rvap_5000masl, rvap_6000masl, ...
                rh_3000masl, rh_4000masl, rh_5000masl, rh_6000masl,...
		upslope_max, upslope_mean, upslope_depth, upslope_fract, ...
                cell_motion_x, cell_motion_y, cellnum);


            done_cdfwrite_time = strcat(CELLID,' -- ',datestr(now))

        end










        % % simple done log:
        %donefile = horzcat(outdir,'DONE',CELLID,'.txt');
        %finished = 1;
        %writematrix(finished,donefile)


        %  check = 9;






        % % save a checkpoint mat file:
        if(savemats == 1)

            matfile = horzcat(outdir,'ZAwrf2_envparms_cell',CELLID,'.mat');


            parsaver_ZAwrf2(matfile, NANdiag, CAPE_mu, LCL_height_mu ,LFC_height_mu, EL_height_mu, ...
                CIN_mu, CAPEacbl_mu, initial_ht_parcel_mu, ...
                CAPE_sfc, LCL_height_sfc, LFC_height_sfc, EL_height_sfc,...
                CIN_sfc, CAPEacbl_sfc, ...
                CAPE_ml, LCL_height_ml, LFC_height_ml, EL_height_ml, CIN_ml, CAPEacbl_ml, ...
                RH_mean_ACBL, shear_mag_bulk_ACBL_mu, Ucrel_mu, Vcrel_mu, thetae_mu, thetae_mean_subcloud_mu, ...
                rvap_mu, thetae_mean_ACBL_mu, shear_mag_bulk_FT_mu, shear_dir_bulk_FT_mu, ...
                shear_mag_bulk_0to1km, shear_dir_bulk_0to1km, ...
                shear_mag_bulk_0to3km, shear_dir_bulk_0to3km, shear_mag_bulk_0to6km, shear_dir_bulk_0to6km, ...
                rvap_850mb, rvap_700mb, rvap_500mb, rh_850mb, rh_700mb, rh_500mb, freezing_ht, PW,...
                Ucrel_mean_ACBL_mu, Vcrel_mean_ACBL_mu,...
                shear_mag_bulk_0to9km, shear_dir_bulk_0to9km,...
                rvap_925mb, rh_925mb, rvap_600mb, rh_600mb,...
                rh_at_muLCLplus1500m,rh_at_muLFCplus1500m,...
                pCAPE_sfc, pCAPE_ml, pCAPE_mu, ...
                tallenough_mu, pCIN_mu, pCAPEacbl_mu, pLFC_height_mu, pEL_height_mu, ptallenough_mu,...
                tallenough_sfc, pCIN_sfc, pCAPEacbl_sfc, pLFC_height_sfc, pEL_height_sfc, ptallenough_sfc,...
                tallenough_ml, pCIN_ml, pCAPEacbl_ml, pLFC_height_ml, pEL_height_ml, ptallenough_ml,...
                U_mu, V_mu, theta_mu, thetav_mu,...
                EL_temp_sfc, LCL_temp_sfc, LFC_temp_sfc, pLFC_temp_sfc, pEL_temp_sfc,...
                EL_temp_mu, LCL_temp_mu, LFC_temp_mu, pLFC_temp_mu, pEL_temp_mu,...
                EL_temp_ml, LCL_temp_ml, LFC_temp_ml, pLFC_temp_ml, pEL_temp_ml,...
                rh_at_pmuLFCplus1500m, RH_mean_ACBL_pmu, shear_mag_bulk_ACBL_pmu, ...
                thetae_mean_ACBL_pmu, Ucrel_mean_ACBL_pmu, Vcrel_mean_ACBL_pmu,...
                rvap_400mb, rvap_300mb, rh_400mb, rh_300mb,...
                rvap_3000masl, rvap_4000masl, rvap_5000masl, rvap_6000masl, ...
                rh_3000masl, rh_4000masl, rh_5000masl, rh_6000masl,...
                cell_motion_x, cell_motion_y, cellnum);
        end

        start_calculate_time = strcat(CELLID,' --' ,start_calculate_time)
        done_calculate_time = strcat(CELLID,' -- ', datestr(now)  )

        labc = strcat('done cell ',num2str(c) );
        disp(labc)

    %end % parforcell track loop

end   %time files
