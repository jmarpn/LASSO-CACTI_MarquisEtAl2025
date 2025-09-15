

clear all


maxnum_guess = 4000;
depth = 38;
X = 500;
Y = 1600;

KEEP_num_updrafts = 1200;

% run = '20181129_gefs09_base_d4';

% run = '20190122_gefs18_base_d4';

% run = '20190125_eda07_base_d4';

% run = '20190129_gefs11_base_d4';

% RUN = ['20181129_gefs09_base_d4';
%        '20190122_gefs18_base_d4';
%        '20190125_eda07_base__d4';
%        '20190129_gefs11_base_d4'];

RUN = ['20181129_gefs09_base_d4';
        '20181204_gefs19_base_d4';
        '20190122_gefs18_base_d4';
        '20190123_gefs18_base_d4';
        '20190125_eda07_base__d4';
        '20190129_gefs11_base_d4';
        '20190208_eda08__base_d4'];

for R = 1:7

   % R = 2;

    run = RUN(R,:)

    %post-processed wID mat files:
    sourcedir = ['/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/',run,'/'] ;
    wIDtimes_list = ls( horzcat(sourcedir,'/wID_ObjectOutput_NSj_100m5min_qcldice0.0001_dbz10_tt*mat') )  ;
%     if(R==2)
%         wIDtimes_list = ls( horzcat(sourcedir,'/wID_ObjectOutput_100m5min_qcld0.0001_dbz10_tt*mat') )  ;        
%     else
%         wIDtimes_list = ls( horzcat(sourcedir,'/TEST_wID_ObjectOutput_100m5min_qcld0.0001_dbz10_tt*mat') )  ;
%     end

    wIDtimes_list = split(wIDtimes_list) ;
    wIDtimes_list(end) = []  ;
    [sa sb] = size(wIDtimes_list); clear sb;

    FwIDprofs_equivarea     = zeros(maxnum_guess,depth,sa);    FwIDprofs_equivarea(:) = NaN;
    FwIDprofs_numij         = zeros(maxnum_guess,depth,sa);    FwIDprofs_numij(:) = NaN;
    FwIDprofs_maxmag        = zeros(maxnum_guess,depth,sa);    FwIDprofs_maxmag(:) = NaN;
    FwIDprofs_massflux      = zeros(maxnum_guess,depth,sa);    FwIDprofs_massflux(:) = NaN;
    FwIDprofs_centlat       = zeros(maxnum_guess,depth,sa);    FwIDprofs_centlat(:) = NaN;
    FwIDprofs_centlon       = zeros(maxnum_guess,depth,sa);    FwIDprofs_centlon(:) = NaN;
    FwIDprofs_centHTagl     = zeros(maxnum_guess,depth,sa);    FwIDprofs_centHTagl(:) = NaN;
    FwIDprofs_centHTasl     = zeros(maxnum_guess,depth,sa);    FwIDprofs_centHTasl(:) = NaN;
    FwIDprofs_northj        = zeros(maxnum_guess,depth,sa);    FwIDprofs_northj(:) = NaN;
    FwIDprofs_southj        = zeros(maxnum_guess,depth,sa);    FwIDprofs_southj(:) = NaN;    
    FwIDprofs_ID            = zeros(maxnum_guess,sa);          FwIDprofs_ID(:) = NaN;
    FwIDprofs_klfcagl       = zeros(maxnum_guess,sa);          FwIDprofs_klfcagl(:) = NaN;

    HHMMSS = [];
    YYMMDD = [];




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% sounding environments:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    FwID_MUCAPE     = zeros(maxnum_guess,sa);               FwID_MUCAPE(:) = NaN;
    FwID_MUCIN     = zeros(maxnum_guess,sa);                FwID_MUCIN(:) = NaN;
    FwID_MULFC     = zeros(maxnum_guess,sa);                FwID_MULFC(:) = NaN;
    FwID_MULCL     = zeros(maxnum_guess,sa);                FwID_MULCL(:) = NaN;
    FwID_MUEL     = zeros(maxnum_guess,sa);                 FwID_MUEL(:) = NaN;
    FwID_MUvap     = zeros(maxnum_guess,sa);                FwID_MUvap(:) = NaN;
    FwID_MUtheta     = zeros(maxnum_guess,sa);              FwID_MUtheta(:) = NaN;
    FwID_MUthetae     = zeros(maxnum_guess,sa);             FwID_MUthetae(:) = NaN;
    FwID_MUACBLCAPE     = zeros(maxnum_guess,sa);           FwID_MUACBLCAPE(:) = NaN;
    FwID_RH_mean_ACBL = zeros(maxnum_guess,sa);             FwID_RH_mean_ACBL(:) = NaN;
    FwID_RH_500mb = zeros(maxnum_guess,sa);                 FwID_RH_500mb(:) = NaN;
    FwID_RH_600mb = zeros(maxnum_guess,sa);                 FwID_RH_600mb(:) = NaN;
    FwID_RH_5kmasl = zeros(maxnum_guess,sa);                FwID_RH_5kmasl(:) = NaN;
    FwID_RH_6kmasl = zeros(maxnum_guess,sa);                FwID_RH_6kmasl(:) = NaN;
    FwID_thetae_mean_ACBL_mu = zeros(maxnum_guess,sa);      FwID_thetae_mean_ACBL_mu(:) = NaN;
    FwID_shear_mag_bulk_FT_mu = zeros(maxnum_guess,sa);     FwID_shear_mag_bulk_FT_mu(:) = NaN;
    FwID_shear_mag_bulk_ACBL_mu = zeros(maxnum_guess,sa);   FwID_shear_mag_bulk_ACBL_mu(:) = NaN;
    FwID_shear_mag_bulk_0to1km = zeros(maxnum_guess,sa);    FwID_shear_mag_bulk_0to1km(:) = NaN;
    FwID_shear_mag_bulk_0to3km = zeros(maxnum_guess,sa);    FwID_shear_mag_bulk_0to3km(:) = NaN;
    FwID_shear_mag_bulk_0to6km = zeros(maxnum_guess,sa);    FwID_shear_mag_bulk_0to6km(:) = NaN;
    FwID_shear_mag_bulk_0to9km = zeros(maxnum_guess,sa);    FwID_shear_mag_bulk_0to9km(:) = NaN;
    FwID_upslope_mean = zeros(maxnum_guess,sa);             FwID_upslope_mean(:) = NaN;
    FwID_upslope_max = zeros(maxnum_guess,sa);              FwID_upslope_max(:) = NaN;
    FwID_upslope_depth = zeros(maxnum_guess,sa);            FwID_upslope_depth(:) = NaN;
    FwID_upslope_fract = zeros(maxnum_guess,sa);            FwID_upslope_fract(:) = NaN;


    %esourcedir = ['/Users/marq789/Downloads/wenvsave/'] ;
    esourcedir = ['/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/',run,'/'] ;
    %wIDenvs_list = ls( horzcat(esourcedir,'/WIDlasso_envparms_*nc') )  ;
    %wIDenvs_list = ls( horzcat(esourcedir,'/WIDlasso_deeper_envparms_*nc') )  ;
    wIDenvs_list = ls( horzcat(esourcedir,'/WIDlasso_deeper_3R_envparms_*nc') )  ;

    wIDenvs_list = split(wIDenvs_list) ;
    wIDenvs_list(end) = []  ;
    [ea eb] = size(wIDenvs_list); clear eb;

        %  ncdisp(char(wIDenvs_list(1) ))


    for t = 1 : sa

        %  t = 8;

        disp(['t=',num2str(t)])
        load(  char(wIDtimes_list(t,:) )  )

        [mt mz] = size(wIDprofs_equivarea);
        %[a1 a2 mt] = size(wIDprofs_2d_masks);
        clear a1 a2
        disp(['mt=',num2str(mt)])

        HHMMSS = vertcat(HHMMSS,HHMMSS_all);
        YYMMDD = vertcat(YYMMDD,YYMMDD_all);

        % assign wID profiles
        FwIDprofs_equivarea(1:mt,1:depth,t)     =   wIDprofs_equivarea;
        FwIDprofs_numij(1:mt,1:depth,t)         =   wIDprofs_numij;
        FwIDprofs_maxmag(1:mt,1:depth,t)        =   wIDprofs_maxmag;
        FwIDprofs_massflux(1:mt,1:depth,t)      =   wIDprofs_massflux;
        FwIDprofs_centlat(1:mt,1:depth,t)       =   wIDprofs_centlat;
        FwIDprofs_centlon(1:mt,1:depth,t)       =   wIDprofs_centlon;
        FwIDprofs_centHTagl(1:mt,1:depth,t)     =   wIDprofs_centHTagl;
        FwIDprofs_centHTasl(1:mt,1:depth,t)     =   wIDprofs_centHTasl;
        FwIDprofs_northj(1:mt,1:depth,t)        =   wIDprofs_northj;
        FwIDprofs_southj(1:mt,1:depth,t)        =   wIDprofs_southj;
        FwIDprofs_ID(1:mt,t)                    =   wIDprofs_ID;
        FwIDprofs_klfcagl(1:mt,t)               =   wIDprofs_klfcagl;
        %FwIDprofs_2d_masks(1:X,1:Y,1:mt,t)      =   wIDprofs_2d_masks;
        FwIDprofs_zagl                          =   wIDprofs_zagl;

        % ncdisp(string(wIDenvs_list(t)))

        % read in & assign wID environments, assigning each the mean around the circle and excluding element 1 (the center pt of the updraft because of contamination):
        var = ncread(string(wIDenvs_list(t)),'CAPE_mu');
        FwID_MUCAPE(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);
        %FwID_MUCAPE(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'CAPE_mu'), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'CIN_mu');
        %FwID_MUCIN(1:mt,t)               =  permute( mean( ncread(string(wIDenvs_list(t)),'CIN_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MUCIN(1:mt,t)               =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'LFC_height_mu');
        %FwID_MULFC(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'LFC_height_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MULFC(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'LCL_height_mu');        
        %FwID_MULCL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'LCL_height_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MULCL(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'EL_height_mu');
        %FwID_MUEL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'EL_height_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MUEL(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'rvap_mu');
        %FwID_MUvap(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'rvap_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MUvap(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'theta_mu');        
        %FwID_MUtheta(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'theta_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MUtheta(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'thetae_mu');        
        %FwID_MUthetae(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'thetae_mu'), 1 ,'omitnan' ), [4 1 2 3]); 
        FwID_MUthetae(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'CAPEacbl_mu');
        %FwID_MUACBLCAPE(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'CAPEacbl_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_MUACBLCAPE(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'RH_mean_ACBL');    
        %FwID_RH_mean_ACBL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'RH_mean_ACBL'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_RH_mean_ACBL(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'rh_500mb');    
        %FwID_RH_mean_ACBL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'RH_mean_ACBL'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_RH_500mb(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'rh_600mb');    
        %FwID_RH_mean_ACBL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'RH_mean_ACBL'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_RH_600mb(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'rh_5000masl');    
        %FwID_RH_mean_ACBL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'RH_mean_ACBL'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_RH_5kmasl(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'rh_6000masl');    
        %FwID_RH_mean_ACBL(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'RH_mean_ACBL'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_RH_6kmasl(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'thetae_mean_ACBL_mu');
        %FwID_thetae_mean_ACBL_mu(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'thetae_mean_ACBL_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_thetae_mean_ACBL_mu(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'shear_mag_bulk_FT_mu');        
        %FwID_shear_mag_bulk_FT_mu(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'shear_mag_bulk_FT_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_shear_mag_bulk_FT_mu(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'shear_mag_bulk_ACBL_mu');         
        %FwID_shear_mag_bulk_ACBL_mu(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'shear_mag_bulk_ACBL_mu'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_shear_mag_bulk_ACBL_mu(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to1km');         
        %FwID_shear_mag_bulk_0to1km(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to1km'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_shear_mag_bulk_0to1km(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to3km'); 
        %FwID_shear_mag_bulk_0to3km(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to3km'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_shear_mag_bulk_0to3km(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to6km'); 
        %FwID_shear_mag_bulk_0to6km(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to6km'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_shear_mag_bulk_0to6km(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to9km'); 
        %FwID_shear_mag_bulk_0to9km(1:mt,t)              =  permute( mean( ncread(string(wIDenvs_list(t)),'shear_mag_bulk_0to9km'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_shear_mag_bulk_0to9km(1:mt,t)              =  permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'upslope_mean'); 
        %FwID_upslope_mean(1:mt,t)               =    permute( mean( ncread(string(wIDenvs_list(t)),'upslope_mean'), 1 ,'omitnan' ), [4 1 2 3]);
        FwID_upslope_mean(1:mt,t)               =    permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'upslope_max');         
        %FwID_upslope_max(1:mt,t)                =    permute( mean( ncread(string(wIDenvs_list(t)),'upslope_max'), 1 ,'omitnan' ), [4 1 2 3]); 
        FwID_upslope_max(1:mt,t)                =    permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        var = ncread(string(wIDenvs_list(t)),'upslope_depth'); 
        %FwID_upslope_depth(1:mt,t)              =    permute( mean( ncread(string(wIDenvs_list(t)),'upslope_depth'), 1 ,'omitnan' ), [4 1 2 3]); 
        FwID_upslope_depth(1:mt,t)              =    permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]); 

        var = ncread(string(wIDenvs_list(t)),'upslope_fract');         
        %FwID_upslope_fract(1:mt,t)              =    permute( mean( ncread(string(wIDenvs_list(t)),'upslope_fract'), 1 ,'omitnan' ), [4 1 2 3]);  
        FwID_upslope_fract(1:mt,t)              =    permute( mean( var(2:end,:,:,:), 1 ,'omitnan' ), [4 1 2 3]);

        clear wIDprofs_equivarea wIDprofs_numij wIDprofs_maxmag wIDprofs_massflux wIDprofs_centlat wIDprofs_centlon
        clear wIDprofs_centHTagl wIDprofs_centHTasl wIDprofs_ID wIDprofs_klfcagl wIDprofs_zagl wIDprofs_northj wIDprofs_southj

    end

    %purge the extra updraft dummy columns
    FwIDprofs_equivarea(KEEP_num_updrafts+1:end,:,:)     =   [];
    FwIDprofs_numij(KEEP_num_updrafts+1:end,:,:)         =   [];
    FwIDprofs_maxmag(KEEP_num_updrafts+1:end,:,:)        =   [];
    FwIDprofs_massflux(KEEP_num_updrafts+1:end,:,:)      =   [];
    FwIDprofs_centlat(KEEP_num_updrafts+1:end,:,:)       =   [];
    FwIDprofs_centlon(KEEP_num_updrafts+1:end,:,:)       =   [];
    FwIDprofs_centHTagl(KEEP_num_updrafts+1:end,:,:)     =   [];
    FwIDprofs_centHTasl(KEEP_num_updrafts+1:end,:,:)     =   [];
    FwIDprofs_ID(KEEP_num_updrafts+1:end,:)              =   [];
    FwIDprofs_klfcagl(KEEP_num_updrafts+1:end,:)         =   [];
    FwIDprofs_northj(KEEP_num_updrafts+1:end,:,:)        =   [];
    FwIDprofs_southj(KEEP_num_updrafts+1:end,:,:)        =   [];

    FwID_MUCAPE(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MUCIN(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MULFC(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MULCL(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MUEL(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MUvap(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MUtheta(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MUthetae(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_MUACBLCAPE(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_RH_mean_ACBL(KEEP_num_updrafts+1:end,:)     =   [];

    FwID_RH_500mb(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_RH_600mb(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_RH_5kmasl(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_RH_6kmasl(KEEP_num_updrafts+1:end,:)     =   [];

    FwID_thetae_mean_ACBL_mu(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_shear_mag_bulk_FT_mu(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_shear_mag_bulk_ACBL_mu(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_shear_mag_bulk_0to1km(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_shear_mag_bulk_0to3km(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_shear_mag_bulk_0to6km(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_shear_mag_bulk_0to9km(KEEP_num_updrafts+1:end,:)     =   [];

    FwID_upslope_mean(KEEP_num_updrafts+1:end,:)      =   [];
    FwID_upslope_max(KEEP_num_updrafts+1:end,:)       =   [];
    FwID_upslope_depth(KEEP_num_updrafts+1:end,:)     =   [];
    FwID_upslope_fract(KEEP_num_updrafts+1:end,:)     =   [];


    outfile = strcat(sourcedir,'wID_merged_deepmetB3R_newtimes_',run,'.mat');
    %outfile = strcat(sourcedir,'wID_merged_',run,'.mat');
    save(outfile,'HHMMSS','YYMMDD','FwIDprofs_equivarea','FwIDprofs_numij','FwIDprofs_maxmag','FwIDprofs_massflux','FwIDprofs_centlat',...
        'FwIDprofs_centlon','FwIDprofs_centHTagl','FwIDprofs_centHTasl','FwIDprofs_ID','FwIDprofs_klfcagl','FwIDprofs_zagl',...
        'FwIDprofs_northj','FwIDprofs_southj','FwID_MUCAPE' ,'FwID_MUCIN' ,'FwID_MULFC' ,'FwID_MULCL' ,'FwID_MUEL', 'FwID_MUvap','FwID_MUACBLCAPE' ,'FwID_RH_mean_ACBL',... 
        'FwID_RH_500mb', 'FwID_RH_600mb', 'FwID_RH_5kmasl', 'FwID_RH_6kmasl', ...
        'FwID_thetae_mean_ACBL_mu' ,'FwID_shear_mag_bulk_FT_mu' ,'FwID_shear_mag_bulk_ACBL_mu' ,'FwID_shear_mag_bulk_0to1km' , ...
        'FwID_shear_mag_bulk_0to3km' ,'FwID_shear_mag_bulk_0to6km','FwID_shear_mag_bulk_0to9km',...
        'FwID_upslope_mean', 'FwID_upslope_max', 'FwID_upslope_depth', 'FwID_upslope_fract', 'FwID_MUtheta', 'FwID_MUthetae' )

end

%%

clear

load('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20181129_gefs09_base_d4/wID_merged_deepmetB3R_newtimes_20181129_gefs09_base_d4.mat')


