function parcdfer_ZAwrf3(ncfile, NANdiag, CAPE_mu, LCL_height_mu ,LFC_height_mu, EL_height_mu, ...
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

        
        % surface-based lifted parcels
        
         [xl yl tl] = size(tallenough_sfc);
         nccreate(ncfile,'tallenough_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'tallenough_sfc',tallenough_sfc) ;
         ncwriteatt(ncfile,'tallenough_sfc','units','unitless');
         ncwriteatt(ncfile,'tallenough_sfc','description','0 = temperature profile not deep enough at this time to capture EL for lifted sfc parcel; 1 = profile is deep enough');
        
% %          [xl yl tl] = size(tallenough_sfc);
% %          nccreate(ncfile,'tallenough_sfc',...
% %              'Dimensions', {'x',xl,'y',yl,'t',tl,'cell',1},...
% %              'FillValue','disable');
% %          ncwrite(ncfile,'tallenough_sfc',tallenough_sfc) ;
% %          ncwriteatt(ncfile,'tallenough_sfc','units','unitless');
% %          ncwriteatt(ncfile,'tallenough_sfc','description','0 = temperature profile not deep enough at this time to capture EL for lifted sfc parcel; 1 = profile is deep enough');
        


        % ncdisp(ncfile)

         [xl yl tl] = size(ptallenough_sfc);
         nccreate(ncfile,'ptallenough_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'ptallenough_sfc',ptallenough_sfc) ;
         ncwriteatt(ncfile,'ptallenough_sfc','units','unitless');
         ncwriteatt(ncfile,'ptallenough_sfc','description','0 = temperature profile not deep enough at this time to capture EL for lifted sfc pseudoadiabatic parcel; 1 = profile is deep enough');


        nccreate(ncfile,'cellnum',...
            'Dimensions', {'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'cellnum',cellnum) ;
        ncwriteatt(ncfile,'cellnum','units','non-dim');
        ncwriteatt(ncfile,'cellnum','description','Cell number in tracking data set');

       
        nccreate(ncfile,'cell_motion_x',...
            'Dimensions', {'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'cell_motion_x',cell_motion_x) ;
        ncwriteatt(ncfile,'cell_motion_x','units','m/s');
        ncwriteatt(ncfile,'cell_motion_x','description','Zonal component of cell motion (mean during first 3 time steps of its track)');
        
        
        nccreate(ncfile,'cell_motion_y',...
            'Dimensions', {'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'cell_motion_y',cell_motion_y) ;
        ncwriteatt(ncfile,'cell_motion_y','units','m/s');
        ncwriteatt(ncfile,'cell_motion_y','description','Merdional component of cell motion (mean during first 3 time steps of its track)');
        
        
        
        [xl yl tl] = size(CAPE_sfc);
        nccreate(ncfile,'CAPE_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CAPE_sfc',CAPE_sfc) ;
        ncwriteatt(ncfile,'CAPE_sfc','units','J/kg');
        ncwriteatt(ncfile,'CAPE_sfc','description','Integration of buoyancy between the moist LFC and EL for lifted sfc parcel (reversible warm moist lift)');
        
        [xl yl tl] = size(CIN_sfc);
        nccreate(ncfile,'CIN_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CIN_sfc',CIN_sfc) ;
        ncwriteatt(ncfile,'CIN_sfc','units','J/kg');
        ncwriteatt(ncfile,'CIN_sfc','description','Integration of negative buoyancy between the sfc and moist LFC for a lifted sfc parcel (reversible warm moist lift)');
        
%         [xl yl tl] = size(CIN_NA_fract_sfc);
%         nccreate(ncfile,'CIN_NA_fract_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CIN_NA_fract_sfc',CIN_NA_fract_sfc) ;
%         ncwriteatt(ncfile,'CIN_NA_fract_sfc','units','unitless');
%         ncwriteatt(ncfile,'CIN_NA_fract_sfc','description','Fraction of the number of data points between the sfc and LFC containing negative buoyancy for a lifted sfc parcel (reversible warm moist lift)');
        
%         [xl yl tl] = size(CIN_IB_sfc);
%         nccreate(ncfile,'CIN_IB_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CIN_IB_sfc',CIN_IB_sfc) ;
%         ncwriteatt(ncfile,'CIN_IB_sfc','units','J/kg');
%         ncwriteatt(ncfile,'CIN_IB_sfc','description','Integration of total parcel buoyancy (positive + negative) between the sfc and LFC');
        
        [xl yl tl] = size(LFC_height_sfc);
        nccreate(ncfile,'LFC_height_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LFC_height_sfc',LFC_height_sfc) ;
        ncwriteatt(ncfile,'LFC_height_sfc','units','m');
        ncwriteatt(ncfile,'LFC_height_sfc','description','Height above ground level of the moist level of free convection for lifted sfc parcel (reversible warm moist lift)');
       
        [xl yl tl] = size(pLFC_height_sfc);
        nccreate(ncfile,'pLFC_height_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pLFC_height_sfc',pLFC_height_sfc) ;
        ncwriteatt(ncfile,'pLFC_height_sfc','units','m');
        ncwriteatt(ncfile,'pLFC_height_sfc','description','Height above ground level of the moist level of free convection for lifted sfc parcel (pseudoadiabatic warm moist lift)');


%         [xl yl tl] = size(LFC_pres_sfc);
%         nccreate(ncfile,'LFC_pres_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'LFC_pres_sfc',LFC_pres_sfc) ;
%         ncwriteatt(ncfile,'LFC_pres_sfc','units','hPa');
%         ncwriteatt(ncfile,'LFC_pres_sfc','description','Pressure level of the level of free convection for lifted sfc parcel');
        
         [xl yl tl] = size(LFC_temp_sfc);
         nccreate(ncfile,'LFC_temp_sfc',...
                  'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'LFC_temp_sfc',LFC_temp_sfc) ;
         ncwriteatt(ncfile,'LFC_temp_sfc','units','K');
         ncwriteatt(ncfile,'LFC_temp_sfc','description','Environmental temperature at the level of free convection for lifted sfc parcel');
        
         [xl yl tl] = size(pLFC_temp_sfc);
         nccreate(ncfile,'pLFC_temp_sfc',...
                  'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'pLFC_temp_sfc',pLFC_temp_sfc) ;
         ncwriteatt(ncfile,'pLFC_temp_sfc','units','K');
         ncwriteatt(ncfile,'pLFC_temp_sfc','description','Environmental temperature at the level of free convection for lifted sfc pseudoadiabatic parcel');


        [xl yl tl] = size(LCL_height_sfc);
        nccreate(ncfile,'LCL_height_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LCL_height_sfc',LCL_height_sfc) ;
        ncwriteatt(ncfile,'LCL_height_sfc','units','m');
        ncwriteatt(ncfile,'LCL_height_sfc','description','Height above ground of condensation level for lifted sfc parcel');
        
%         [xl yl tl] = size(LCL_pres_sfc);
%         nccreate(ncfile,'LCL_pres_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'LCL_pres_sfc',LCL_pres_sfc) ;
%         ncwriteatt(ncfile,'LCL_pres_sfc','units','hPa');
%         ncwriteatt(ncfile,'LCL_pres_sfc','description','Pressure level of condensation level for lifted sfc parcel');
        
        [xl yl tl] = size(LCL_temp_sfc);
         nccreate(ncfile,'LCL_temp_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'LCL_temp_sfc',LCL_temp_sfc) ;
         ncwriteatt(ncfile,'LCL_temp_sfc','units','K');
         ncwriteatt(ncfile,'LCL_temp_sfc','description','Environmental temperature at the condensation level for lifted sfc parcel');
        
        [xl yl tl] = size(EL_height_sfc);
        nccreate(ncfile,'EL_height_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'EL_height_sfc',EL_height_sfc) ;
        ncwriteatt(ncfile,'EL_height_sfc','units','m');
        ncwriteatt(ncfile,'EL_height_sfc','description','Height above ground of equilibrium level for lifted sfc parcel (reversible warm moist lift)');
        
        [xl yl tl] = size(pEL_height_sfc);
        nccreate(ncfile,'pEL_height_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pEL_height_sfc',pEL_height_sfc) ;
        ncwriteatt(ncfile,'pEL_height_sfc','units','m');
        ncwriteatt(ncfile,'pEL_height_sfc','description','Height above ground of equilibrium level for lifted sfc parcel (pseudoadiabatic warm moist lift)');

%         [xl yl tl] = size(EL_pres_sfc);
%         nccreate(ncfile,'EL_pres_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EL_pres_sfc',EL_pres_sfc) ;
%         ncwriteatt(ncfile,'EL_pres_sfc','units','hPa');
%         ncwriteatt(ncfile,'EL_pres_sfc','description','Pressure level of equilibrium level for lifted sfc parcel');
        
        [xl yl tl] = size(EL_temp_sfc);
        nccreate(ncfile,'EL_temp_sfc',...
                  'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
        ncwrite(ncfile,'EL_temp_sfc',EL_temp_sfc) ;
        ncwriteatt(ncfile,'EL_temp_sfc','units','K');
        ncwriteatt(ncfile,'EL_temp_sfc','description','Environmental temperature at the equilibrium level for lifted sfc parcel');
        
        
        [xl yl tl] = size(pEL_temp_sfc);
        nccreate(ncfile,'pEL_temp_sfc',...
                  'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
        ncwrite(ncfile,'pEL_temp_sfc',pEL_temp_sfc) ;
        ncwriteatt(ncfile,'pEL_temp_sfc','units','K');
        ncwriteatt(ncfile,'pEL_temp_sfc','description','Environmental temperature at the equilibrium level for pseudoadiabatic lifted sfc parcel');


        [xl yl tl] = size(CAPEacbl_sfc);
        nccreate(ncfile,'CAPEacbl_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CAPEacbl_sfc',CAPEacbl_sfc) ;
        ncwriteatt(ncfile,'CAPEacbl_sfc','units','J/kg');
        ncwriteatt(ncfile,'CAPEacbl_sfc','description','Integration of parcel buoyancy in the active cloud bearing layer (reversible warm moist lift)');
       

        [xl yl tl] = size(pCAPEacbl_sfc);
        nccreate(ncfile,'pCAPEacbl_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCAPEacbl_sfc',pCAPEacbl_sfc) ;
        ncwriteatt(ncfile,'pCAPEacbl_sfc','units','J/kg');
        ncwriteatt(ncfile,'pCAPEacbl_sfc','description','Integration of parcel buoyancy in the active cloud bearing layer (pseudoadiabatic warm moist lift)');

%         [xl yl tl] = size(CAPElcl_IB_sfc);
%         nccreate(ncfile,'CAPElcl_IB_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CAPElcl_IB_sfc',CAPElcl_IB_sfc) ;
%         ncwriteatt(ncfile,'CAPElcl_IB_sfc','units','J/kg');
%         ncwriteatt(ncfile,'CAPElcl_IB_sfc','description','Integration of total parcel buoyancy (positive + negative) between the LCL (for a lifted sfc parcel) and 2 km above it');
        

        [xl yl tl] = size(pCAPE_sfc);
        nccreate(ncfile,'pCAPE_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCAPE_sfc',pCAPE_sfc) ;
        ncwriteatt(ncfile,'pCAPE_sfc','units','J/kg');
        ncwriteatt(ncfile,'pCAPE_sfc','description','Integration of buoyancy between the moist LFC and EL for lifted sfc parcel (warm moist pseudoadiabatic lift, no entrainemnt)');
        
        [xl yl tl] = size(pCIN_sfc);
        nccreate(ncfile,'pCIN_sfc',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCIN_sfc',pCIN_sfc) ;
        ncwriteatt(ncfile,'pCIN_sfc','units','J/kg');
        ncwriteatt(ncfile,'pCIN_sfc','description','Integration of negative buoyancy between the sfc and LFC for a lifted sfc parcel (warm moist pseudoadiabatic lift, no entrainemnt)');












        
        % lifted 100-hPa mean layer parcel metrics
        
        
        [xl yl tl] = size(tallenough_ml);
        nccreate(ncfile,'tallenough_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'tallenough_ml',tallenough_ml) ;
        ncwriteatt(ncfile,'tallenough_ml','units','unitless');
        ncwriteatt(ncfile,'tallenough_ml','description','0 = temperature profile not deep enough at this time to capture EL for lifted ML parcel; 1 = profile is deep enough');
        
         [xl yl tl] = size(ptallenough_ml);
         nccreate(ncfile,'ptallenough_ml',...
              'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'ptallenough_ml',ptallenough_ml) ;
         ncwriteatt(ncfile,'ptallenough_ml','units','unitless');
         ncwriteatt(ncfile,'ptallenough_ml','description','0 = temperature profile not deep enough at this time to capture EL for lifted ml pseudoadiabatic parcel; 1 = profile is deep enough');


        [xl yl tl] = size(CAPE_ml);
        nccreate(ncfile,'CAPE_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CAPE_ml',CAPE_ml) ;
        ncwriteatt(ncfile,'CAPE_ml','units','J/kg');
        ncwriteatt(ncfile,'CAPE_ml','description','Integration of buoyancy between the LFC and EL for lifted ML parcel (reversible warm moist lift)');
        
        [xl yl tl] = size(CIN_ml);
        nccreate(ncfile,'CIN_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CIN_ml',CIN_ml) ;
        ncwriteatt(ncfile,'CIN_ml','units','J/kg');
        ncwriteatt(ncfile,'CIN_ml','description','Integration of negative buoyancy between the origin height and moist LFC for a lifted ML parcel (reversible warm moist lift)');
        
%         [xl yl tl] = size(CIN_NA_fract_ml);
%         nccreate(ncfile,'CIN_NA_fract_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CIN_NA_fract_ml',CIN_NA_fract_ml) ;
%         ncwriteatt(ncfile,'CIN_NA_fract_ml','units','unitless');
%         ncwriteatt(ncfile,'CIN_NA_fract_ml','description','Fraction of the number of data points between the origin height and LFC (for a lifted ML parcel) containing negative buoyancy (reversible warm moist lift)');
        
%         [xl yl tl] = size(CIN_IB_ml);
%         nccreate(ncfile,'CIN_IB_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CIN_IB_ml',CIN_IB_ml) ;
%         ncwriteatt(ncfile,'CIN_IB_ml','units','J/kg');
%         ncwriteatt(ncfile,'CIN_IB_ml','description','Integration of total parcel buoyancy (positive + negative) between the origin height and LFC for a lifted ML parcel');
        
        [xl yl tl] = size(LFC_height_ml);
        nccreate(ncfile,'LFC_height_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LFC_height_ml',LFC_height_ml) ;
        ncwriteatt(ncfile,'LFC_height_ml','units','m');
        ncwriteatt(ncfile,'LFC_height_ml','description','Height above ground of the moist level of free convection for lifted ML parcel (reversible warm moist lift)');
       
        [xl yl tl] = size(pLFC_height_ml);
        nccreate(ncfile,'pLFC_height_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pLFC_height_ml',pLFC_height_ml) ;
        ncwriteatt(ncfile,'pLFC_height_ml','units','m');
        ncwriteatt(ncfile,'pLFC_height_ml','description','Height above ground of the moist level of free convection for lifted ML parcel (pseudoadiabtaic warm moist lift)');


%         [xl yl tl] = size(LFC_pres_ml);
%         nccreate(ncfile,'LFC_pres_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'LFC_pres_ml',LFC_pres_ml) ;
%         ncwriteatt(ncfile,'LFC_pres_ml','units','hPa');
%         ncwriteatt(ncfile,'LFC_pres_ml','description','Pressure level of the level of free convection for lifted ML parcel');
        
        [xl yl tl] = size(LFC_temp_ml);
        nccreate(ncfile,'LFC_temp_ml',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LFC_temp_ml',LFC_temp_ml) ;
        ncwriteatt(ncfile,'LFC_temp_ml','units','K');
        ncwriteatt(ncfile,'LFC_temp_ml','description','Environmental temperature at the level of free convection for lifted ML parcel');
        
        [xl yl tl] = size(pLFC_temp_ml);
        nccreate(ncfile,'pLFC_temp_ml',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pLFC_temp_ml',pLFC_temp_ml) ;
        ncwriteatt(ncfile,'pLFC_temp_ml','units','K');
        ncwriteatt(ncfile,'pLFC_temp_ml','description','Environmental temperature at the level of free convection for pseudoadiabatic lifted ML parcel');        


        [xl yl tl] = size(LCL_height_ml);
        nccreate(ncfile,'LCL_height_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LCL_height_ml',LCL_height_ml) ;
        ncwriteatt(ncfile,'LCL_height_ml','units','m');
        ncwriteatt(ncfile,'LCL_height_ml','description','Height above ground of condensation level for lifted ML parcel');
        
%         [xl yl tl] = size(LCL_pres_ml);
%         nccreate(ncfile,'LCL_pres_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'LCL_pres_ml',LCL_pres_ml) ;
%         ncwriteatt(ncfile,'LCL_pres_ml','units','hPa');
%         ncwriteatt(ncfile,'LCL_pres_ml','description','Pressure level of condensation level for lifted ML parcel');
        
        
        [xl yl tl] = size(LCL_temp_ml);
        nccreate(ncfile,'LCL_temp_ml',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LCL_temp_ml',LCL_temp_ml) ;
        ncwriteatt(ncfile,'LCL_temp_ml','units','K');
        ncwriteatt(ncfile,'LCL_temp_ml','description','Environmental temperature at the condensation level for lifted ML parcel');
        
        
        [xl yl tl] = size(EL_height_ml);
        nccreate(ncfile,'EL_height_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'EL_height_ml',EL_height_ml) ;
        ncwriteatt(ncfile,'EL_height_ml','units','m');
        ncwriteatt(ncfile,'EL_height_ml','description','Height above ground of equilibrium level for lifted ML parcel (reversible warm moist lift)');
        
        [xl yl tl] = size(pEL_height_ml);
        nccreate(ncfile,'pEL_height_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pEL_height_ml',pEL_height_ml) ;
        ncwriteatt(ncfile,'pEL_height_ml','units','m');
        ncwriteatt(ncfile,'pEL_height_ml','description','Height above ground of equilibrium level for lifted ML parcel (pseudoadiabatic warm moist lift)');


%         [xl yl tl] = size(EL_pres_ml);
%         nccreate(ncfile,'EL_pres_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EL_pres_ml',EL_pres_ml) ;
%         ncwriteatt(ncfile,'EL_pres_ml','units','hPa');
%         ncwriteatt(ncfile,'EL_pres_ml','description','Pressure level of equilibrium level for lifted ML parcel');
        
        [xl yl tl] = size(EL_temp_ml);
        nccreate(ncfile,'EL_temp_ml',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'EL_temp_ml',EL_temp_ml) ;
        ncwriteatt(ncfile,'EL_temp_ml','units','K');
        ncwriteatt(ncfile,'EL_temp_ml','description','Environmental temperature at the equilibrium level for lifted ML parcel');
        
        [xl yl tl] = size(pEL_temp_ml);
        nccreate(ncfile,'pEL_temp_ml',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pEL_temp_ml',pEL_temp_ml) ;
        ncwriteatt(ncfile,'pEL_temp_ml','units','K');
        ncwriteatt(ncfile,'pEL_temp_ml','description','Environmental temperature at the equilibrium level for pseudoadiabatic lifted ML parcel');


        [xl yl tl] = size(CAPEacbl_ml);
        nccreate(ncfile,'CAPEacbl_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CAPEacbl_ml',CAPEacbl_ml) ;
        ncwriteatt(ncfile,'CAPEacbl_ml','units','J/kg');
        ncwriteatt(ncfile,'CAPEacbl_ml','description','Integration of parcel buoyancy in the layer between the moist LFC (of lifted ML parcel) and 1.5 km above it (reversible warm moist lift)');
       
        [xl yl tl] = size(pCAPEacbl_ml);
        nccreate(ncfile,'pCAPEacbl_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCAPEacbl_ml',pCAPEacbl_ml) ;
        ncwriteatt(ncfile,'pCAPEacbl_ml','units','J/kg');
        ncwriteatt(ncfile,'pCAPEacbl_ml','description','Integration of parcel buoyancy in the layer between the moist LFC (of lifted ML parcel) and 1.5 km above it (pseudoadiabatic warm moist lift)');
	
%         [xl yl tl] = size(CAPElcl_IB_ml);
%         nccreate(ncfile,'CAPElcl_IB_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CAPElcl_IB_ml',CAPElcl_IB_ml) ;
%         ncwriteatt(ncfile,'CAPElcl_IB_ml','units','J/kg');
%         ncwriteatt(ncfile,'CAPElcl_IB_ml','description','Integration of total parcel buoyancy (positive + negative) between the LCL (for a lifted ML parcel) and 2 km above it');
        
        [xl yl tl] = size(pCAPE_ml);
        nccreate(ncfile,'pCAPE_ml',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCAPE_ml',pCAPE_ml) ;
        ncwriteatt(ncfile,'pCAPE_ml','units','J/kg');
        ncwriteatt(ncfile,'pCAPE_ml','description','Integration of buoyancy between the moist LFC and EL for lifted ml parcel (warm moist pseudoadiabatic lift, no entrainemnt)');
        
        [xl yl tl] = size(pCIN_ml);
        nccreate(ncfile,'pCIN_ml',...		
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
	     'FillValue','disable');
        ncwrite(ncfile,'pCIN_ml',pCIN_ml) ;
        ncwriteatt(ncfile,'pCIN_ml','units','J/kg');
        ncwriteatt(ncfile,'pCIN_ml','description','Integration of negative buoyancy between the sfc and LFC for a lifted ml parcel (warm moist pseudoadiabatic lift, no entrainemnt)');






        
        
        %%%%%%%%%%% MU parcel
        
        [xl yl tl] = size(tallenough_mu);
        nccreate(ncfile,'tallenough_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'tallenough_mu',tallenough_mu) ;
        ncwriteatt(ncfile,'tallenough_mu','units','unitless');
        ncwriteatt(ncfile,'tallenough_mu','description','0 = temperature profile not deep enough at this time to capture EL for lifted MU parcel; 1 = profile is deep enough');
        
        
         [xl yl tl] = size(ptallenough_mu);
         nccreate(ncfile,'ptallenough_mu',...
              'Dimensions', {'location',xl,'time',yl,'cell',1},...
             'FillValue','disable');
         ncwrite(ncfile,'ptallenough_mu',ptallenough_mu) ;
         ncwriteatt(ncfile,'ptallenough_mu','units','unitless');
         ncwriteatt(ncfile,'ptallenough_mu','description','0 = temperature profile not deep enough at this time to capture EL for lifted mu pseudoadiabatic parcel; 1 = profile is deep enough');

        [xl yl tl] = size(CAPE_mu);
        nccreate(ncfile,'CAPE_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CAPE_mu',CAPE_mu) ;
        ncwriteatt(ncfile,'CAPE_mu','units','J/kg');
        ncwriteatt(ncfile,'CAPE_mu','description','Integration of buoyancy between the moist LFC and EL for lifted MU parcel (reversible warm moist lift)');
        
        [xl yl tl] = size(CIN_mu);
        nccreate(ncfile,'CIN_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CIN_mu',CIN_mu) ;
        ncwriteatt(ncfile,'CIN_mu','units','J/kg');
        ncwriteatt(ncfile,'CIN_mu','description','Integration of negative buoyancy between the origin height and moist LFC for a lifted MU parcel (reversible warm moist lift)');
        
%         [xl yl tl] = size(CIN_NA_fract_mu);
%         nccreate(ncfile,'CIN_NA_fract_mu',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CIN_NA_fract_mu',CIN_NA_fract_mu) ;
%         ncwriteatt(ncfile,'CIN_NA_fract_mu','units','unitless');
%         ncwriteatt(ncfile,'CIN_NA_fract_mu','description','Fraction of the number of data points between the origin height and LFC (for a lifted MU parcel) containing negative buoyancy (reversible warm moist lift)');
        
%         [xl yl tl] = size(CIN_IB_mu);
%         nccreate(ncfile,'CIN_IB_mu',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CIN_IB_mu',CIN_IB_mu) ;
%         ncwriteatt(ncfile,'CIN_IB_mu','units','J/kg');
%         ncwriteatt(ncfile,'CIN_IB_mu','description','Integration of total parcel buoyancy (positive + negative) between the origin height and LFC for a lifted MU parcel');
        
        [xl yl tl] = size(LFC_height_mu);
        nccreate(ncfile,'LFC_height_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LFC_height_mu',LFC_height_mu) ;
        ncwriteatt(ncfile,'LFC_height_mu','units','m');
        ncwriteatt(ncfile,'LFC_height_mu','description','Height above ground of the moist level of free convection for lifted MU parcel (reversible warm moist lift)');
       
        [xl yl tl] = size(pLFC_height_mu);
        nccreate(ncfile,'pLFC_height_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pLFC_height_mu',pLFC_height_mu) ;
        ncwriteatt(ncfile,'pLFC_height_mu','units','m');
        ncwriteatt(ncfile,'pLFC_height_mu','description','Height above ground of the moist level of free convection for lifted MU parcel (pseudoadiabatic warm moist lift)');

%         [xl yl tl] = size(LFC_pres_mu);
%         nccreate(ncfile,'LFC_pres_mu',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'LFC_pres_mu',LFC_pres_mu) ;
%         ncwriteatt(ncfile,'LFC_pres_mu','units','hPa');
%         ncwriteatt(ncfile,'LFC_pres_mu','description','Pressure level of the level of free convection for lifted MU parcel');
        
        [xl yl tl] = size(LFC_temp_mu);
        nccreate(ncfile,'LFC_temp_mu',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LFC_temp_mu',LFC_temp_mu) ;
        ncwriteatt(ncfile,'LFC_temp_mu','units','K');
        ncwriteatt(ncfile,'LFC_temp_mu','description','Environmental temperature at the level of free convection for lifted MU parcel');
        
        [xl yl tl] = size(pLFC_temp_mu);
        nccreate(ncfile,'pLFC_temp_mu',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pLFC_temp_mu',pLFC_temp_mu) ;
        ncwriteatt(ncfile,'pLFC_temp_mu','units','K');
        ncwriteatt(ncfile,'pLFC_temp_mu','description','Environmental temperature at the level of free convection for pseudoadiabatic lifted MU parcel');


        % [xl yl tl] = size(dZstar_mu);
        % nccreate(ncfile,'dZstar_mu',...
        %         'Dimensions', {'x',xl,'y',yl,'t',tl,'cell',1},...
        %     'FillValue','disable');
        % ncwrite(ncfile,'dZstar_mu',dZstar_mu) ;
        % ncwriteatt(ncfile,'dZstar_mu','units','m');
        % ncwriteatt(ncfile,'dZstar_mu','description','Difference between the LFC and the origin height of a MU parcel');
        
        
        [xl yl tl] = size(LCL_height_mu);
        nccreate(ncfile,'LCL_height_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LCL_height_mu',LCL_height_mu) ;
        ncwriteatt(ncfile,'LCL_height_mu','units','m');
        ncwriteatt(ncfile,'LCL_height_mu','description','Height above ground of condensation level for lifted MU parcel');
        
%         [xl yl tl] = size(LCL_pres_mu);
%         nccreate(ncfile,'LCL_pres_mu',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'LCL_pres_mu',LCL_pres_mu) ;
%         ncwriteatt(ncfile,'LCL_pres_mu','units','hPa');
%         ncwriteatt(ncfile,'LCL_pres_mu','description','Pressure level of condensation level for lifted MU parcel');
        
        [xl yl tl] = size(LCL_temp_mu);
        nccreate(ncfile,'LCL_temp_mu',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'LCL_temp_mu',LCL_temp_mu) ;
        ncwriteatt(ncfile,'LCL_temp_mu','units','K');
        ncwriteatt(ncfile,'LCL_temp_mu','description','Environmental temperature at the condensation level for lifted MU parcel');
        
        [xl yl tl] = size(EL_height_mu);
        nccreate(ncfile,'EL_height_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'EL_height_mu',EL_height_mu) ;
        ncwriteatt(ncfile,'EL_height_mu','units','m');
        ncwriteatt(ncfile,'EL_height_mu','description','Height above ground of equilibrium level for lifted MU parcel (reversible warm moist lift)');
       
        [xl yl tl] = size(pEL_height_mu);
        nccreate(ncfile,'pEL_height_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pEL_height_mu',pEL_height_mu) ;
        ncwriteatt(ncfile,'pEL_height_mu','units','m');
        ncwriteatt(ncfile,'pEL_height_mu','description','Height above ground of equilibrium level for lifted MU parcel (pseudoadiabatic warm moist lift)');

%         [xl yl tl] = size(EL_pres_mu);
%         nccreate(ncfile,'EL_pres_mu',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EL_pres_mu',EL_pres_mu) ;
%         ncwriteatt(ncfile,'EL_pres_mu','units','hPa');
%         ncwriteatt(ncfile,'EL_pres_mu','description','Pressure level of equilibrium level for lifted MU parcel');
        
        [xl yl tl] = size(EL_temp_mu);
        nccreate(ncfile,'EL_temp_mu',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'EL_temp_mu',EL_temp_mu) ;
        ncwriteatt(ncfile,'EL_temp_mu','units','K');
        ncwriteatt(ncfile,'EL_temp_mu','description','Environmental temperature at the equilibrium level for lifted MU parcel');
        
        [xl yl tl] = size(pEL_temp_mu);
        nccreate(ncfile,'pEL_temp_mu',...
                 'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pEL_temp_mu',pEL_temp_mu) ;
        ncwriteatt(ncfile,'pEL_temp_mu','units','K');
        ncwriteatt(ncfile,'pEL_temp_mu','description','Environmental temperature at the equilibrium level for pseudoadiabatic lifted MU parcel');


        [xl yl tl] = size(CAPEacbl_mu);
        nccreate(ncfile,'CAPEacbl_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'CAPEacbl_mu',CAPEacbl_mu) ;
        ncwriteatt(ncfile,'CAPEacbl_mu','units','J/kg');
        ncwriteatt(ncfile,'CAPEacbl_mu','description','Integration of parcel buoyancy in the layer between the moist LFC (of lifted MU parcel) and 1.5 km above it (reversible warm moist lift)');
        
        [xl yl tl] = size(pCAPEacbl_mu);
        nccreate(ncfile,'pCAPEacbl_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCAPEacbl_mu',pCAPEacbl_mu) ;
        ncwriteatt(ncfile,'pCAPEacbl_mu','units','J/kg');
        ncwriteatt(ncfile,'pCAPEacbl_mu','description','Integration of parcel buoyancy in the layer between the moist LFC (of lifted MU parcel) and 1.5 km above it (pseudoadiabatic warm moist lift)');

%         [xl yl tl] = size(CAPElcl_IB_mu);
%         nccreate(ncfile,'CAPElcl_IB_mu',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'CAPElcl_IB_mu',CAPElcl_IB_mu) ;
%         ncwriteatt(ncfile,'CAPElcl_IB_mu','units','J/kg');
%         ncwriteatt(ncfile,'CAPElcl_IB_mu','description','Integration of total parcel buoyancy (positive + negative) between the LCL (for a lifted MU parcel) and 2 km above it');
        

        [xl yl tl] = size(pCAPE_mu);
        nccreate(ncfile,'pCAPE_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCAPE_mu',pCAPE_mu) ;
        ncwriteatt(ncfile,'pCAPE_mu','units','J/kg');
        ncwriteatt(ncfile,'pCAPE_mu','description','Integration of buoyancy between the moist LFC and EL for lifted mu parcel (warm moist pseudoadiabatic lift, no entrainemnt)');
        
        [xl yl tl] = size(pCIN_mu);
        nccreate(ncfile,'pCIN_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'pCIN_mu',pCIN_mu) ;
        ncwriteatt(ncfile,'pCIN_mu','units','J/kg');
        ncwriteatt(ncfile,'pCIN_mu','description','Integration of negative buoyancy between the sfc and LFC for a lifted mu parcel (warm moist pseudoadiabatic lift, no entrainemnt)');




        [xl yl tl] = size(initial_ht_parcel_mu);
        nccreate(ncfile,'initial_ht_parcel_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'initial_ht_parcel_mu',initial_ht_parcel_mu) ;
        ncwriteatt(ncfile,'initial_ht_parcel_mu','units','m');
        ncwriteatt(ncfile,'initial_ht_parcel_mu','description','Height above ground level of the MU parcel origin');
        



        [xl yl tl] = size(theta_mu);
        nccreate(ncfile,'theta_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'theta_mu',theta_mu) ;
        ncwriteatt(ncfile,'theta_mu','units','K');
        ncwriteatt(ncfile,'theta_mu','description','Potential temperature at the origin height of a MU parcel');

	[xl yl tl] = size(thetav_mu);
        nccreate(ncfile,'thetav_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'thetav_mu',thetav_mu) ;
        ncwriteatt(ncfile,'thetav_mu','units','K');
        ncwriteatt(ncfile,'thetav_mu','description','Virtual potential temperature at the origin height of a MU parcel');



	[xl yl tl] = size(thetae_mu);
        nccreate(ncfile,'thetae_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'thetae_mu',thetae_mu) ;
        ncwriteatt(ncfile,'thetae_mu','units','K');
        ncwriteatt(ncfile,'thetae_mu','description','Equivalent potential temperature at the origin height of a MU parcel');
        
        [xl yl tl] = size(thetae_mean_subcloud_mu);
        nccreate(ncfile,'thetae_mean_subcloud_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'thetae_mean_subcloud_mu',thetae_mean_subcloud_mu) ;
        ncwriteatt(ncfile,'thetae_mean_subcloud_mu','units','K');
        ncwriteatt(ncfile,'thetae_mean_subcloud_mu','description','Mean equivalent potential temperature below the LCL for a MU parcel');
        
        [xl yl tl] = size(thetae_mean_ACBL_mu);
        nccreate(ncfile,'thetae_mean_ACBL_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'thetae_mean_ACBL_mu',thetae_mean_ACBL_mu) ;
        ncwriteatt(ncfile,'thetae_mean_ACBL_mu','units','K');
        ncwriteatt(ncfile,'thetae_mean_ACBL_mu','description','Mean equivalent potential temperature within the active cloud bearing layer for a MU parcel');
       

        [xl yl tl] = size(U_mu);
        nccreate(ncfile,'U_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'U_mu',U_mu) ;
        ncwriteatt(ncfile,'U_mu','units','m/s');
        ncwriteatt(ncfile,'U_mu','description',' Ground-relative zonal wind at the origin height of the MU parcel');

        [xl yl tl] = size(V_mu);
        nccreate(ncfile,'V_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'V_mu',V_mu) ;
        ncwriteatt(ncfile,'V_mu','units','m/s');
        ncwriteatt(ncfile,'V_mu','description',' Ground-relative meridional wind at the origin height of the MU parcel');	
        
        [xl yl tl] = size(Ucrel_mu);
        nccreate(ncfile,'Ucrel_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'Ucrel_mu',Ucrel_mu) ;
        ncwriteatt(ncfile,'Ucrel_mu','units','m/s');
        ncwriteatt(ncfile,'Ucrel_mu','description',' Cloud-relative zonal wind at the origin height of the MU parcel');
        
        [xl yl tl] = size(Vcrel_mu);
        nccreate(ncfile,'Vcrel_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'Vcrel_mu',Vcrel_mu) ;
        ncwriteatt(ncfile,'Vcrel_mu','units','m/s');
        ncwriteatt(ncfile,'Vcrel_mu','description',' Cloud-relative meridional wind at the origin height of the MU parcel');
        
        [xl yl tl] = size(RH_mean_ACBL);
        nccreate(ncfile,'RH_mean_ACBL',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'RH_mean_ACBL',RH_mean_ACBL) ;
        ncwriteatt(ncfile,'RH_mean_ACBL','units','0-100 0-100 percent');
        ncwriteatt(ncfile,'RH_mean_ACBL','description','Mean relative humidity within the active cloud bearing layer for a MU parcel.');
        

        [xl yl tl] = size(rh_at_muLCLplus1500m);
        nccreate(ncfile,'rh_at_muLCLplus1500m',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_at_muLCLplus1500m',rh_at_muLCLplus1500m) ;
        ncwriteatt(ncfile,'rh_at_muLCLplus1500m','units','0-100 0-100 percent');
        ncwriteatt(ncfile,'rh_at_muLCLplus1500m','description','Mean relative humidity at z = 1.5km above the muLCL.');


        [xl yl tl] = size(rh_at_muLFCplus1500m);
        nccreate(ncfile,'rh_at_muLFCplus1500m',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_at_muLFCplus1500m',rh_at_muLFCplus1500m) ;
        ncwriteatt(ncfile,'rh_at_muLFCplus1500m','units','0-100 0-100 percent');
        ncwriteatt(ncfile,'rh_at_muLFCplus1500m','description','Mean relative humidity at z = 1.5km above the muLFC.');

        [xl yl tl] = size(rh_at_pmuLFCplus1500m);
        nccreate(ncfile,'rh_at_pmuLFCplus1500m',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_at_pmuLFCplus1500m',rh_at_pmuLFCplus1500m) ;
        ncwriteatt(ncfile,'rh_at_pmuLFCplus1500m','units','0-100 0-100 percent');
        ncwriteatt(ncfile,'rh_at_pmuLFCplus1500m','description','Mean relative humidity at z = 1.5km above the pseudoadiabatic muLFC.');



        [xl yl tl] = size(rvap_mu);
        nccreate(ncfile,'rvap_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_mu',rvap_mu) ;
        ncwriteatt(ncfile,'rvap_mu','units','g/kg');
        ncwriteatt(ncfile,'rvap_mu','description','Vaper mixing ratio at the origin height of a MU parcel');
        
        [xl yl tl] = size(shear_mag_bulk_FT_mu);
        nccreate(ncfile,'shear_mag_bulk_FT_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_FT_mu',shear_mag_bulk_FT_mu) ;
        ncwriteatt(ncfile,'shear_mag_bulk_FT_mu','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_FT_mu','description','Magnitude of the difference of horizontal wind vectors spanning the LFC and EL for a mu parcel');
        
        [xl yl tl] = size(shear_dir_bulk_FT_mu);
        nccreate(ncfile,'shear_dir_bulk_FT_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_dir_bulk_FT_mu',shear_dir_bulk_FT_mu) ;
        ncwriteatt(ncfile,'shear_dir_bulk_FT_mu','units','degrees');
        ncwriteatt(ncfile,'shear_dir_bulk_FT_mu','description','Direction [met convention - dir from, where 0deg = N] of the difference of horizontal wind vectors spanning the LFC and EL for MU parcel');
        
        [xl yl tl] = size(shear_mag_bulk_ACBL_mu);
        nccreate(ncfile,'shear_mag_bulk_ACBL_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_ACBL_mu',shear_mag_bulk_ACBL_mu) ;
        ncwriteatt(ncfile,'shear_mag_bulk_ACBL_mu','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_ACBL_mu','description','Magnitude of the difference of horizontal wind vectors spanning the top and bottom of the active cloud bearing layer for a MU parcel');
        
        
        
        [xl yl tl] = size(shear_mag_bulk_ACBL_pmu);
        nccreate(ncfile,'shear_mag_bulk_ACBL_pmu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_ACBL_pmu',shear_mag_bulk_ACBL_pmu) ;
        ncwriteatt(ncfile,'shear_mag_bulk_ACBL_pmu','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_ACBL_pmu','description','Magnitude of the difference of horizontal wind vectors spanning the top and bottom of the active cloud bearing layer for a psuedoadiabatic MU parcel');
        

        [xl yl tl] = size(RH_mean_ACBL_pmu);
        nccreate(ncfile,'RH_mean_ACBL_pmu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'RH_mean_ACBL_pmu',RH_mean_ACBL_pmu) ;
        ncwriteatt(ncfile,'RH_mean_ACBL_pmu','units','0-100 percent');
        ncwriteatt(ncfile,'RH_mean_ACBL_pmu','description','Mean relative humidity within the active cloud bearing layer for a pseudoadiabatic MU parcel.');


        [xl yl tl] = size(thetae_mean_ACBL_pmu);
        nccreate(ncfile,'thetae_mean_ACBL_pmu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'thetae_mean_ACBL_pmu',thetae_mean_ACBL_pmu) ;
        ncwriteatt(ncfile,'thetae_mean_ACBL_pmu','units','K');
        ncwriteatt(ncfile,'thetae_mean_ACBL_pmu','description','Mean equivalent potential temperature within the active cloud bearing layer for a pseudoadiabatic MU parcel.');

        [xl yl tl] = size(Ucrel_mean_ACBL_pmu);
        nccreate(ncfile,'Ucrel_mean_ACBL_pmu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'Ucrel_mean_ACBL_pmu',Ucrel_mean_ACBL_pmu) ;
        ncwriteatt(ncfile,'Ucrel_mean_ACBL_pmu','units','m/s');
        ncwriteatt(ncfile,'Ucrel_mean_ACBL_pmu','description','Mean cloud-rel zonal wind within the active cloud bearing layer for a pseudoadiabatic MU parcel.');
      
        [xl yl tl] = size(Vcrel_mean_ACBL_pmu);
        nccreate(ncfile,'Vcrel_mean_ACBL_pmu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'Vcrel_mean_ACBL_pmu',Vcrel_mean_ACBL_pmu) ;
        ncwriteatt(ncfile,'Vcrel_mean_ACBL_pmu','units','m/s');
        ncwriteatt(ncfile,'Vcrel_mean_ACBL_pmu','description','Mean cloud-rel meridional wind within the active cloud bearing layer for a pseudoadiabatic MU parcel.');        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %      thermodynamics/moisture:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [xl yl tl] = size(freezing_ht);
        nccreate(ncfile,'freezing_ht',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'freezing_ht',freezing_ht) ;
        ncwriteatt(ncfile,'freezing_ht','units','m');
        ncwriteatt(ncfile,'freezing_ht','description','Height above ground of T = 0 deg C');
        
        [xl yl tl] = size(PW);
        nccreate(ncfile,'PW',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'PW',PW) ;
        ncwriteatt(ncfile,'PW','units','cm');
        ncwriteatt(ncfile,'PW','description','Precipitable water calculated from profile vapor mixing ratio profile');
        
        
        
        [xl yl tl] = size(rvap_925mb);
        nccreate(ncfile,'rvap_925mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_925mb',rvap_925mb) ;
        ncwriteatt(ncfile,'rvap_925mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_925mb','description','Vapor mixing ratio at 925 hPa pressure level');        
        
        [xl yl tl] = size(rh_925mb);
        nccreate(ncfile,'rh_925mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_925mb',rh_925mb) ;
        ncwriteatt(ncfile,'rh_925mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_925mb','description','Relative humidity at 925 hPa pressure level');
        
        [xl yl tl] = size(rvap_850mb);
        nccreate(ncfile,'rvap_850mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_850mb',rvap_850mb) ;
        ncwriteatt(ncfile,'rvap_850mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_850mb','description','Vapor mixing ratio at 850 hPa pressure level');
        
        [xl yl tl] = size(rh_850mb);
        nccreate(ncfile,'rh_850mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_850mb',rh_850mb) ;
        ncwriteatt(ncfile,'rh_850mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_850mb','description','Relative humidity at 850 hPa pressure level');
        
        [xl yl tl] = size(rvap_700mb);
        nccreate(ncfile,'rvap_700mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_700mb',rvap_700mb) ;
        ncwriteatt(ncfile,'rvap_700mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_700mb','description','Vapor mixing ratio at 700 hPa pressure level');
        
        [xl yl tl] = size(rh_700mb);
        nccreate(ncfile,'rh_700mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_700mb',rh_700mb) ;
        ncwriteatt(ncfile,'rh_700mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_700mb','description','Relative humidity at 700 hPa pressure level');
        
        
        [xl yl tl] = size(rvap_600mb);
        nccreate(ncfile,'rvap_600mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_600mb',rvap_600mb) ;
        ncwriteatt(ncfile,'rvap_600mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_600mb','description','Vapor mixing ratio at 600 hPa pressure level');
        
        [xl yl tl] = size(rh_600mb);
        nccreate(ncfile,'rh_600mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_600mb',rh_600mb) ;
        ncwriteatt(ncfile,'rh_600mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_600mb','description','Relative humidity at 600 hPa pressure level');
        
        
        [xl yl tl] = size(rvap_500mb);
        nccreate(ncfile,'rvap_500mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_500mb',rvap_500mb) ;
        ncwriteatt(ncfile,'rvap_500mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_500mb','description','Vapor mixing ratio at 500 hPa pressure level');
        
        [xl yl tl] = size(rh_500mb);
        nccreate(ncfile,'rh_500mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_500mb',rh_500mb) ;
        ncwriteatt(ncfile,'rh_500mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_500mb','description','Relative humidity at 500 hPa pressure level');
        
        
        [xl yl tl] = size(rvap_400mb);
        nccreate(ncfile,'rvap_400mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_400mb',rvap_400mb) ;
        ncwriteatt(ncfile,'rvap_400mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_400mb','description','Vapor mixing ratio at 400 hPa pressure level');
        
        [xl yl tl] = size(rh_400mb);
        nccreate(ncfile,'rh_400mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_400mb',rh_400mb) ;
        ncwriteatt(ncfile,'rh_400mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_400mb','description','Relative humidity at 400 hPa pressure level');


        [xl yl tl] = size(rvap_300mb);
        nccreate(ncfile,'rvap_300mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_300mb',rvap_300mb) ;
        ncwriteatt(ncfile,'rvap_300mb','units','g/kg');
        ncwriteatt(ncfile,'rvap_300mb','description','Vapor mixing ratio at 300 hPa pressure level');
        
        [xl yl tl] = size(rh_300mb);
        nccreate(ncfile,'rh_300mb',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_300mb',rh_300mb) ;
        ncwriteatt(ncfile,'rh_300mb','units','0-100 percent');
        ncwriteatt(ncfile,'rh_300mb','description','Relative humidity at 300 hPa pressure level');



        [xl yl tl] = size(rvap_3000masl);
        nccreate(ncfile,'rvap_3000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_3000masl',rvap_3000masl) ;
        ncwriteatt(ncfile,'rvap_3000masl','units','g/kg');
        ncwriteatt(ncfile,'rvap_3000masl','description','Vapor mixing ratio at 3km asl');
        
        [xl yl tl] = size(rh_3000masl);
        nccreate(ncfile,'rh_3000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_3000masl',rh_3000masl) ;
        ncwriteatt(ncfile,'rh_3000masl','units','0-100 percent');
        ncwriteatt(ncfile,'rh_3000masl','description','Relative humidity at 3km asl');


        [xl yl tl] = size(rvap_4000masl);
        nccreate(ncfile,'rvap_4000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_4000masl',rvap_4000masl) ;
        ncwriteatt(ncfile,'rvap_4000masl','units','g/kg');
        ncwriteatt(ncfile,'rvap_4000masl','description','Vapor mixing ratio at 4km asl');
        

        [xl yl tl] = size(rh_4000masl);
        nccreate(ncfile,'rh_4000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_4000masl',rh_4000masl) ;
        ncwriteatt(ncfile,'rh_4000masl','units','0-100 percent');
        ncwriteatt(ncfile,'rh_4000masl','description','Relative humidity at 4km asl');


        [xl yl tl] = size(rvap_5000masl);
        nccreate(ncfile,'rvap_5000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_5000masl',rvap_5000masl) ;
        ncwriteatt(ncfile,'rvap_5000masl','units','g/kg');
        ncwriteatt(ncfile,'rvap_5000masl','description','Vapor mixing ratio at 5km asl');
        
        
        [xl yl tl] = size(rh_5000masl);
        nccreate(ncfile,'rh_5000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_5000masl',rh_5000masl) ;
        ncwriteatt(ncfile,'rh_5000masl','units','0-100 percent');
        ncwriteatt(ncfile,'rh_5000masl','description','Relative humidity at 5km asl');


        [xl yl tl] = size(rvap_6000masl);
        nccreate(ncfile,'rvap_6000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rvap_6000masl',rvap_6000masl) ;
        ncwriteatt(ncfile,'rvap_6000masl','units','g/kg');
        ncwriteatt(ncfile,'rvap_6000masl','description','Vapor mixing ratio at 6km asl');
        
        
        [xl yl tl] = size(rh_6000masl);
        nccreate(ncfile,'rh_6000masl',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'rh_6000masl',rh_6000masl) ;
        ncwriteatt(ncfile,'rh_6000masl','units','0-100 percent');
        ncwriteatt(ncfile,'rh_6000masl','description','Relative humidity at 6km asl');



%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%  Effective Inflow Layer data
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         
%         [xl yl tl] = size(EIL10_bot_height );
%         nccreate(ncfile,'EIL10_bot_height',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EIL10_bot_height',EIL10_bot_height ) ;
%         ncwriteatt(ncfile,'EIL10_bot_height','units','m');
%         ncwriteatt(ncfile,'EIL10_bot_height','description','Height above ground of the bottom of the layer composed of parcels with CAPE  > 100 and CIN_IB < 10 J/kg');
%         
%         [xl yl tl] = size(EIL10_top_height );
%         nccreate(ncfile,'EIL10_top_height',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EIL10_top_height',EIL10_top_height ) ;
%         ncwriteatt(ncfile,'EIL10_top_height','units','m');
%         ncwriteatt(ncfile,'EIL10_top_height','description','Height above ground of the top of the layer composed of parcels with CAPE  > 100 and CIN_IB < 10 J/kg');
%         
%         [xl yl tl] = size(EIL25_bot_height );
%         nccreate(ncfile,'EIL25_bot_height',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EIL25_bot_height',EIL25_bot_height ) ;
%         ncwriteatt(ncfile,'EIL25_bot_height','units','m');
%         ncwriteatt(ncfile,'EIL25_bot_height','description','Height above ground of the bottom of the layer composed of parcels with CAPE  > 100 and CIN_IB < 25 J/kg');
%         
%         [xl yl tl] = size(EIL25_top_height );
%         nccreate(ncfile,'EIL25_top_height',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EIL25_top_height',EIL25_top_height ) ;
%         ncwriteatt(ncfile,'EIL25_top_height','units','m');
%         ncwriteatt(ncfile,'EIL25_top_height','description','Height above ground of the top of the layer composed of parcels with CAPE  > 100 and CIN_IB < 25 J/kg');
%         
%         [xl yl tl] = size(EIL50_bot_height );
%         nccreate(ncfile,'EIL50_bot_height',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EIL50_bot_height',EIL50_bot_height ) ;
%         ncwriteatt(ncfile,'EIL50_bot_height','units','m');
%         ncwriteatt(ncfile,'EIL50_bot_height','description','Height above ground of the bottom of the layer composed of parcels with CAPE  > 100 and CIN_IB < 50 J/kg');
%         
%         [xl yl tl] = size(EIL50_top_height );
%         nccreate(ncfile,'EIL50_top_height',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'EIL50_top_height',EIL50_top_height ) ;
%         ncwriteatt(ncfile,'EIL50_top_height','units','m');
%         ncwriteatt(ncfile,'EIL50_top_height','description','Height above ground of the top of the layer composed of parcels with CAPE  > 100 and CIN_IB < 50 J/kg');
%         
%         [xl yl tl] = size(rvap_mean_EIL10);
%         nccreate(ncfile,'rvap_mean_EIL10',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'rvap_mean_EIL10',rvap_mean_EIL10) ;
%         ncwriteatt(ncfile,'rvap_mean_EIL10','units','g/kg');
%         ncwriteatt(ncfile,'rvap_mean_EIL10','description','Mean vapor mixing ratio within the EIL with CIN < -10 J/kg');
%         
%         [xl yl tl] = size(rvap_mean_EIL25);
%         nccreate(ncfile,'rvap_mean_EIL25',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'rvap_mean_EIL25',rvap_mean_EIL25) ;
%         ncwriteatt(ncfile,'rvap_mean_EIL25','units','g/kg');
%         ncwriteatt(ncfile,'rvap_mean_EIL25','description','Mean vapor mixing ratio within the EIL with CIN < -25 J/kg');
%         
%         [xl yl tl] = size(rvap_mean_EIL50);
%         nccreate(ncfile,'rvap_mean_EIL50',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'rvap_mean_EIL50',rvap_mean_EIL50) ;
%         ncwriteatt(ncfile,'rvap_mean_EIL50','units','g/kg');
%         ncwriteatt(ncfile,'rvap_mean_EIL50','description','Mean vapor mixing ratio within the EIL with CIN < -50 J/kg');
  
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        % etc wind and shear
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [xl yl tl] = size(shear_mag_bulk_0to1km);
        nccreate(ncfile,'shear_mag_bulk_0to1km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_0to1km',shear_mag_bulk_0to1km) ;
        ncwriteatt(ncfile,'shear_mag_bulk_0to1km','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_0to1km','description','Magnitude of the difference of horizontal wind vectors at 1 km and 10 m');
        
        [xl yl tl] = size(shear_dir_bulk_0to1km);
        nccreate(ncfile,'shear_dir_bulk_0to1km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_dir_bulk_0to1km',shear_dir_bulk_0to1km) ;
        ncwriteatt(ncfile,'shear_dir_bulk_0to1km','units','degrees');
        ncwriteatt(ncfile,'shear_dir_bulk_0to1km','description','Direction [met convention - dir from, where 0deg = N] of the difference of horizontal wind vectors at 1 km and 10 m');
        
        
        [xl yl tl] = size(shear_mag_bulk_0to3km);
        nccreate(ncfile,'shear_mag_bulk_0to3km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_0to3km',shear_mag_bulk_0to3km) ;
        ncwriteatt(ncfile,'shear_mag_bulk_0to3km','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_0to3km','description','Magnitude of the difference of horizontal wind vectors at 3 km and 10 m');
        
        [xl yl tl] = size(shear_dir_bulk_0to3km);
        nccreate(ncfile,'shear_dir_bulk_0to3km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_dir_bulk_0to3km',shear_dir_bulk_0to3km) ;
        ncwriteatt(ncfile,'shear_dir_bulk_0to3km','units','degrees');
        ncwriteatt(ncfile,'shear_dir_bulk_0to3km','description','Direction [met convention - dir from, where 0deg = N] of the difference of horizontal wind vectors at 3 km and 10 m');
        
        [xl yl tl] = size(shear_mag_bulk_0to6km);
        nccreate(ncfile,'shear_mag_bulk_0to6km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_0to6km',shear_mag_bulk_0to6km) ;
        ncwriteatt(ncfile,'shear_mag_bulk_0to6km','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_0to6km','description','Magnitude of the difference of horizontal wind vectors at 6 km and 10 m');
        
        [xl yl tl] = size(shear_dir_bulk_0to6km);
        nccreate(ncfile,'shear_dir_bulk_0to6km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_dir_bulk_0to6km',shear_dir_bulk_0to6km) ;
        ncwriteatt(ncfile,'shear_dir_bulk_0to6km','units','degrees');
        ncwriteatt(ncfile,'shear_dir_bulk_0to6km','description','Direction [met convention - dir from, where 0deg = N] of the difference of horizontal wind vectors at 6 km and 10 m');
        

        [xl yl tl] = size(shear_mag_bulk_0to9km);
        nccreate(ncfile,'shear_mag_bulk_0to9km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_mag_bulk_0to9km',shear_mag_bulk_0to9km) ;
        ncwriteatt(ncfile,'shear_mag_bulk_0to9km','units','m/s');
        ncwriteatt(ncfile,'shear_mag_bulk_0to9km','description','Magnitude of the difference of horizontal wind vectors at 9 km and 10 m');
        
        [xl yl tl] = size(shear_dir_bulk_0to9km);
        nccreate(ncfile,'shear_dir_bulk_0to9km',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'shear_dir_bulk_0to9km',shear_dir_bulk_0to9km) ;
        ncwriteatt(ncfile,'shear_dir_bulk_0to9km','units','degrees');
        ncwriteatt(ncfile,'shear_dir_bulk_0to9km','description','Direction [met convention - dir from, where 0deg = N] of the difference of horizontal wind vectors at 9 km and 10 m');    
        
        
%         [al bl] = size(Ucrel_mean_EIL10);
%         nccreate(ncfile,'Ucrel_mean_EIL10',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_EIL10',Ucrel_mean_EIL10) ;
%         ncwriteatt(ncfile,'Ucrel_mean_EIL10','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_EIL10','description','Mean zonal cloud-relative wind within the EIL with CIN < 10 J/kg');
%         
%         [al bl] = size(Vcrel_mean_EIL10);
%         nccreate(ncfile,'Vcrel_mean_EIL10',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_EIL10',Vcrel_mean_EIL10) ;
%         ncwriteatt(ncfile,'Vcrel_mean_EIL10','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_EIL10','description','Mean meridional cloud-relative wind within the EIL with CIN < 10 J/kg');
%         
%         
%         [al bl] = size(Ucrel_mean_EIL25);
%         nccreate(ncfile,'Ucrel_mean_EIL25',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_EIL25',Ucrel_mean_EIL25) ;
%         ncwriteatt(ncfile,'Ucrel_mean_EIL25','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_EIL25','description','Mean zonal cloud-relative wind within the EIL with CIN < 25 J/kg');
%         
%         [al bl] = size(Vcrel_mean_EIL25);
%         nccreate(ncfile,'Vcrel_mean_EIL25',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_EIL25',Vcrel_mean_EIL25) ;
%         ncwriteatt(ncfile,'Vcrel_mean_EIL25','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_EIL25','description','Mean meridional cloud-relative wind within the EIL with CIN < 25 J/kg');
%         
%         
%         [al bl] = size(Ucrel_mean_EIL50);
%         nccreate(ncfile,'Ucrel_mean_EIL50',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_EIL50',Ucrel_mean_EIL50) ;
%         ncwriteatt(ncfile,'Ucrel_mean_EIL50','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_EIL50','description','Mean zonal cloud-relative wind within the EIL with CIN < 50 J/kg');
%         
%         [al bl] = size(Vcrel_mean_EIL50);
%         nccreate(ncfile,'Vcrel_mean_EIL50',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_EIL50',Vcrel_mean_EIL50) ;
%         ncwriteatt(ncfile,'Vcrel_mean_EIL50','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_EIL50','description','Mean meridional cloud-relative wind within the EIL with CIN < 50 J/kg');
        
        
%         [al bl] = size(Ucrel_mean_ACBL_sfc);
%         nccreate(ncfile,'Ucrel_mean_ACBL_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_ACBL_sfc',Ucrel_mean_ACBL_sfc) ;
%         ncwriteatt(ncfile,'Ucrel_mean_ACBL_sfc','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_ACBL_sfc','description',' Mean of the cloud-relative zonal wind in the active cloud bearing layer for a sfc parcel ');
%         
%         [al bl] = size(Vcrel_mean_ACBL_sfc);
%         nccreate(ncfile,'Vcrel_mean_ACBL_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_ACBL_sfc',Vcrel_mean_ACBL_sfc) ;
%         ncwriteatt(ncfile,'Vcrel_mean_ACBL_sfc','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_ACBL_sfc','description',' Mean of the cloud-relative meridional wind in the active cloud bearing layer for a sfc parcel ');
%         
%         [al bl] = size(Ucrel_mean_0toACBL_sfc);
%         nccreate(ncfile,'Ucrel_mean_0toACBL_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_0toACBL_sfc',Ucrel_mean_0toACBL_sfc) ;
%         ncwriteatt(ncfile,'Ucrel_mean_0toACBL_sfc','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_0toACBL_sfc','description',' Mean of the cloud-relative zonal wind between the surface and top of the active cloud bearing layer for a sfc parcel ');
%         
%         [al bl] = size(Vcrel_mean_0toACBL_sfc);
%         nccreate(ncfile,'Vcrel_mean_0toACBL_sfc',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_0toACBL_sfc',Vcrel_mean_0toACBL_sfc) ;
%         ncwriteatt(ncfile,'Vcrel_mean_0toACBL_sfc','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_0toACBL_sfc','description',' Mean of the cloud-relative meridional wind between the surface and top of the active cloud bearing layer for a sfc parcel ');
        
%         [al bl] = size(Ucrel_mean_ACBL_ml);
%         nccreate(ncfile,'Ucrel_mean_ACBL_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_ACBL_ml',Ucrel_mean_ACBL_ml) ;
%         ncwriteatt(ncfile,'Ucrel_mean_ACBL_ml','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_ACBL_ml','description',' Mean of the cloud-relative zonal wind in the active cloud bearing layer for a 100-hPa mean layer parcel ');
%         
%         [al bl] = size(Vcrel_mean_ACBL_ml);
%         nccreate(ncfile,'Vcrel_mean_ACBL_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_ACBL_ml',Vcrel_mean_ACBL_ml) ;
%         ncwriteatt(ncfile,'Vcrel_mean_ACBL_ml','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_ACBL_ml','description',' Mean of the cloud-relative meridional wind in the active cloud bearing layer for a 100-hPa mean layer parcel ');
%         
%         
%         [al bl] = size(Ucrel_mean_0toACBL_ml);
%         nccreate(ncfile,'Ucrel_mean_0toACBL_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Ucrel_mean_0toACBL_ml',Ucrel_mean_0toACBL_ml) ;
%         ncwriteatt(ncfile,'Ucrel_mean_0toACBL_ml','units','m/s');
%         ncwriteatt(ncfile,'Ucrel_mean_0toACBL_ml','description',' Mean of the cloud-relative zonal wind from the surface to below the active cloud bearing layer for a 100-hPa mean layer parcel ');
%         
%         [al bl] = size(Vcrel_mean_0toACBL_ml);
%         nccreate(ncfile,'Vcrel_mean_0toACBL_ml',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'Vcrel_mean_0toACBL_ml',Vcrel_mean_0toACBL_ml) ;
%         ncwriteatt(ncfile,'Vcrel_mean_0toACBL_ml','units','m/s');
%         ncwriteatt(ncfile,'Vcrel_mean_0toACBL_ml','description',' Mean of the cloud-relative zonal wind from the surface to below the active cloud bearing layer for a 100-hPa mean layer parcel ');
        
        
        [al bl] = size(Ucrel_mean_ACBL_mu);
        nccreate(ncfile,'Ucrel_mean_ACBL_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'Ucrel_mean_ACBL_mu',Ucrel_mean_ACBL_mu) ;
        ncwriteatt(ncfile,'Ucrel_mean_ACBL_mu','units','m/s');
        ncwriteatt(ncfile,'Ucrel_mean_ACBL_mu','description',' Mean of the cloud-relative zonal wind in the active cloud bearing layer for a MU parcel ');
        
        [al bl] = size(Vcrel_mean_ACBL_mu);
        nccreate(ncfile,'Vcrel_mean_ACBL_mu',...
             'Dimensions', {'location',xl,'time',yl,'cell',1},...
            'FillValue','disable');
        ncwrite(ncfile,'Vcrel_mean_ACBL_mu',Vcrel_mean_ACBL_mu) ;
        ncwriteatt(ncfile,'Vcrel_mean_ACBL_mu','units','m/s');
        ncwriteatt(ncfile,'Vcrel_mean_ACBL_mu','description',' Mean of the cloud-relative meridional wind in the active cloud bearing layer for a MU parcel ');
        
        %[al bl] = size(Ucrel_mean_0toACBL_mu);
        %nccreate(ncfile,'Ucrel_mean_0toACBL_mu',...
        %    'Dimensions', {'x',xl,'y',yl,'t',tl,'cell',1},...
        %    'FillValue','disable');
        %ncwrite(ncfile,'Ucrel_mean_0toACBL_mu',Ucrel_mean_0toACBL_mu) ;
        %ncwriteatt(ncfile,'Ucrel_mean_0toACBL_mu','units','m/s');
        %ncwriteatt(ncfile,'Ucrel_mean_0toACBL_mu','description',' Mean of the cloud-relative zonal wind between the initial ht and top of the active cloud bearing layer for a MU parcel ');
        %
        %[al bl] = size(Vcrel_mean_0toACBL_mu);
        %nccreate(ncfile,'Vcrel_mean_0toACBL_mu',...
        %    'Dimensions', {'x',xl,'y',yl,'t',tl,'cell',1},...
        %    'FillValue','disable');
        %ncwrite(ncfile,'Vcrel_mean_0toACBL_mu',Vcrel_mean_0toACBL_mu) ;
        %ncwriteatt(ncfile,'Vcrel_mean_0toACBL_mu','units','m/s');
        %ncwriteatt(ncfile,'Vcrel_mean_0toACBL_mu','description',' Mean of the cloud-relative meridional wind between the initial ht and top of the active cloud bearing layer for a MU parcel ');
        

%         [al bl] = size(U_4000m);
%         nccreate(ncfile,'U_4000m',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'U_4000m',U_4000m) ;
%         ncwriteatt(ncfile,'U_4000m','units','m/s');
%         ncwriteatt(ncfile,'U_4000m','description',' U at 4 km AGL ');
%         
%         
%         [al bl] = size(V_4000m);
%         nccreate(ncfile,'V_4000m',...
%              'Dimensions', {'location',xl,'time',yl,'cell',1},...
%             'FillValue','disable');
%         ncwrite(ncfile,'V_4000m',V_4000m) ;
%         ncwriteatt(ncfile,'V_4000m','units','m/s');
%         ncwriteatt(ncfile,'V_4000m','description',' V at 4 km AGL ');        
        
end
