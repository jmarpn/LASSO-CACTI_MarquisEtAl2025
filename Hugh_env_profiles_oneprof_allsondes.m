
%use this to dumb down the 9 wid 3denv profiels into the aray of 8 ring
%profiles, each done one at a time



% STEP 1, run for each date:
for loc = 2:9


    clearvars -except loc



    hhmmss = ['130000';'130500';'131000';'131500';'132000';'132500';'133000';'133500';'134000';'134500';'135000';'135500';...
        '140000';'140500';'141000';'141500';'142000';'142500';'143000';'143500';'144000';'144500';'145000';'145500';...
        '150000';'150500';'151000';'151500';'152000';'152500';'153000';'153500';'154000';'154500';'155000';'155500';...
        '160000';'160500';'161000';'161500';'162000';'162500';'163000';'163500';'164000';'164500';'165000';'165500';...
        '170000';'170500';'171000';'171500';'172000';'172500';'173000';'173500';'174000';'174500';'175000';'175500';...
        '180000';'180500';'181000';'181500';'182000';'182500';'183000'];

    envdir = '/Users/marq789/Downloads/widenvs/29Nov/';   date = '20181129' ;  ens = '_gefs09_subset.mat'; R=1; d2 = '29Nov100m';
    psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20181129_gefs09_base_d4/';

%     envdir = '/Users/marq789/Downloads/widenvs/04Dec/';   date = '20181204'; ens = 'tindlookback3_gefs19_base_d4.mat'; R=2;
%     psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20181204_gefs19_base_d4/';

%     envdir = '/Users/marq789/Downloads/widenvs/22Jan/';   date = '20190122'; ens = '_gefs18_subset.mat'; R=3;  d2 = '22Jan100m';
%     psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20190122_gefs18_base_d4/';

%     envdir = '/Users/marq789/Downloads/widenvs/23Jan/';   date = '20190123'; ens = '_gefs18_subset.mat'; R=4;  d2 = '23Jan100m';
%     psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20190123_gefs18_base_d4/';

%     envdir = '/Users/marq789/Downloads/widenvs/25Jan/';   date = '20190125'; ens = '_eda07_subset.mat'; R=5;  d2 = '25Jan100m';
%     psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20190125_eda07_base__d4/';

%     envdir = '/Users/marq789/Downloads/widenvs/29Jan/';   date = '20190129'; ens = '_gefs11_subset.mat'; R=6;  d2 = '29Jan100m';
%     psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20190129_gefs11_base_d4/';

%     envdir = '/Users/marq789/Downloads/widenvs/08Feb/';   date = '20190208'; ens = 'tindlookback3_gefs09_base_d4.mat'; R=7;
%     psourcedir = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/20190208_eda08__base_d4/';





    % date =  filedir(end-5:end-1) ;

    envlist = ls([envdir,'wID_3Denv_3R_*']) ;
    envlist = split(envlist) ;
    envlist(end,:) = []

    listnum  = length(envlist) ;


    %%% new vert profile resolution:
    %interp_htasl  =  [250:10:21000] ;
    interp_htasl  =  [250:20:3000, 3050:50:5000, 5250:250:21000 ] ;

    maxwids = 1200;  %this is actaully a little too short for some cases - 29N and 8F had 1214 and 1245 wids respectively.
    numwid  = 2000;  %taking a stab at max number of updrafts ID'ed per case & per time in each case
    alltimes  = 75;  %taking a stab at max num of times per case


    % numwid, time, height %
    oneprof_temper      = zeros(numwid, alltimes, length(interp_htasl));  oneprof_temper(:) = NaN;
    oneprof_qvapor      = zeros(numwid, alltimes, length(interp_htasl));  oneprof_qvapor(:) = NaN;
    oneprof_press       = zeros(numwid, alltimes, length(interp_htasl));  oneprof_press(:) = NaN;
    oneprof_rh          = zeros(numwid, alltimes, length(interp_htasl));  oneprof_rh(:) = NaN;
    oneprof_u           = zeros(numwid, alltimes, length(interp_htasl));  oneprof_u(:) = NaN;
    oneprof_v           = zeros(numwid, alltimes, length(interp_htasl));  oneprof_v(:) = NaN;
    oneprof_htasl       = zeros(numwid, alltimes, length(interp_htasl));  oneprof_htasl(:) = NaN;
    oneprof_terr  = zeros(numwid, alltimes);                              oneprof_terr(:) = NaN;

    oneprof_MUCAPE      = zeros(numwid, alltimes);
    oneprof_MULFCagl    = zeros(numwid, alltimes);
    oneprof_MULCLagl    = zeros(numwid, alltimes);
    oneprof_meanRHacbl  = zeros(numwid, alltimes);
    oneprof_06shear     = zeros(numwid, alltimes);

    % start up read in of already calc'ed sounding params:
    %psourcedir = ['/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/',run,'/'] ;
    wIDparam_list = ls( horzcat(psourcedir,'/WIDlasso_deeper_3R_envparms_*nc') )  ;
    wIDparam_list = split(wIDparam_list) ; wIDparam_list(end) = []  ;
    [ea eb] = size(wIDparam_list); clear eb;


    %loop thru time files for this day
    for TT = 1 : listnum

        %   TT = 20;



        %  load the env profiles
        envfile = char( envlist(TT) ) ;
        load(envfile) ;
        [ba bs bd] = size(WID_3Dtemper);


        % read in & assign wID environments, assigning each the mean around the circle and excluding element 1 (the center pt of the updraft because of contamination):
        MUCAPE = ncread(string(wIDparam_list(TT)),'CAPE_mu');                MUCAPE = permute(MUCAPE,[4 1 2 3]);
        MULFC  = ncread(string(wIDparam_list(TT)),'LFC_height_mu');          MULFC = permute(MULFC,[4 1 2 3]);
        MULCL  = ncread(string(wIDparam_list(TT)),'LCL_height_mu');          MULCL = permute(MULCL,[4 1 2 3]);
        MURHA  = ncread(string(wIDparam_list(TT)),'RH_mean_ACBL');           MURHA = permute(MURHA,[4 1 2 3]);
        S06    = ncread(string(wIDparam_list(TT)),'shear_mag_bulk_0to6km');  S06 = permute(S06,[4 1 2 3]);


        %{
    %tag the profile closest to the median cape for each wid:
    medianloc = zeros(ba,1);
    for w = 1:ba
        % w = 1;
        difmedc = MUCAPE(w,2:9) - median(MUCAPE(w,2:9),'omitnan')  ;
        medc = find( abs(difmedc) ==  min(abs(difmedc))  ) + 1 ;  % the plus one is because you're looking in 2:9 rather than 1:9
        if( isempty(medc) )
            medianloc(w) = NaN;
        else
            medianloc(w) = medc(1);
        end
    end

    %[ca cs] = size(MUCAPE)
        %}

        % regularly look at each sounding location (2-9) in the cloud ring


        for w = 1:ba

            %  w = 200;
            %grab the target sounding

            prof_htasl  = permute(WID_3Dhtasl(w,loc,:),[1 3 2]) ;

            if(  length( find(isnan(prof_htasl)) ) < 5     )
                
                prof_temper = permute(WID_3Dtemper(w,loc,:),[3 1 2]) ;
                prof_qvapor = permute(WID_3Dqvapor(w,loc,:),[3 1 2]) ;
                prof_press  = permute(WID_3Dpress(w,loc,:),[3 1 2]) ;
                prof_rh     = permute(WID_3Drh(w,loc,:),[3 1 2]) ;
                prof_u      = permute(WID_3Du(w,loc,:),[3 1 2]) ;
                prof_v      = permute(WID_3Dv(w,loc,:),[3 1 2]) ;

                interped_htasl = interp_htasl ;
                interp_temper = interp1( prof_htasl, prof_temper, interp_htasl) ;
                interp_qvapor = interp1( prof_htasl, prof_qvapor, interp_htasl) ;
                interp_press  = interp1( prof_htasl, prof_press, interp_htasl) ;
                interp_rh     = interp1( prof_htasl, prof_rh, interp_htasl) ;
                interp_u      = interp1( prof_htasl, prof_u, interp_htasl) ;
                interp_v      = interp1( prof_htasl, prof_v, interp_htasl) ;

                prof_MUCAPE      =  MUCAPE(w,loc);
                prof_MULFCagl    =  MULFC(w,loc);
                prof_MULCLagl    =  MULCL(w,loc);
                prof_meanRHacbl  =  MURHA(w,loc);
                prof_06shear     =  S06(w,loc);

            else

                prof_htasl  = zeros(bd,1); prof_htasl(:) = NaN;
                prof_temper = zeros(bd,1); prof_temper(:) = NaN;
                prof_qvapor = zeros(bd,1); prof_qvapor(:) = NaN;
                prof_rh     = zeros(bd,1); prof_rh(:) = NaN;
                prof_press  = zeros(bd,1); prof_press(:) = NaN;
                prof_u      = zeros(bd,1); prof_u(:) = NaN;
                prof_v      = zeros(bd,1); prof_v(:) = NaN;

                interped_htasl  = zeros(length(interp_htasl),1);   interped_htasl(:) = NaN;
                interp_temper = zeros(length(interp_htasl),1);   interp_temper(:) = NaN;
                interp_qvapor = zeros(length(interp_htasl),1);   interp_qvapor(:) = NaN;
                interp_rh     = zeros(length(interp_htasl),1);   interp_rh(:) = NaN;
                interp_press  = zeros(length(interp_htasl),1);   interp_press(:) = NaN;
                interp_u      = zeros(length(interp_htasl),1);   interp_u(:) = NaN;
                interp_v      = zeros(length(interp_htasl),1);   interp_v(:) = NaN;

                prof_MUCAPE      =  NaN;
                prof_MULFCagl    =  NaN;
                prof_MULCLagl    =  NaN;
                prof_meanRHacbl  =  NaN;
                prof_06shear     =  NaN;

            end

            %log your targetd interped profile here:
            oneprof_temper(w,TT,:)       =   interp_temper' ;
            oneprof_qvapor(w,TT,:)       =   interp_qvapor' ;
            oneprof_press(w,TT,:)        =   interp_press' ;
            oneprof_rh(w,TT,:)           =   interp_rh' ;
            oneprof_u(w,TT,:)            =   interp_u' ;
            oneprof_v(w,TT,:)            =   interp_v' ;
            oneprof_htasl(w,TT,:)        =   interped_htasl' ;
            if( isnan(loc)==0 )
                oneprof_terr(w,TT)           = WID_3Dterr(w,loc);
            else
                oneprof_terr(w,TT)           = NaN;
            end

            oneprof_MUCAPE(w,TT)      = prof_MUCAPE;
            oneprof_MULFCagl(w,TT)    = prof_MULFCagl;
            oneprof_MULCLagl(w,TT)    = prof_MULCLagl;
            oneprof_meanRHacbl(w,TT)  = prof_meanRHacbl;
            oneprof_06shear(w,TT)     = prof_06shear;


        end  %w

        %     clear WID_3Dtemper WID_3Dqvapor WID_3Dpress WID_3Drh WID_3Du WID_3Dv
        %     clear interp_temper interp_qvapor interp_press interp_rh interp_u interp_v
        %     clear prof_temper prof_qvapor prof_press prof_rh prof_u prof_v

    end %t

    % blah = oneprof_htasl(:,:,250) ;
    % blah = oneprof_temper(:,:,250) ;
    % blah = oneprof_qvapor(:,:,250) ;
    % blah = oneprof_press(:,:,250) ;
    % blah = oneprof_rh(:,:,250) ;
    % blah = oneprof_u(:,:,150) ;
    % blah = oneprof_v(:,:,250) ;
    % blah = oneprof_terr(:,:) ;



    %define agl frame:
    oneprof_htagl = oneprof_htasl;
    %oneprof_pertime_htagl = oneprof_pertime_htasl;

    for TT = 1 : listnum

        %kill the dirt in htagl profile
        for w = 1:ba
            % w = 100;  TT = 20;
            %find first k index with data to call the ground
            profz = find( isnan( oneprof_qvapor(w,TT,:) )==0 )  ;
            if( isempty( profz)==0  )
                %find first k index with data to call the ground
                oneprof_htagl(w,TT,:) = oneprof_htasl(w,TT,:) -  oneprof_htasl(w,TT,profz(1)) ;
            end
        end

        % %kill the dirt in htagl profile
        % for w = 1:ba
        %     % w = 100;  TT = 20;
        %     %find first k index with data to call the ground
        %     profz = find( isnan( oneprof_pertime_qvapor(w,1,:) )==0 )  ;
        %     if( isempty( profz)==0  )
        %         %find first k index with data to call the ground
        %         oneprof_pertime_htagl(w,1,:) = oneprof_pertime_htasl(w,1,:) -  oneprof_pertime_htasl(w,1,profz(1)) ;
        %     end
        % end
        %
        % oneprof_pertime_htagl(oneprof_pertime_htagl<0) = NaN;
        %
        % oneprof_pertime_temper((maxwids+1):end,:,:)     =   [];
        % oneprof_pertime_qvapor((maxwids+1):end,:,:)     =   [];
        % oneprof_pertime_press((maxwids+1):end,:,:)      =   [];
        % oneprof_pertime_rh((maxwids+1):end,:,:)         =   [];
        % oneprof_pertime_u((maxwids+1):end,:,:)          =   [];
        % oneprof_pertime_v((maxwids+1):end,:,:)          =   [];
        % oneprof_pertime_htasl((maxwids+1):end,:,:)      =   [];
        % oneprof_pertime_htagl((maxwids+1):end,:,:)      =   [];
        % oneprof_pertime_centerterr((maxwids+1):end,:)      =   [];


        %     if(R==2 | R==7)
        %         save([filedir,'WID_OneProfHugh_aslmean_',date,'_',hhmmss(TT,:),'_',ens],'oneprof_pertime_temper','oneprof_pertime_qvapor',...
        %             'oneprof_pertime_press','oneprof_pertime_rh','oneprof_pertime_u','oneprof_pertime_v',...
        %             'oneprof_pertime_htasl','oneprof_pertime_htagl','oneprof_pertime_centerterr')
        %     else
        %         save([filedir,'WID_OneProfHugh_aslmean_',d2,'_',hhmmss(TT,:),'_tindlookback3_',date,ens],'oneprof_pertime_temper','oneprof_pertime_qvapor',...
        %             'oneprof_pertime_press','oneprof_pertime_rh','oneprof_pertime_u','oneprof_pertime_v',...
        %             'oneprof_pertime_htasl','oneprof_pertime_htagl','oneprof_pertime_centerterr')
        %     end
        %
        %



    end%t


    % blah = oneprof_htagl(:,:,150) ;


    oneprof_htagl(oneprof_htagl<0) = NaN;

    %%%%% prune off to 1:1200 wids - even though this kills a few samples in 29N and 8F
    oneprof_temper((maxwids+1):end,:,:)     =   [];
    oneprof_qvapor((maxwids+1):end,:,:)     =   [];
    oneprof_press((maxwids+1):end,:,:)      =   [];
    oneprof_rh((maxwids+1):end,:,:)         =   [];
    oneprof_u((maxwids+1):end,:,:)          =   [];
    oneprof_v((maxwids+1):end,:,:)          =   [];
    oneprof_htasl((maxwids+1):end,:,:)      =   [];
    oneprof_htagl((maxwids+1):end,:,:)      =   [];
    oneprof_terr((maxwids+1):end,:)         =   [];

    oneprof_MUCAPE((maxwids+1):end,:)       =   [];
    oneprof_MULFCagl((maxwids+1):end,:)     =   [];
    oneprof_MULCLagl((maxwids+1):end,:)     =   [];
    oneprof_meanRHacbl((maxwids+1):end,:)   =   [];
    oneprof_06shear((maxwids+1):end,:)      =   [];



    % blah0 = oneprof_htagl(:,:,200) ;
    % blah1 = oneprof_htagl(79,15,:) ;
    % blah2 = oneprof_htasl(79,15,:) ;
    % blah3 = oneprof_qvapor(79,15,:) ;
    % blah4 = oneprof_press(79,15,:) ;



    save([envdir,'WID_OneProf_3R_',date,'_loc',num2str(loc),'.mat'],'oneprof_temper','oneprof_qvapor',...
        'oneprof_press','oneprof_rh','oneprof_u','oneprof_v','oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear','loc')


end %loc loop

%%% END OF STEP 1








%%%%%%%%%%%%%%%%%%%%%%%% START STEP 2


for loc = 2:9

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%    ALLCASES glued together
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %  loc = 9

    clearvars -except loc

    envdir = '/Users/marq789/Downloads/widenvs/';

    wID_T = [];
    wID_QV = [];
    wID_P = [];
    wID_RH = [];
    wID_U = [];
    wID_V = [];
    wID_HTasl = [];
    wID_HTagl = [];
    wID_terr = [];
    wID_CAPE = [];
    wID_LFCagl = [];
    wID_LCLagl = [];
    wID_Sh06 = [];
    wID_RHa = [];

    filein = ['/Users/marq789/Downloads/widenvs/29Nov/WID_OneProf_3R_20181129_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    filein = ['/Users/marq789/Downloads/widenvs/04Dec/WID_OneProf_3R_20181204_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    filein = ['/Users/marq789/Downloads/widenvs/22Jan/WID_OneProf_3R_20190122_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    filein = ['/Users/marq789/Downloads/widenvs/23Jan/WID_OneProf_3R_20190123_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    filein = ['/Users/marq789/Downloads/widenvs/25Jan/WID_OneProf_3R_20190125_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    filein = ['/Users/marq789/Downloads/widenvs/29Jan/WID_OneProf_3R_20190129_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    filein = ['/Users/marq789/Downloads/widenvs/08Feb/WID_OneProf_3R_20190208_loc',num2str(loc),'.mat']
    load(filein,...
        'oneprof_temper','oneprof_qvapor','oneprof_press','oneprof_rh','oneprof_u','oneprof_v',...
        'oneprof_htasl','oneprof_htagl','oneprof_terr',...
        'oneprof_MUCAPE','oneprof_MULFCagl','oneprof_MULCLagl','oneprof_meanRHacbl','oneprof_06shear')
    wID_T       = cat(4,wID_T,oneprof_temper);
    wID_QV      = cat(4,wID_QV,oneprof_qvapor);
    wID_P       = cat(4,wID_P,oneprof_press);
    wID_RH      = cat(4,wID_RH,oneprof_rh);
    wID_U       = cat(4,wID_U,oneprof_u);
    wID_V       = cat(4,wID_V,oneprof_v);
    wID_HTasl   = cat(4,wID_HTasl,oneprof_htasl);
    wID_HTagl   = cat(4,wID_HTagl,oneprof_htagl);
    wID_terr    = cat(3,wID_terr,oneprof_terr);
    wID_CAPE    = cat(3,wID_CAPE,oneprof_MUCAPE);
    wID_LFCagl  = cat(3,wID_LFCagl,oneprof_MULFCagl);
    wID_LCLagl  = cat(3,wID_LCLagl,oneprof_MULCLagl);
    wID_Sh06    = cat(3,wID_Sh06,oneprof_06shear);
    wID_RHa     = cat(3,wID_RHa,oneprof_meanRHacbl);
    clear oneprof_temper oneprof_qvapor oneprof_press oneprof_rh oneprof_u oneprof_v oneprof_htasl oneprof_htagl oneprof_terr
    clear oneprof_MUCAPE oneprof_MULFCagl oneprof_MULCLagl oneprof_meanRHacbl oneprof_06shear


    wID_T = permute(wID_T,[3 2 1 4]);
    wID_QV = permute(wID_QV,[3 2 1 4]);
    wID_P = permute(wID_P,[3 2 1 4]);
    wID_RH = permute(wID_RH,[3 2 1 4]);
    wID_U = permute(wID_U,[3 2 1 4]);
    wID_V = permute(wID_V,[3 2 1 4]);
    wID_HTasl = permute(wID_HTasl,[3 2 1 4]);
    wID_HTagl = permute(wID_HTagl,[3 2 1 4]);

    wID_terr = permute(wID_terr,[2 1 3]);
    wID_LFCagl = permute(wID_LFCagl,[2 1 3]);
    wID_LCLagl = permute(wID_LCLagl,[2 1 3]);
    wID_CAPE = permute(wID_CAPE,[2 1 3]);
    wID_RHa = permute(wID_RHa,[2 1 3]);
    wID_Sh06 = permute(wID_Sh06,[2 1 3]);


    [az at aw ar] = size(wID_RH) ;
    %%%%%%%%%%%%% load the wID's



    load('/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/processed_ALLCASES/SwathIDthermoOutput_c1p11h1p13t_wraw_wmag025_env3R_v12b_SDCsmj1p5_topf.mat',...
        'widLFCAR_profs_maxmag_ALLCASES','widLFCAR_profs_area_ALLCASES','wid_maxk_topasl_ALLCASES','dist_wid_and_swath_ALLCASES')

    CLDupdraft_area       =  max(widLFCAR_profs_area_ALLCASES(:,:,:,:),[],2) ;   %
    CLDupdraft_area = permute(CLDupdraft_area ,[1 3 4 2]);
    CLDupdraft_Radius     =  ((CLDupdraft_area/3.14159).^0.5)/1000;      %circ equiv rad [km]


    % CLDupdraft_MaxMagnitude = max(widLFCAR_profs_maxmag_ALLCASES(:,:,:,:),[],2,'omitnan') ;
    %     CLDupdraft_MaxMagnitude = permute(CLDupdraft_MaxMagnitude,[1 3 4 2]);

    CLDupdraft_topZasl    =  wid_maxk_topasl_ALLCASES(:,:,:);


    %filter on distance from swath
    TIMLG = 4;
    distswathwid_thresh = 3.0 ;
    kill = find(  abs( dist_wid_and_swath_ALLCASES(:,:,TIMLG,:) )  >  distswathwid_thresh  );
    CLDupdraft_area(kill) = NaN;
    %CLDupdraft_MaxMagnitude(kill) = NaN;
    CLDupdraft_Radius(kill) = NaN;
    CLDupdraft_topZasl(kill) = NaN;

    %mask to deal with the LFCAR frame not yet handled in the cloud top field,
    %it now has dist from swath info in it too:
    cht_mask = CLDupdraft_area;   cht_mask( isnan(CLDupdraft_area)==0 )=1;


    CLDupdraft_topZasl = CLDupdraft_topZasl .* cht_mask ;


    CLDupdraft_area    = permute(CLDupdraft_area,[2 1 3]);
    CLDupdraft_Radius  = permute(CLDupdraft_Radius,[2 1 3]);
    CLDupdraft_topZasl = permute(CLDupdraft_topZasl,[2 1 3]);

    %make mask of correct dimensions:
    CHT_MASK = permute(cht_mask,[2 1 3]) ;
    CHT_MASK_4D = zeros(at, aw, ar, az);
    for k = 1:az
        CHT_MASK_4D(:,:,:,k) = CHT_MASK;
    end
    CHT_MASK_4D = permute(CHT_MASK_4D,[4 1 2 3]);

    % hit all with joint wid filters:
    wID_terr_filt   = wID_terr   .* CHT_MASK;
    wID_LFCagl_filt = wID_LFCagl .* CHT_MASK;
    wID_LCLagl_filt = wID_LCLagl .* CHT_MASK;
    wID_CAPE_filt   = wID_CAPE   .* CHT_MASK;
    wID_RHa_filt    = wID_RHa    .* CHT_MASK;
    wID_Sh06_filt   = wID_Sh06   .* CHT_MASK;

    wID_LFCasl_filt = wID_LFCagl_filt + wID_terr_filt ;
    wID_LCLasl_filt = wID_LCLagl_filt + wID_terr_filt ;

    wID_T_filt       = wID_T .* CHT_MASK_4D ;
    wID_QV_filt      = wID_QV .* CHT_MASK_4D ;
    wID_P_filt       = wID_P .* CHT_MASK_4D ;
    wID_RH_filt      = wID_RH .* CHT_MASK_4D ;
    wID_U_filt       = wID_U .* CHT_MASK_4D ;
    wID_V_filt       = wID_V .* CHT_MASK_4D ;
    wID_HTasl_filt   = wID_HTasl .* CHT_MASK_4D ;
    wID_HTagl_filt   = wID_HTagl .* CHT_MASK_4D ;





    %{

    length(find(isnan( CLDupdraft_area )==0))

    length(find(isnan( wID_LFCagl_filt )==0))

    %}





    %%% up to here are full profiles in Htasl/agl frame
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   reorient to LFC-relative frame for Hugh code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    [ pz pt pw pr ] = size( wID_HTasl_filt );

    lfcrel_htasl = wID_HTasl_filt ;
    for t = 1:pt
        for w = 1:pw
            for r = 1:pr
                lfcrel_htasl(:,t,w,r) = wID_HTasl_filt(:,t,w,r) - wID_LCLasl_filt(t,w,r);
            end
        end
    end

    %reinterp to regular lfc-rel ht:

    Zasl_lfcrel = [0:100:18000] ;
    Temp_lfcrel = zeros(length(Zasl_lfcrel),pt,pw,pr);
    Qvap_lfcrel = zeros(length(Zasl_lfcrel),pt,pw,pr);
    P_lfcrel    = zeros(length(Zasl_lfcrel),pt,pw,pr);
    U_lfcrel    = zeros(length(Zasl_lfcrel),pt,pw,pr);
    V_lfcrel    = zeros(length(Zasl_lfcrel),pt,pw,pr);

    for t = 1:pt
        for w = 1:pw
            for r = 1:pr

                %  t=5; w=9; r=1;

                lfcrel_htasl = wID_HTasl_filt(:,t,w,r) - wID_LFCasl_filt(t,w,r);

                if(     length(find(isnan( wID_T_filt(:,t,w,r) )))  <  length( wID_T_filt(:,t,w,r) )   &  ...
                           length(find(isnan(lfcrel_htasl))) <  length(find(lfcrel_htasl))    )   %there can be cases without an LFC (surprisingly)

%                     %lfcrel_htasl(:,t,w,r) = wID_HTasl_filt(:,t,w,r) - wID_LCLasl_filt(t,w,r);
%                     lfcrel_htasl = wID_HTasl_filt(:,t,w,r) - wID_LFCasl_filt(t,w,r);

                    Temp_lfcrel(:,t,w,r) = interp1( lfcrel_htasl, wID_T_filt(:,t,w,r),  Zasl_lfcrel ) ;
                    Qvap_lfcrel(:,t,w,r) = interp1( lfcrel_htasl, wID_QV_filt(:,t,w,r), Zasl_lfcrel ) ;
                    P_lfcrel(:,t,w,r) = interp1( lfcrel_htasl, wID_P_filt(:,t,w,r), Zasl_lfcrel ) ;
                    U_lfcrel(:,t,w,r) = interp1( lfcrel_htasl, wID_U_filt(:,t,w,r), Zasl_lfcrel ) ;
                    V_lfcrel(:,t,w,r) = interp1( lfcrel_htasl, wID_V_filt(:,t,w,r), Zasl_lfcrel ) ;

                else

                    Temp_lfcrel(:,t,w,r) = NaN ;
                    Qvap_lfcrel(:,t,w,r) = NaN ;
                    P_lfcrel(:,t,w,r)    = NaN ;
                    U_lfcrel(:,t,w,r)    = NaN ;
                    V_lfcrel(:,t,w,r)    = NaN ;

                end

            end
        end
    end


    % t=5; w=9; r=1;
    % t=17; w=13; r=1;

    % wID_LCLasl_filt(:,:,1)

    % Temp_lfcrel(:,t,w,r)
    % Qvap_lfcrel(:,t,w,r)*1000
    % P_lfcrel(:,t,w,r)




    %%%%%%%%%%%% condense by culling the nan profiles:

    culled_T  = zeros(length(Zasl_lfcrel),75,100,7) ;  culled_T(:) = NaN;
    culled_QV = zeros(length(Zasl_lfcrel),75,100,7) ;  culled_QV(:) = NaN;
    culled_P  = zeros(length(Zasl_lfcrel),75,100,7) ;  culled_P(:) = NaN;
    culled_U  = zeros(length(Zasl_lfcrel),75,100,7) ;  culled_U(:) = NaN;
    culled_V  = zeros(length(Zasl_lfcrel),75,100,7) ;  culled_V(:) = NaN;

    culled_MUCAPE = zeros(75,100,7) ;  culled_MUCAPE(:) = NaN;
    culled_LCLagl = zeros(75,100,7) ;  culled_LCLagl(:) = NaN;
    culled_LCLasl = zeros(75,100,7) ;  culled_LCLasl(:) = NaN;
    culled_LFCagl = zeros(75,100,7) ;  culled_LCLagl(:) = NaN;
    culled_LFCasl = zeros(75,100,7) ;  culled_LFCasl(:) = NaN;
    culled_Sh06   = zeros(75,100,7) ;  culled_Sh06(:) = NaN;
    culled_RHa    = zeros(75,100,7) ;  culled_RHa(:) = NaN;

    culled_CLDupdraft_area       = zeros(75,100,7) ;  culled_CLDupdraft_area(:) = NaN;
    culled_CLDupdraft_topZasl    = zeros(75,100,7) ;  culled_CLDupdraft_topZasl(:) = NaN;

    for R = 1:7
        for t = 1:75

            % R = 1; t = 1;

            keepw = find(  isnan(     wID_LCLasl_filt(t,:,R)     )==0  ) ;
            for c = 1:length(keepw)
                culled_MUCAPE(t,c,R)  =  wID_CAPE_filt( t,keepw(c),R ) ;
                culled_LCLagl(t,c,R)  =  wID_LCLagl_filt( t,keepw(c),R ) ;
                culled_LCLasl(t,c,R)  =  wID_LCLasl_filt( t,keepw(c),R ) ;
                culled_LFCagl(t,c,R)  =  wID_LFCagl_filt( t,keepw(c),R ) ;
                culled_LFCasl(t,c,R)  =  wID_LFCasl_filt( t,keepw(c),R ) ;
                culled_Sh06(t,c,R)    =  wID_Sh06_filt( t,keepw(c),R ) ;
                culled_RHa(t,c,R)     =  wID_RHa_filt( t,keepw(c),R ) ;

                culled_CLDupdraft_area(t,c,R)        =  CLDupdraft_area( t, keepw(c), R ) ;
                culled_CLDupdraft_topZasl(t,c,R)     =  CLDupdraft_topZasl( t, keepw(c), R ) ;

                culled_T(:,t,c,R)     =  Temp_lfcrel(:,t,keepw(c),R ) ;
                culled_QV(:,t,c,R)    =  Qvap_lfcrel(:,t,keepw(c),R ) ;
                culled_P(:,t,c,R)     =  P_lfcrel(:,t,keepw(c),R ) ;
                culled_U(:,t,c,R)     =  U_lfcrel(:,t,keepw(c),R ) ;
                culled_V(:,t,c,R)     =  V_lfcrel(:,t,keepw(c),R ) ;

            end
        end
    end


    save([envdir,'WID_OneProf_3R_alldates_loc', num2str(loc) ,'.mat'],...
        'culled_T','culled_QV','culled_P','culled_U','culled_V','Zasl_lfcrel',...
        ...
        'culled_CLDupdraft_area','culled_CLDupdraft_topZasl',...
        ...
        'culled_MUCAPE','culled_LCLagl','culled_LCLasl','culled_LFCagl','culled_LFCasl','culled_Sh06','culled_RHa',...
        ...
        'loc')


end %loc loop




% clear 
% load('/Users/marq789/Downloads/widenvs/WID_OneProf_3R_alldates_loc3.mat')




%%%%%%%%%%%%%%%%%%%%%% STEP 3: append all locs


envdir = '/Users/marq789/Downloads/widenvs/';

updraft_area    = [];
updraft_topZasl = [];
% %updraft_terrht  = wID_terr_filt ;
updraft_LFCagl = [];
updraft_LCLagl = [];
updraft_LFCasl = [];
updraft_LCLasl = [];
updraft_MUCAPE = [];
updraft_acblRH = [];
updraft_Shear06 = [];
prof_T = [];
prof_Qvap = [];
prof_P = [];
prof_Z = [];
prof_U = [];
prof_V = [];

for loc = 2:9


    load([envdir,'WID_OneProf_3R_alldates_loc', num2str(loc) ,'.mat'])

    updraft_area    = culled_CLDupdraft_area;
    updraft_topZasl = culled_CLDupdraft_topZasl;

    updraft_LFCagl = cat(4,updraft_LFCagl,culled_LFCagl);
    updraft_LCLagl = cat(4,updraft_LCLagl,culled_LCLagl);
    updraft_LFCasl = cat(4,updraft_LFCasl,culled_LFCasl);
    updraft_LCLasl = cat(4,updraft_LCLasl,culled_LCLasl);
    updraft_MUCAPE = cat(4,updraft_MUCAPE,culled_MUCAPE);
    updraft_acblRH = cat(4,updraft_acblRH,culled_RHa);
    updraft_Shear06 = cat(4,updraft_Shear06,culled_Sh06);
   
    prof_T      = cat(5,prof_T,culled_T);
    prof_Qvap   = cat(5,prof_Qvap,culled_QV);
    prof_P      = cat(5,prof_P,culled_P);
    prof_Z      = Zasl_lfcrel;
    prof_U      = cat(5,prof_U,culled_U);
    prof_V      = cat(5,prof_V,culled_V);
    % %updraft_terrht  = wID_terr_filt ;



end






%{    
% zl tl nl rl ll
up = 10
tt = 20
    figure;   plot(prof_Z,'k')
    figure;   plot(   prof_Qvap(:,tt,up,1,1) ,'k')
%}







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FINAL STEP
% write everything to netcdf
% write everything to netcdf
% write everything to netcdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncfile  = ['/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/processed_ALLCASES/ProfilesForHugh_8sites.nc' ] ;
delete(ncfile)

[tl nl rl] = size(updraft_area) ;
nccreate(ncfile,'updraft_area',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_area',updraft_area) ;
ncwriteatt(ncfile,'updraft_area','units','m^2');
ncwriteatt(ncfile,'updraft_area','description','maximum horizontal updraft area at a height within 2km-deep layer whose bottom bound is LFC');

% [tl nl rl] = size(CLDupdraft_Radius) ;
% nccreate(ncfile,'CLDupdraft_Radius',...
%     'Dimensions', {'time',tl,'updraft',nl,'case',rl},...
%     'FillValue','disable');
% ncwrite(ncfile,'CLDupdraft_Radius',CLDupdraft_Radius) ;
% ncwriteatt(ncfile,'CLDupdraft_Radius','units','km');
% ncwriteatt(ncfile,'CLDupdraft_Radius','description','circular equiv radius of CLDupdraft_area');


[tl nl rl] = size(updraft_topZasl)
nccreate(ncfile,'updraft_topZasl',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_topZasl',updraft_topZasl) ;
ncwriteatt(ncfile,'updraft_topZasl','units','m asl');
ncwriteatt(ncfile,'updraft_topZasl','description','Top height of the cloudy updraft');



[tl nl rl ll] = size(updraft_MUCAPE) ;
nccreate(ncfile,'updraft_MUCAPE',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_MUCAPE',updraft_MUCAPE) ;
ncwriteatt(ncfile,'updraft_MUCAPE','units','J/kg');
ncwriteatt(ncfile,'updraft_MUCAPE','description','MU CAPE (deep layer, reversible, no ice, no entrainment assumed)');

[tl nl rl ll] = size(updraft_acblRH) ;
nccreate(ncfile,'updraft_acblRH',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_acblRH',updraft_acblRH) ;
ncwriteatt(ncfile,'updraft_acblRH','units','0-100%');
ncwriteatt(ncfile,'updraft_acblRH','description','mean RH in vert layer between MU LFC ht and 1.5km above it');

% [tl nl rl] = size(updraft_Shear06) ;
% nccreate(ncfile,'updraft_Shear06',...
%     'Dimensions', {'time',tl,'updraft',nl,'case',rl},...
%     'FillValue','disable');
% ncwrite(ncfile,'updraft_Shear06',updraft_Shear06) ;
% ncwriteatt(ncfile,'updraft_Shear06','units','m/s');
% ncwriteatt(ncfile,'updraft_Shear06','description','Bulk wind difference between 6km agl and 0 km');

[tl nl rl ll] = size(updraft_LFCagl) ;
nccreate(ncfile,'updraft_LFCagl',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_LFCagl',updraft_LFCagl) ;
ncwriteatt(ncfile,'updraft_LFCagl','units','m agl');
ncwriteatt(ncfile,'updraft_LFCagl','description','MU LFC height above ground level');

[tl nl rl ll] = size(updraft_LCLagl) ;
nccreate(ncfile,'updraft_LCLagl',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_LCLagl',updraft_LCLagl) ;
ncwriteatt(ncfile,'updraft_LCLagl','units','m agl');
ncwriteatt(ncfile,'updraft_LCLagl','description','MU LCL height above ground level');

[tl nl rl ll] = size(updraft_LFCasl) ;
nccreate(ncfile,'updraft_LFCasl',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_LFCasl',updraft_LFCasl) ;
ncwriteatt(ncfile,'updraft_LFCasl','units','m asl');
ncwriteatt(ncfile,'updraft_LFCasl','description','MU LFC height above sea level');

[tl nl rl ll] = size(updraft_LCLasl) ;
nccreate(ncfile,'updraft_LCLasl',...
    'Dimensions', {'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'updraft_LCLasl',updraft_LCLasl) ;
ncwriteatt(ncfile,'updraft_LCLasl','units','m asl');
ncwriteatt(ncfile,'updraft_LCLasl','description','MU LCL height above sea level');


% prof_T      = culled_T;
% prof_Qvap   = culled_QV;
% prof_P      = culled_P;
% prof_Z      = Zasl_lfcrel;
% prof_U      = culled_U;
% prof_V      = culled_V;


[zl tl nl rl ll] = size(prof_T) ;
nccreate(ncfile,'prof_T',...
    'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'prof_T',prof_T) ;
ncwriteatt(ncfile,'prof_T','units','K');
ncwriteatt(ncfile,'prof_T','description','Temperature profile above LFC ht per updraft');

[zl tl nl rl ll] = size(prof_Qvap) ;
nccreate(ncfile,'prof_Qvap',...
    'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'prof_Qvap',prof_Qvap) ;
ncwriteatt(ncfile,'prof_Qvap','units','kg/kg');
ncwriteatt(ncfile,'prof_Qvap','description','Vapor mixing ratio profile above LFC ht per updraft');

[zl tl nl rl ll] = size(prof_P) ;
nccreate(ncfile,'prof_P',...
    'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'prof_P',prof_P) ;
ncwriteatt(ncfile,'prof_P','units','hPa');
ncwriteatt(ncfile,'prof_P','description','Pressure profile above LFC ht per updraft');

% [zl tl nl rl] = size(prof_RH) ;
% nccreate(ncfile,'prof_RH',...
%     'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl},...
%     'FillValue','disable');
% ncwrite(ncfile,'prof_RH',prof_RH) ;
% ncwriteatt(ncfile,'prof_RH','units','0-100%');
% ncwriteatt(ncfile,'prof_RH','description','Relative humidity profile per updraft');

[zl tl nl rl ll] = size(prof_U) ;
nccreate(ncfile,'prof_U',...
    'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'prof_U',prof_U) ;
ncwriteatt(ncfile,'prof_U','units','m/s');
ncwriteatt(ncfile,'prof_U','description','Zonal component of wind profile per updraft');

[zl tl nl rl ll] = size(prof_V) ;
nccreate(ncfile,'prof_V',...
    'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl,'loc',ll},...
    'FillValue','disable');
ncwrite(ncfile,'prof_V',prof_V) ;
ncwriteatt(ncfile,'prof_V','units','m/s');
ncwriteatt(ncfile,'prof_V','description','Meridional component of wind profile per updraft');

zl = length(prof_Z')
nccreate(ncfile,'prof_Z',...
    'Dimensions', {'height',zl},...
    'FillValue','disable');
ncwrite(ncfile,'prof_Z',prof_Z) ;
ncwriteatt(ncfile,'prof_Z','units','m ');
ncwriteatt(ncfile,'prof_Z','description','Altitude profile above LFC ht per updraft (mental note: originally calculated from from asl frame)');

% [zl tl nl rl] = size(prof_Zagl) ;
% nccreate(ncfile,'prof_Zagl',...
%     'Dimensions', {'height',zl,'time',tl,'updraft',nl,'case',rl},...
%     'FillValue','disable');
% ncwrite(ncfile,'prof_Zagl',prof_Zagl) ;
% ncwriteatt(ncfile,'prof_Zagl','units','m agl');
% ncwriteatt(ncfile,'prof_Zagl','description','Height above ground lvl profile per updraft');











%post diagnotiscs:
file = '/Users/marq789/Documents/PROJECTS/ICLASS/LASSO_stuff/wID_swathID_files/processed_ALLCASES/ProfilesForHugh_8sites.nc';


ncdisp(file)







