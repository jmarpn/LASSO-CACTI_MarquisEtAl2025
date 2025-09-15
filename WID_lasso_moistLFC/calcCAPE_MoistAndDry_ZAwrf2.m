function [moiCAPE, dryCAPE, CIN, moiLFCh, moiELh, CAPEacbl, ELtooshort] = calcCAPE_MoistAndDry_ZAwrf2(envT,lifted_parcel_T,h,Zo,LCLi,SatLog);


%%%%% input
%  envT = your chosen environmental temperature profile (suggest either Tv or Trho - granted, env. Trho probably = Tv, so Tv is prob fine; K)
%  lifted_parcel_T = your chosen lifted parcel temperature profile (Tv or Trho)
%  h = height above ground level of sounding (m)
%  Zo = height index to lift parcel from (k)
%  SatLog = binary profile of 0 (lifted parcel = dry) or 1 (lifted parce =
%       saturated) - JNM added this is output from John's lifted parcel
%       algorithm

%%%%% output
% CAPE      -  (J/kg)
% LFCh      -  lfc height (m)
% ELh       -  equilibrium level (m)
% negsum    -  just negative area CIN (l/kg)
% CAPEacbl  -  CAPE in just the ACBL (J/kg)
% cinrat    -  CIN 'depth', i.e., fraction of # of vertical levels below LFC with negative area.

%strip%        LFCp = NaN;
%strip%        ELp = NaN;
%strip%        CAPElcl = NaN;


%   envT = Tv;  lifted_parcel_T = revT_rho_lif_sfc;  h = h;  Zo = Zo;  LCLi = sfcLCLi;  SatLog = IsSat_JNM; 

%   envT = Tv;  lifted_parcel_T = revT_rho_lif_mu;  h = h;  Zo = Zo;  LCLi = muLCLi;  SatLog = IsSat_JNM; 

%   envT = Tv_4ml;  lifted_parcel_T = revT_rho_lif_ml;  h = h;  Zo = Zo;  LCLi = mlLCLi;  SatLog = IsSat_JNM; 



%hacks to adapt to old code that was hardcoded to Tv, but not Trho might be the better variable(?), so you can set it here instead of Tv. Clean this up later:
Tv = envT;
parcelTvprof = lifted_parcel_T;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% now determine EL, LFC:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



diffTv = Tv; diffTv(:) = NaN;
for k = length(Tv):-1:1
    diffTv(k) = parcelTvprof(k)'-Tv(k);
end



% %f/g-correction
% diffTv(Zo) = 0;  %had to set this because for some reason, I found an example where T'(Zo) was not 0.0
%
%
% if(diffTv(end) > 0) % sounding temp profile is not deep enough to find EL
%     %tallenough = vertcat(tallenough,0);
%     tallenough = 0 ;
% elseif(diffTv(end) < 0) % sounding temp profile is deep enough to find EL
%     %tallenough = vertcat(tallenough,1);
%     tallenough = 1;
% end

%strip% LCLp = p(topdry)/100 ;   %mb
%LCLh = h(topdry);        %meters (agl)

% define full (dry & moist) B profile for integration:
B = 9.8.*(diffTv)./Tv;
B(  find( isnan(SatLog)  )  ) = NaN ;
%now break into dry and moist portions using SatLog :
Bdry = B;
Bdry(SatLog == 1  |  isnan(SatLog) ) = 0;
Bmoist = B;
Bmoist(SatLog == 0  |  isnan(SatLog) ) = 0;



ELtooshort = 1;
%find moist LFC (lowest point where Bmoist > 0):
PosBmoist = find( Bmoist > 0) ;
%set lfc to bottom Bmoist>0 height list (if there is Bmoist > 0)

if( isempty(PosBmoist)==0  & length(PosBmoist)>1 )  %if there is any even thin layer of moist Buoyancy

    %final answer for LFC:
    k_moiLFC = PosBmoist(1);  %lfc is lowest point in column
    moiLFCh = h(k_moiLFC);
    %moiLFCp = p(k_moiLFC)/100

    %find locations in profile where B full "flips" sign, which you will
    %use later for EL and CAPE.
    diffsignT = sign(Bmoist);  %info on just the + or - sign of profile
    flipsTmoi = [];
    for k = Zo+1:1:length(diffTv)   %g-correction
        if(diffsignT(k)-diffsignT(k-1) ~= 0)
            flipsTmoi = vertcat(flipsTmoi,k);
        end
    end
    %get rid of sub-LFC flips (if Bmoist < 0)
    kill = find(flipsTmoi < k_moiLFC);
    if(isempty(kill)==0)
        flipsTmoi( kill ) = [];    %now, the first cross over should be the intended moistLFC
    end
    %add top of profile as last flip entry (for complete cape(z) and because sometimes profile cuts too short to capture EL)
    flipsTmoi = vertcat(flipsTmoi,length(Tv));

else   %if no positive Bmoist

    moiLFCh = NaN;
    flipsTmoi = NaN;
    moiCAPE = NaN;
    CIN = NaN;
    moiELh = NaN;
    CAPEacbl = NaN;

end




%if there is at least 1 Tenv-Tparcel cross over (minimally, the LFC)
if( isempty(moiLFCh)==0 & isnan(moiLFCh)==0 & length(isnan(flipsTmoi)==0)>0 )  

    %go thru all layers bounded by Tenv-Tparcel cross overs
    MINICAPES = [];
    flipsTabove = flipsTmoi;
    Bb = Bmoist;
    for i = 1:length(flipsTabove)-1   
        minicape = [];
        %for k = flipsTabove(i):1:flipsTabove(i+1)  
        for k = flipsTabove(i):1:flipsTabove(i+1)  
            clear Bchunk minisum
            Bchunk = 0.5*(Bb(k)+Bb(k-1))*(h(k)-h(k-1));
            minicape = vertcat(minicape,Bchunk);
            minisum = sum(minicape); %J/kg
        end
        MINICAPES = vertcat(MINICAPES,minisum);   %this is the CAPE as f(k) in each layer bounded by flipsTmoi(k)-flipsTmoi(k+1)
    end  %i




    %find the series of layers bound by the various flipTmoi k indices that maximize cape
    SUMMINICAPES = zeros(length(MINICAPES),1);
    for k = 1:length(MINICAPES) 
        SUMMINICAPES(k) = sum(MINICAPES(1:k));
    end
    kmaxcapelayer = find( SUMMINICAPES == max(SUMMINICAPES) )  ;  %i'm not yet sure if this is a guaranteed CAPE>0... i *think* so?


    if( length(kmaxcapelayer) > 1 )
            kmaxcapelayer = kmaxcapelayer(end);
    end


    %EL is the height at the top of the maximized cape layer
    k_moiEL  = flipsTmoi(kmaxcapelayer+1);
    moiELh   = h(k_moiEL); 
    % EL_pres  = p(k_moiEL)/100
    % p(flipsTmoi)/100

    %if the top of this layer is the top of the soudning, then flag the profile for being too short to fully capture the EL and CAPE
    if ( k_moiEL == length(Tv) )
        ELtooshort = 0;
    end

moiCAPE = SUMMINICAPES(kmaxcapelayer);

else  %no Tenv - Tparcel crossovers

    moiLFCh = NaN;
    moiELh = NaN;
    moiCAPE = NaN;

end




% %%%% find MOIST LFC/EL, using first/last k of Bmoist > 0:
% 
% %vert incidces in full profile where parcel is saturated
% drylayer = find(SatLog == 0 );
% moistlayer = find(SatLog == 1 );
% 
% % %diagnostics:
% % disp(['Zo',num2str(Zo)])
% % disp('  ')
% % disp(['drylayer'])
% % disp(drylayer)
% % disp(['moistlayer'])
% % disp(moistlayer)
% % disp('satlog')
% % disp(SatLog)
% % disp('   ')
% 
% ELtooshort = 1;  %log for if profile is tall enough to capture EL
% if( isempty(moistlayer)==0 & length(find( Bmoist > 0))> 1 )   %if there is a moist layer (I'm guessing there will always be one) and it has at leas a modicem of positive Buoy
% 
%     moistbot = moistlayer(1);  %basically just the LCL
% 
%     %look for Bmoist > 0 layer
%     PosBmoist = find( Bmoist > 0) ;
%     %set lfc and el to be bottom and top of the Bmoist>0 height list
%     k_moiLFC = PosBmoist(1);  %final answer for LFC
%     %k_moiEL = PosBmoist(end);  %first guess of EL
% 
% %     % if the top of the PosBmoist is the top of the profile, then flag the profile for not being deep enough to fully capture the EL
% %     if( k_moiEL == length(moistlayer) )
% %         ELtooshort = 0; %sounding not tall enough to capture EL
% %     end
% 
% %     moiLFCh = h(k_moiLFC);    %agl
% %     moiELh  = h(k_moiEL);     %agl
% 
% else
% 
%     moiCAPE = NaN;
%     dryCAPE = NaN;
%     CIN = NaN;
%     moiLFCh = NaN;
%     moiELh = NaN;
%     CAPEacbl = NaN;
% 
% end









% when LFC last first becomes B>0
% ELtooshort = 1;  %tall enought to capture EL
% if( isempty(moistlayer)==0 )   %if there is a moist lyyer (I'm guessing there will always be one)
%
%     moistbot = moistlayer(1);  %basically just the LCL
%
%     %look for Bmoist > 0 layer with max Bmoist, call the top of that layer the moiEL and bottom moistLFC
%     k_maxBm = find( max(Bmoist)==Bmoist & Bmoist > 0) ;
%
%     if(isempty(k_maxBm)==0) %if there is B>0 layer and a max ID'ed
%         k_maxBm = k_maxBm(1);
%
%         %find EL
%         elflag = 0;
%         for k = k_maxBm : length(Bmoist)
%             if(Bmoist(k) <= 0.0  & elflag == 0)
%                 k_EL = k;
%                 elflag = 1;
%             end
%             %if you get to the top of the moist layer and its still
%             %buoyant, just tag the EL = top  and set a TEwarning
%             if( k == length(Bmoist)  &  Bmoist(k) > 0.0  &  elflag == 0 )
%                 k_EL = k;
%                 ELtooshort = 0;  %too short to capture EL
%                 elflag = 1;
%             end
%         end
%
%         %find LFC
%         lfcflag = 0;
%         for k = k_maxBm : -1 : moistbot
%             if(Bmoist(k) <= 0.0  & lfcflag == 0)
%                 k_moiLFC = k;  %set LFC = k that crosses env profile while moist
%                 lfcflag = 1;
%             end
%
%             %if you get to the bottom of the moist layer and its still
%             %buoyant, just tag the moistLFC = LCL  [e.g., abs unstable
%             %parcel with near grnd superadiab layer]
%             if( k == moistbot  &  Bmoist(k) > 0.0  &  lfcflag == 0 )
%                 k_moiLFC = k;
%                 lfcflag = 1;
%             end
%
%         end
%
%         moiLFCh = h(k_moiLFC);    %agl
%         moiELh  = h(k_EL);     %agl
%
%     else  % no Bmoist > 0 layer
%
%         moiLFCh = NaN;
%         moiELh = NaN;
%         %     moiLFCp = NaN;
%         %     moiELp = NaN;
%
%     end












% %%%% find DRY "LFC", if you really care about it
% drybuoy = find(Bdry(1:drytop)>0);  % find k in dry layer with B>0
% if( isempty(drybuoy)==0  )  %if there is a dry B>0 region
%     k_dryLFC = drybuoy(1);
%     dryLFCh = h(k_dryLFC)    %agl
% end


if( isempty(moiLFCh)==0 & isnan(moiLFCh)==0 )  %if there's a moist LFC

    % find k of acbl top; and LCL+2km layer:
    kacbltop =  find( abs(h - (moiLFCh+1500)) == min(abs(h - (moiLFCh+1500))  ) ) ;    kacbltop = kacbltop(1);    %newvars
    %strip%       klclplus2km =  find( abs(h - (LCLh+2000)) == min(abs(h - (LCLh+2000))  ) ) ;   klclplus2km = klclplus2km(1);   %newvars

    if( k_moiLFC - Zo > 0 )  %if the parcel lift starts at least a point below the moiLFC, calc CIN
        cin = [];
        for k = Zo+1:1:k_moiLFC
            clear Bchunk
            Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1)) ;
            cin = vertcat(cin,Bchunk);
        end
        %strip%        CIN = sum(cin);  %J/kg  - this is just CIN below LFC

        negs = find( cin < 0 );   %the data points in the profile with negative B below the LFC
        CIN = sum( cin(negs) ); %total negative area ONLY below the LFC (i.e., "CIN")

    else  %parcel starts lifting @ about the moiLFC, then no CIN

        CIN = 0.0;

    end


%     % moiCAPE = integral of B above the moist LFC... now done above
%     cape = [];
%     for k = k_moiLFC+1:1:k_moiEL
%         clear Bchunk
%         Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
%         cape = vertcat(cape,Bchunk);
%     end
%     %poss = find( cape > 0 );   %the data points in the profile with + B above the LFC
%     %moiCAPE = sum(cape(poss)); %J/kg ?   - this is CAPE only between LFC and EL
%     moiCAPE = sum(cape); %J/kg


    %if there is a drylayer, calc dry cape:
    drylayer = find(SatLog == 0 );
    if( isempty(drylayer)==0 & drylayer(end)-Zo > 0 ) %if there's a dry layer at least 2 k-points deep - sometimes there may not be if Zo is saturated
       
        drytop = drylayer(end);  %k = one pt below LCL
        % dryCAPE = integral of positive area only in the lift parce's sub-saturated layer
        cape = [];
        for k = Zo+1:1:drytop
            clear Bchunk
            Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
            cape = vertcat(cape,Bchunk);
        end
        poss = find( cape > 0 );   %the data points in the profile with + B in the dry layer
        dryCAPE = sum(cape(poss)); %J/kg ?   - this is + area only in dry layer

    else  % if no dry layer (i.e., parcel is sat at Zo)
        dryCAPE = 0.0;
    end

    %strip%
    %         % CAPElcl = integral of positive and negative area between LCL and LCL+2km
    %         % Lock and Houston (2014) variable
    %         capelcl = [];
    %         for k = LCLi+1:1:klclplus2km
    %             clear Bchunk
    %             Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
    %             capelcl = vertcat(capelcl,Bchunk);
    %         end
    %         CAPElcl = sum(capelcl); %J/kg
    %strip%

    % CAPEacbl = integral of positive and negative area between moist lfc and lfc+1.5km
    % inspired by Lock and Houston (2014)
    capeacbl = [];
    for k = k_moiLFC+1:1:kacbltop
        clear Bchunk
        Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
        capeacbl = vertcat(capeacbl,Bchunk);
    end
    %poss = find( capeacbl > 0 );   %the data points in the profile with + B in the dry layer
    %CAPEacbl = sum(capeacbl(poss)); %J/kg ?   - this is + area only in dry layer
    CAPEacbl = sum(capeacbl); %J/kg (+ and -)

    %%J%%        kLFCplus1 = find( abs(h - (LFCh+1000)) == min(abs(h - (LFCh+1000))  ) ) ;    kLFCplus1 = kLFCplus1(1);    %newvars
    %%J%%        kLFCplus2 = find( abs(h - (LFCh+2000)) == min(abs(h - (LFCh+2000))  ) ) ;    kLFCplus2 = kLFCplus2(1);    %newvars
    %%J%%        kLFCplus3 = find( abs(h - (LFCh+3000)) == min(abs(h - (LFCh+3000))  ) ) ;    kLFCplus3 = kLFCplus3(1);    %newvars
    %%J%%        kLFCplus4 = find( abs(h - (LFCh+4000)) == min(abs(h - (LFCh+4000))  ) ) ;    kLFCplus4 = kLFCplus4(1);    %newvars
    %%J%%
    %%J%%        if(isempty(kLFCplus1) == 0 & kLFCplus1 < ELi)
    %%J%%            % CAPE = integral of positive and negative area between LFC and
    %%J%%            % 1 km above it (if that +1km ht is below the EL)
    %%J%%            cape = [];
    %%J%%            for k = LFCi+1:1:kLFCplus1
    %%J%%                clear Bchunk
    %%J%%                Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
    %%J%%                cape = vertcat(cape,Bchunk);
    %%J%%            end
    %%J%%            CAPE_L01 = sum(cape); %J/kg ?   - this is CAPE + CIN between LFC and EL
    %%J%%        else
    %%J%%            CAPE_L01 = NaN;
    %%J%%        end
    %%J%%
    %%J%%        if( isempty(kLFCplus1) == 0 & kLFCplus1 < ELi & isempty(kLFCplus2) == 0 & kLFCplus2 < ELi )
    %%J%%            % CAPE = integral of positive and negative area between LFC + 1km  and
    %%J%%            % LFC + 2 km above it (if that +1km ht is below the EL)
    %%J%%            cape = [];
    %%J%%            for k = kLFCplus1:1:kLFCplus2
    %%J%%                clear Bchunk
    %%J%%                Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
    %%J%%                cape = vertcat(cape,Bchunk);
    %%J%%            end
    %%J%%            CAPE_L12 = sum(cape); %J/kg ?   - this is CAPE + CIN between LFC and EL
    %%J%%        else
    %%J%%            CAPE_L12 = NaN;
    %%J%%        end
    %%J%%
    %%J%%
    %%J%%        if( isempty(kLFCplus2) == 0 & kLFCplus2 < ELi & isempty(kLFCplus3) == 0 & kLFCplus3 < ELi )
    %%J%%            % CAPE = integral of positive and negative area between LFC + 2km  and
    %%J%%            % LFC + 3 km above it (if that +1km ht is below the EL)
    %%J%%            cape = [];
    %%J%%            for k = kLFCplus2:1:kLFCplus3
    %%J%%                clear Bchunk
    %%J%%                Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
    %%J%%                cape = vertcat(cape,Bchunk);
    %%J%%            end
    %%J%%            CAPE_L23 = sum(cape); %J/kg ?   - this is CAPE + CIN between LFC and EL
    %%J%%        else
    %%J%%            CAPE_L23 = NaN;
    %%J%%        end
    %%J%%
    %%J%%        if( isempty(kLFCplus3) == 0 & kLFCplus3 < ELi & isempty(kLFCplus4) == 0 & kLFCplus4 < ELi )
    %%J%%            % CAPE = integral of positive and negative area between LFC + 3km  and
    %%J%%            % LFC + 4 km above it (if that +1km ht is below the EL)
    %%J%%            cape = [];
    %%J%%            for k = kLFCplus3:1:kLFCplus4
    %%J%%                clear Bchunk
    %%J%%                Bchunk = 0.5*(B(k)+B(k-1))*(h(k)-h(k-1));
    %%J%%                cape = vertcat(cape,Bchunk);
    %%J%%            end
    %%J%%            CAPE_L34 = sum(cape); %J/kg ?   - this is CAPE + CIN between LFC and EL
    %%J%%        else
    %%J%%            CAPE_L34 = NaN;
    %%J%%        end


    % I dont fully recall what glitches these statements below were designed to catch:

    %strip%
    %     if( abs(ELp - p(Zo)) < 2000 )
    %         CAPE = -9999;
    %         CIN  = -9999;
    % %strip%        LFCp = -9999;
    %         LFCh = -9999;
    % %strip%        ELp = -9999;
    %     end
    %strip%


    %strip%
    %     %  kill parcel T profile above EL (if it exists)
    %     %  if EL doesnt exist, set "EL" = 150 mb and flag CIN = -999, CAPE = -999,
    %     %LFC = -999
    %     if( abs(ELp - p(Zo)) > 2000 )
    %         %parcelTprof(ELi+1:1:length(T)) = NaN;
    %         parcelTVprof(ELi+1:1:length(T)) = NaN;
    %         %parcelTprof(1:1:Zo-1) = NaN;
    %         parcelTVprof(1:1:Zo-1) = NaN;
    %     end
    %strip%

    %strip%     parcelTvprof(ELi+1:end) = NaN;

else   %not a moistLFC

    moiCAPE = NaN;
    dryCAPE = NaN;
    CIN = NaN;
    moiLFCh = NaN;
    moiELh = NaN;
    CAPEacbl = NaN;

end %isempty(moiLFCh)==0 check


