function [sfcCAPE_moi, sfc_CIN, sfcLCL, sfcLFC, sfcEL, sfc_CAPEACBL, sfcCAPEpse_moi, sfc_CINpse, sfcLFCpse, sfcELpse, sfc_CAPEACBLpse,...
          muCAPE_moi,  mu_CIN,  muLCL,  muLFC,  muEL,  mu_CAPEACBL,  muCAPEpse_moi,  mu_CINpse,  muLFCpse,  muELpse,  mu_CAPEACBLpse,...
          mlCAPE_moi,  ml_CIN,  mlLCL,  mlLFC,  mlEL,  ml_CAPEACBL,  mlCAPEpse_moi,  ml_CINpse,  mlLFCpse,  mlELpse,  ml_CAPEACBLpse, ...
          hMU,sfcELtooshort,muELtooshort,mlELtooshort,sfcELtooshortpse,muELtooshortpse,mlELtooshortpse] = SOUNDING_SFC_MU_ML_ZAwrf2(pres, tdry, TDc, h);

%function [sfcCAPE, sfc_negCINsum, sfc_negCINfract, sfcLCL, sfcLFC, sfcEL, sfc_CAPEACBL, sfcCAPEpse, sfc_negCINsumpse, muCAPE,  mu_negCINsum,  mu_negCINfract,  muLCL, muLFC,  muEL,  mu_CAPEACBL,  muCAPEpse,  mu_negCINsumpse, mlCAPE,  ml_negCINsum,  ml_negCINfract, mlLCL,  mlLFC,  mlEL,  ml_CAPEACBL,  mlCAPEpse,  ml_negCINsumpse, hMU] = SOUNDING_SFC_MU_ML_jdebug(pres, tdry, TDc, h);
      
      


%some history:

%g: same as "_ftest", but trims off sounding above 28 km, corrected the flipsT
%to look only from Zo and up.

%h: same as g, but added in a factor to control the mysterious
%discontinuity of the parcel temperature curve where it saturates

%j: calls John Peter's parcel lifting algorithm

%k: uses calcCAPE_MoistAndDry.m to deal with superadiabitic layer nonsense
%so it eturns only moist CAPE and LFC = moistLFC. I calculate dry CAPE, but
%am not currently passing it through (trivial to do if you want to do it
%later)


% 'sfc' = calcualted for surface parcel
% 'mu'  = calculated for most unstable parcel (the most cape)
% 'ml'  = for parcel representative of the mean of atmosphere between sfc
%           and 100 mb above it

% hMU = altitude (m AGL) of the mu parcel origin (Zo)
% LCL (m AGL)
% LFC = height (m AGL) of the parcel crossing the env Tv profile at the bottom of the layer
%         which contains the largest positive buoyancy
% EL the height (m AGL) of the parcel crossing the env Tv profile at the top of the layer
%         which contains the largest positive buoyancy
% CAPE = sum of parcel buoyancy (+ and - area) between the LFC and EL -
%       using parcel and environment Tv
% CIN = sum of parcel buoyancy (+ and - area) between Zo
%       and EL - using parcel and environment Tv
% CAPEacbl = CAPE measured between LFC and 1.5km above it
% CAPElcl = CAPE measured between LCL and 2 km above it
% negCINfract = fraction of k levels between Zo and LFC which contain
%       negative (area) buoyancy
% negCINsum = 'CIN' that only sums the negative area between Zo - LFC
% EILXX = bottom and top of effective inflow layer [height_of_bottom
%       height_of_top] = the layer for which there is > 100 J/kg CAPE and <
%       XX j/kg of CIN (CIN = sum of + and - area)

% metrics with 'p' suffices are the pressure levels (mb) rather than height.

%teSFC,ML,MU = 1 log of soundings that are/arent deep enough to track full
                    %parcel. 


%%


% %{  diagnostics
% 
% 
% clear 
% %load Cell16Tester.mat
% load Cell3Tester.mat
% 
% %enoch's test indices
% x = 7 +1; y = 9 +1; t = 120 +1;
% 
% height = wrf_hasl(x,y,:,t) - wrf_hasl(x,y,1,t);         %m AGL
% pres = wrf_pres(x,y,:,t);           %mb
% tdry = wrf_T(x,y,:,t) - 273.15;      %C
% mr = wrf_mr(x,y,:,t) ;              %kg/kg
% U = wrf_u(x,y,:,t);                 %m/s
% V = wrf_v(x,y,:,t);                 %m/s
% rh = wrf_rh(x,y,:,t);               %  0-100%
% 
% %fix index order via permutation:
% U = permute(U,[3 2 1]);
% V = permute(V,[3 2 1]);
% ht = permute(height,[3 2 1]);
% pres = permute(pres,[3 2 1]);
% tdry = permute(tdry,[3 2 1]);
% mr = permute(mr,[3 2 1]);
% rh = permute(rh,[3 2 1]);
% 
% [sHt sTi] = size(tdry);
% 
% clear wrf_hasl wrf_mr wrf_pres wrf_rh wrf_T wrf_times wrf_u wrf_v height
% 
% 
% %%%%% chop out stray mr < 0 quirk, replace with interpolation:
% deadpixels = find (  mr <= 0 );
% mr(deadpixels) = NaN;
% for nm = 1:length(deadpixels)
%     mr(deadpixels(nm)) = mean ( mean(  mr(deadpixels(nm)-1:deadpixels(nm)+1)  ,'omitnan') ) ;
% end
% %%%%%
% 
% %calc some other things:
% p = pres*100;
% E = 0.622;
% T = tdry + 273.15;       % temperature K
% f = mr / 0.622;
% e = (p.*f)./(1+f);   %pa
% sphu  = mr./(1+mr);  %specific humidity
% 
% % calculate potential temp:
% theta = T.*(100000./p).^0.286; %K
% 
% %emanual dew point:
% eee = e/100; %vap press [mb]
% TdC =273.15+  243.5./( (17.67./log(eee./6.112)) -1) - 273.15;    %  deg C      K.Emanuel
% 
% 
% %%% define some constants:
% cpd = 1006; % J/kg K
% g = 9.8;    % m/s2
% cpv = 1870; %J/kg K
% E = 0.622;
% lv0 = 2501000;
% Rv = 461.5;
% Rd = 287.04;
% cw = 4190;             % heat capacity of water
% cc = 2320;
% cl = 4200;
% cvv = 1410;
% cvd = 719;
% 
% alv = lv0 - cc*(T-273.15); %J/kg
% 
% Tv    = T.*(1+0.608.*mr);  %K
% thetav = Tv.*(100000./p).^0.286; %K
% 
% teA = T.*(100000./p).^(Rd./(cpd + cl.*mr));
% teB = exp( (lv0.*mr)./((cpd + mr*cl).*T)   );
% thetae = teA .* teB;
% 
% 
% 
% clear x y 
% 
% 
% % % % % % % % % % % debuggers
% %  h = ht; TDc = TdC(:,t);  pres = pres(:,t);  tdry = tdry(:,t);
% %  h = ht; TDc = TdC;  pres = pres;  tdry = tdry;
% %  h_preserve = h;    TDc_preserve = TDc;   pres_preserve = pres;  tdry_preserve = tdry; 
% % 
% % 
% % 
% % % 
% % hhh = HH(t);
% % mmm = MM(t);
% % sss = SS(t);
% % TTT = t;
% % % % % % % % % %
% 
% %}






sfcELtooshort = 2;  %tall enought to capture EL - set to 3 to start. It will most likely turn to 0 or 1 in calcCAPE
muELtooshort = 2;  %tall enought to capture EL - set to 3 to start. It will most likely turn to 0 or 1 in calcCAPE
mlELtooshort = 2;  %tall enought to capture EL - set to 3 to start. It will most likely turn to 0 or 1 in calcCAPE

sfcELtooshortpse = 2;  %tall enought to capture EL - set to 3 to start. It will most likely turn to 0 or 1 in calcCAPE
muELtooshortpse = 2;  %tall enought to capture EL - set to 3 to start. It will most likely turn to 0 or 1 in calcCAPE
mlELtooshortpse = 2;  %tall enought to capture EL - set to 3 to start. It will most likely turn to 0 or 1 in calcCAP




%trim off sounding above 20 km:
kill = find(h > 20000);
h(kill) = [];
TDc(kill) = [];
tdry(kill) = [];
pres(kill) = [];



% hhh = sprintf( '%02d', hhh ) ;
% mmm = sprintf( '%02d', mmm ) ;
% sss = sprintf( '%02d', sss ) ;


[aa bb] = size(pres);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%chop out NaNs, they mess up things:

deadpixels = find (  isnan(tdry) == 1 );
h(deadpixels) = [];
tdry(deadpixels) = [];
pres(deadpixels) = [];
TDc(deadpixels) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TERRAIN_HT = Alt(1);    % meters - altitude of the ground above sea level at sonde launch site
% h = Alt-TERRAIN_HT;   % m agl


%look for repeated h's:
repeats = [];
for i = 1:length(h)-1
    if(h(i+1) == h(i))
        repeats = vertcat(repeats,i);
    end
end
if( isempty(repeats) == 0 )
    h(repeats) = [];
    tdry(repeats) = [];
    %RH(repeats) = [];
    TDc(repeats) = [];
    pres(repeats) = [];
end


% if Tdew > Tdry, "fix" it:
soaked = find( TDc >= tdry );
TDc(soaked) = tdry(soaked) - 0.1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    if the sounding isn't deep, just punt on this one
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [sfcCAPE_moi, sfc_CIN, sfcLCL, sfcLFC, sfcEL, sfc_CAPEACBL, sfcCAPEpse_moi, ...
%           muCAPE_moi, mu_CIN, muLCL, muLFC, muEL, mu_CAPEACBL, muCAPEpse_moi, ...
%           mlCAPE_moi, ml_CIN, mlLCL, mlLFC, mlEL, ml_CAPEACBL, mlCAPEpse_moi, ...
%           hMU] 



if(h(end) < 8000)
    
%     sfcCAPE = NaN; sfcLCL = NaN; sfcLFC = NaN; sfcEL = NaN; 
%     sfc_negCINfract = NaN; sfc_negCINsum = NaN; sfc_CAPEACBL = NaN; 
%     sfcCAPEpse = NaN; sfc_negCINsumpse = NaN;
%     %strip% sfcLCLp = NaN; sfcLFCp = NaN; sfcELp = NaN; sfc_CAPELCL = NaN; sfcCIN = NaN;
    
    sfcCAPE_moi=NaN; sfcCAPE_dry=NaN; sfc_CIN=NaN; sfcLFC=NaN; sfcEL=NaN; sfc_CAPEACBL=NaN;
    sfcCAPEpse_moi=NaN; sfcCAPEpse_dry=NaN; sfcLCL=NaN;

    sfc_CINpse=NaN; sfcLFCpse=NaN; sfcELpse=NaN; sfc_CAPEACBLpse=NaN;
 

%     muCAPE = NaN; muLCL = NaN; muLFC = NaN; muEL = NaN; 
%     mu_negCINfract = NaN; mu_negCINsum = NaN; mu_CAPEACBL = NaN; 
%     muCAPEpse = NaN; mu_negCINsumpse = NaN;
%     %strip% muLCLp = NaN; muLFCp = NaN; muELp = NaN; mu_CAPELCL = NaN; muCIN = NaN; 
    
    muCAPE_moi=NaN; muCAPE_dry=NaN; mu_CIN=NaN; muLFC=NaN; muEL=NaN; mu_CAPEACBL=NaN;
    muCAPEpse_moi=NaN; muCAPEpse_dry=NaN; muLCL=NaN;
    kMU = NaN; hMU  = -98989;
    mu_CINpse=NaN; muLFCpse=NaN; muELpse=NaN; mu_CAPEACBLpse=NaN;

%     mlCAPE = NaN; mlLCL = NaN; mlLFC = NaN; mlEL = NaN; 
%     ml_negCINfract = NaN; ml_negCINsum = NaN; ml_CAPEACBL = NaN; 
%     mlCAPEpse = NaN; ml_negCINsumpse = NaN;
%    %mlLCLp = NaN; mlLFCp = NaN; mlELp = NaN; ml_CAPELCL = NaN; mlCIN = NaN; 

    mlCAPE_moi=NaN; mlCAPE_dry=NaN; ml_CIN=NaN; mlLFC=NaN; mlEL=NaN; ml_CAPEACBL=NaN;
    mlCAPEpse_moi=NaN; mlCAPEpse_dry=NaN; mlLCL=NaN;
    ml_CINpse=NaN; mlLFCpse=NaN; mlELpse=NaN; ml_CAPEACBLpse=NaN;    
    
%strip% EIL250 = [NaN NaN]; EIL100 = [NaN NaN]; 
%strip% EIL50 = [NaN NaN]; EIL25 = [NaN NaN]; EIL10 = [NaN NaN];
    
%     teSFC = 0;      
%     teMU = 0;      
%     teML = 0; 
    
    % mu_CAPEL01 = NaN; mu_CAPEL12 = NaN; mu_CAPEL23 = NaN; mu_CAPEL34 = NaN;
    
%     disp('           ')
%     disp('short sounding (< 12 km)... punt!')
%     disp('           ')
%     
    return
    
end









% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % interp base state to a Zii-m vert resolution
% % to pinpoint desired ht.
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zii = 20;
h2 = [h(1):Zii:h(end)]';

TDc2 = interp1(h,TDc,h2);
pres2 = interp1(h,pres,h2);
tdry2 = interp1(h,tdry,h2);

% hold on
% plot(h,'bx')

h = h2;
TDc = TDc2;
pres = pres2;
tdry = tdry2;
clear h2 TDc2 pres2 tdry2
% %%%%%%%%% done with bonus interpolation






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% calculate sounding metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% some derived fields:
p = pres*100;         % convert pressure (pa)
T = tdry + 273.15;       % temperature K
%theta = T.*(1000./pres).^0.286;  %potential temperature K
%e = 6.112*exp(17.67.*TDc/243.5);  % vapor pressure (mb)
%mr = (287.04/461.5)*(e./(pres-e));  %mixing ratio (kg/kg)
%Tv = T.*(1+0.608*mr);           % virtual Temp K


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


e =   6.112*exp( (17.67*TDc )./(243.5+TDc ) )*100; %pa
esat = 6.112*exp( (17.67*tdry)./(243.5+tdry) )*100; % pa
mrsat = E*esat./(p-esat); %kg/kg
RH = e./esat;
mr = 0.622*(e./(p-e));  %kg/kg
alv = lv0 - cc*(T-273.15); %J/kg

tdry  = T - 273.15;
Tv    = T.*(1+0.608.*mr);
theta = T.*(100000./p).^0.286;
thetav = Tv.*(100000./p).^0.286;

teA = T.*(100000./p).^(Rd./(cpd + cl.*mr));
teB = exp( (lv0.*mr)./((cpd + mr*cl).*T));
teC = RH.^(  (-mr*Rv) ./ (cpd + cl.*mr)   );
thetae = teA .* teB .* teC;



% Tv_preserve = Tv;
% Tt_preserve = T;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ok, now start lifing parcels:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h4km = 4000 ;
k4km  =  find( abs(h - h4km) == min(abs(h-h4km)  ) ) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   sfc parcel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Zo = 1;


clear T0 mr0
mr0 = mr(Zo);
T0 = T(Zo);

% % I think this is moot with John's code:
% %entropy "slope"
% sl =    ( cpd  +  mr0*cw+  (alv.*alv.*mrsat)./ (Rv*T.*T) )./T;  %for reversible
% % slp = ( cpd + mrsat*cw + (alv.*alv.*mrsat)./(Rv*T.*T) ) ./T ; % check


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% find LCL for height level = Zo   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%from Emanuel
TSAT = 55 + 2840/( 3.5*log(T0)-log(e(Zo)/100)-4.805 );
LCL = h(Zo) + (cpd/g)*(T0-TSAT)*(1+mr0*(cpv/cpd))/(1+mr0); %LCL height (mASL???  - actually, I think AGL?)
ch = T0/(1669 - 122*RH(Zo) - T0);

% drylayer = find(h < LCL);
% sfcLCLi = drylayer(end);
% sfcLCL = h(sfcLCLi)
% sfcLCLp = p(sfcLCLi)


sphum  = mr./(1+mr);  %specific humidity [mass vapor/mass(vapor+dry)]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define the Tv profile of the lifted parcel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%reversible 'warm' moist lifting:
[revT_rho_lif_sfc,IsSat_JNM] = lift_parcel_adiabatic_jpeters_ZAwrf2(T, pres*100, sphum , Zo, 0, 0, h, 3, 2);
%find LCL from johns lifting saturation
sfcLCLi = find(IsSat_JNM==1);  sfcLCLi = sfcLCLi(1);
sfcLCL = h(sfcLCLi);
%  sfcLCLp = p(sfcLCLi)/100

[sfcCAPE_moi, sfcCAPE_dry, sfc_CIN, sfcLFC, sfcEL, sfc_CAPEACBL, sfcELtooshort] = calcCAPE_MoistAndDry_ZAwrf2(Tv,revT_rho_lif_sfc,h,Zo,sfcLCLi,IsSat_JNM);   %chokes here
%note, this LFC is the reversible "moist LFC", which you will use henceforth

%daignostic print:
%  lfci = find(sfcLFC >= h) ; lfci = lfci(end); sfcLFCp = p(lfci)/100
%  eli = find(sfcEL >= h) ; eli = eli(end); sfcELp = p(eli)/100 

% moist,warm psuedoadiabatic, just for the cape for now:
[pseT_rho_lif_sfc,IsSat_JNM] = lift_parcel_adiabatic_jpeters_ZAwrf2(T, pres*100, sphum, Zo, 0, 1 /  max( gradient(h) ), h, 3, 2);   %(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
[sfcCAPEpse_moi, sfcCAPEpse_dry, sfc_CINpse, sfcLFCpse, sfcELpse, sfc_CAPEACBLpse, sfcELtooshortpse] = calcCAPE_MoistAndDry_ZAwrf2(Tv,pseT_rho_lif_sfc,h,Zo,sfcLCLi,IsSat_JNM);
clear blah2 blah3 blah4 blah5 blah6







%strip% dist = 400.0;
%%J%% %define ThompsonEtAl(2007)-esque  EIL(h_bot,h_top) - where cape > 100 and cin > -250
%%J%% %requires two heights to meet cape/cin criteria
%%J%% %requires the mean distance between heights to be < dist [meters]
%%J%% dist = 400.0;
%%J%%
%%J%% keilay250 = find(CAPEtv > 100.0 & CINtv > -250.0);
%%J%% if(isempty(keilay250) == 0 & length(keilay250)>= 2 )
%%J%%     heilay250 = h(keilay250); dheilay = heilay250; dheilay(end) = [];
%%J%%     for ke = 1:length(dheilay)
%%J%%         dheilay(ke) =  heilay250(ke+1) - heilay250(ke);
%%J%%     end
%%J%%     if(mean(dheilay) < dist)
%%J%%         EIL250 = [heilay250(1), heilay250(end)];
%%J%%         EIL250p = [ pres(keilay250(1)), pres(keilay250(end))];
%%J%%     else
%%J%%         EIL250 = [-9999 -9999];   %a flag to see if there are depths that meet cape/cin crit, but not the dz crit.
%%J%%     end
%%J%% else
%%J%%     EIL250 = [NaN NaN];
%%J%% end


%%J%% %define ThompsonEtAl(2007)-esque  EIL(h_bot,h_top) - where cape > 100 and cin > -100
%%J%% %requires two heights to meet cape/cin criteria
%%J%% %requires the mean distance between heights to be < 500m
%%J%% keilay100 = find(CAPEtv > 100.0 & CINtv > -100.0);
%%J%% if(isempty(keilay100) == 0 & length(keilay100)>= 2 )
%%J%%     heilay100 = h(keilay100); dheilay = heilay100; dheilay(end) = [];
%%J%%     for ke = 1:length(dheilay)
%%J%%         dheilay(ke) =  heilay100(ke+1) - heilay100(ke);
%%J%%     end
%%J%%     if(mean(dheilay) < dist)
%%J%%         EIL100 = [heilay100(1), heilay100(end)];
%%J%%         EIL100p = [ pres(keilay100(1)), pres(keilay100(end))];
%%J%%     else
%%J%%         EIL100 = [-9999 -9999];   %a flag to see if there are depths that meet cape/cin crit, but not the dz crit.
%%J%%     end
%%J%% else
%%J%%     EIL100 = [NaN NaN];
%%J%% end


% %strip%
% %define ThompsonEtAl(2007)-esque  EIL(h_bot,h_top) - where cape > 100 and cin > -50
% %requires two heights to meet cape/cin criteria
% %requires the mean distance between heights to be < 500m
% keilay50 = find(CAPEtv > 100.0 & CINtv > -50.0);
% if(isempty(keilay50) == 0 & length(keilay50)>= 2 )
%     heilay50 = h(keilay50); dheilay = heilay50; dheilay(end) = [];
%     for ke = 1:length(dheilay)
%         dheilay(ke) =  heilay50(ke+1) - heilay50(ke);
%     end
%     if(mean(dheilay) < dist)
%         EIL50 = [heilay50(1), heilay50(end)];
%         EIL50p = [ pres(keilay50(1)), pres(keilay50(end))];
%     else
%         EIL50 = [-9999 -9999];   %a flag to see if there are depths that meet cape/cin crit, but not the dz crit.
%     end
% else
%     EIL50 = [NaN NaN];
% end
% %strip%


% %strip%
% %define ThompsonEtAl(2007)-esque  EIL(h_bot,h_top) - where cape > 100 and cin > -25
% %requires two heights to meet cape/cin criteria
% %requires the mean distance between heights to be < 500m
% keilay25 = find(CAPEtv > 100.0 & CINtv > -25.0);
% if(isempty(keilay25) == 0 & length(keilay25)>= 2 )
%     heilay25 = h(keilay25); dheilay = heilay25; dheilay(end) = [];
%     for ke = 1:length(dheilay)
%         dheilay(ke) =  heilay25(ke+1) - heilay25(ke);
%     end
%     if(mean(dheilay) < dist)
%         EIL25 = [heilay25(1), heilay25(end)];
%         EIL25p = [ pres(keilay25(1)), pres(keilay25(end))];
%     else
%         EIL25 = [-9999 -9999];   %a flag to see if there are depths that meet cape/cin crit, but not the dz crit.
%     end
% else
%     EIL25 = [NaN NaN];
% end
% %strip%


% %strip%
% %define ThompsonEtAl(2007)-esque  EIL(h_bot,h_top) - where cape > 100 and
% %cin = low
% %requires two heights to meet cape/cin criteria
% %requires the mean distance between heights to be < 500m
% keilay10 = find(CAPEtv > 100.0 & CINtv > -10.0);
% if(isempty(keilay10) == 0 & length(keilay10)>= 2 )
%     heilay10 = h(keilay10); dheilay = heilay10; dheilay(end) = [];
%     for ke = 1:length(dheilay)
%         dheilay(ke) =  heilay10(ke+1) - heilay10(ke);
%     end
%     if(mean(dheilay) < dist)
%         EIL10 = [heilay10(1), heilay10(end)];
%         EIL10p = [ pres(keilay10(1)), pres(keilay10(end))];
%     else
%         EIL10 = [-9999 -9999];   %a flag to see if there are depths that meet cape/cin crit, but not the dz crit.
%     end
% else
%     EIL10 = [NaN NaN];
% end
% %strip%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% most unstable parcel lifting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


kMU = find( thetae(1:k4km) == max(thetae(1:k4km))  );
kMU = kMU(1);


%kMU = kMU(1);  %I put this here because if I set cape = 0 if lfc doesn't
%exist, there are lots of 0.0

if(isempty(kMU)) % if no mu CAPE
    kMU = NaN;
    hMU = NaN;
    muCAPE = NaN;
    muLCL = NaN;
    muLFC = NaN;
    muEL = NaN;

    mu_negCINfract = NaN;
    mu_negCINsum = NaN;
    mu_CAPEACBL = NaN;
    
%strip%    mu_CAPELCL = NaN;
%strip%    muLCLp = NaN;
%strip%    muLFCp = NaN;
%strip%    muELp = NaN;
%strip%    muCIN = NaN;    
    
    %%J%%     mu_CAPEL01 =  NaN;     %newvars
    %%J%%     mu_CAPEL12 =  NaN;       %newvars
    %%J%%     mu_CAPEL23 =  NaN;       %newvars
    %%J%%     mu_CAPEL34 =  NaN;
    
%strip%    parcelTprof_mu =  pres; parcelTprof_mu(:) = NaN;
%strip%    parcelTVprof_mu =  pres; parcelTVprof_mu(:) = NaN;
    
%strip%        teMU = NaN;
    
else % if mu CAPE
    
    hMU = h(kMU) ;
    
    
    Zo = kMU;
    
    
    clear T0 mr0
    mr0 = mr(Zo);
    T0 = T(Zo);
    
    % % I think this is moot with John's code:
    % %entropy "slope"
    % sl =    ( cpd  +  mr0*cw+  (alv.*alv.*mrsat)./ (Rv*T.*T) )./T;  %for reversible
    % % slp = ( cpd + mrsat*cw + (alv.*alv.*mrsat)./(Rv*T.*T) ) ./T ; % check
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% find LCL for height level = Zo   %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %from Emanuel
    TSAT = 55 + 2840/( 3.5*log(T0)-log(e(Zo)/100)-4.805 );
    LCL = h(Zo) + (cpd/g)*(T0-TSAT)*(1+mr0*(cpv/cpd))/(1+mr0); %LCL height (mASL???  - actually, I think AGL?)
    ch = T0/(1669 - 122*RH(Zo) - T0);

%     drylayer = find(h < LCL);
%     muLCLi = drylayer(end);
%     muLCL = h(muLCLi);
%    muLCLp = p(muLCLi)    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% define the Tv profile of the lifted parcel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    %reversible 'warm' moist lifting:
    [revT_rho_lif_mu,IsSat_JNM] = lift_parcel_adiabatic_jpeters_ZAwrf2(T, pres*100, sphum, Zo, 0, 0, h, 3, 2);
    %find LCL from johns lifting saturation
    muLCLi = find(IsSat_JNM==1);  muLCLi = muLCLi(1);
    muLCL = h(muLCLi);
    %muLCLp = p(muLCLi)/100
    
    [muCAPE_moi, muCAPE_dry, mu_CIN, muLFC, muEL, mu_CAPEACBL, muELtooshort] = calcCAPE_MoistAndDry_ZAwrf2(Tv,revT_rho_lif_mu,h,Zo,muLCLi,IsSat_JNM);  
    %note, this LFC is the reversible "moist LFC", which you will use henceforth

%     lfci = find(muLFC >= h) ; lfci = lfci(end); muLFCp = p(lfci)/100
%     eli = find(muEL >= h) ; eli = eli(end); muELp = p(eli)/100


    % moist,warm psuedoadiabatic, just for the cape for now:
    [pseT_rho_lif_mu,IsSat_JNM] = lift_parcel_adiabatic_jpeters_ZAwrf2(T, pres*100, sphum, Zo, 0, 1/max( gradient(h) ), h, 3, 2);   %(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
    [muCAPEpse_moi, muCAPEpse_dry, mu_CINpse, muLFCpse, muELpse, mu_CAPEACBLpse, muELtooshortpse] = calcCAPE_MoistAndDry_ZAwrf2(Tv,pseT_rho_lif_mu,h,Zo,muLCLi,IsSat_JNM);
    clear blah2 blah3 blah4 blah5 blah
    
    
end    
    
  






    
%strip%        mu_CAPELCL = CAPELCL(kMU);
    
%strip%    muCIN = CINtv(kMU);
%strip%    muLCLp = LCLptv(kMU);
%strip%    muLFCp = LFCptv(kMU);
%strip%    muELp = ELptv(kMU);
    
    
    %%J%%     mu_CAPEL01 =  CAPEL01(kMU);     %newvars
    %%J%%     mu_CAPEL12 =  CAPEL12(kMU);      %newvars
    %%J%%     mu_CAPEL23 =  CAPEL23(kMU);      %newvars
    %%J%%     mu_CAPEL34 =  CAPEL34(kMU);
    
    
%strip%        teMU = tallenough(kMU);
    %kill =  find( abs(pres - LCLptv(kMU))  == min(abs(pres - LCLptv(kMU))  ) ) ;
    
%strip%        parcelTVprof_mu =  parcelTVprof_zos(kMU,:);
%strip%        parcelTprof_mu =  parcelTprof_zos(kMU,:);
    %parcelTVprof_mu(1:kMU) = NaN;
    
    %strip%end
    

    







%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ML parcel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Zo = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if you want to do average of lowest [PDEPTH] mb:

Tv_preserve = Tv;
T_preserve = T;
TDc_preserve = TDc;
mr_preserve = mr;
sphum_preserve = sphum;
thetav_preserve = thetav;
theta_preserve = theta;
e_preserve = e;
RH_preserve = RH;
esat_preserve = esat;

%make a raw profile with Zo modded to the ml value for ml parcel lifting
Tv_4ml = Tv;
T_4ml = T;
% TDc_4ml = TDc;
mr_4ml = mr;
sphum_4ml = sphum; 
% thetav_4ml = thetav;
% theta_4ml = theta;
% e_4ml = e;
% RH_4ml = RH;
% esat_4ml = esat;



PDEPTH = 100;   % lowest # mb layer of sounding
playtop = pres(1) - PDEPTH;
layerdepth = find( pres - playtop >= 0 );
layertop = layerdepth(end);
meanthetav = mean(thetav(1:layertop));
meantheta = mean(theta(1:layertop));
meanmr = mean(mr(1:layertop));

% % % rewrite lower layer of sonde with mean theta and mr:
% thetav(1:layertop) = meanthetav;
% theta(1:layertop) = meantheta;
% mr(1:layertop) = meanmr;

Tv_4ml(Zo) = meanthetav./(100000./p(Zo)).^0.286 ;
T_4ml(Zo)  = meantheta./(100000./p(Zo)).^0.286 ;
mr_4ml(Zo) = meanmr;
sphum_4ml = mr_4ml./(1+mr_4ml);


% theta = T.*(100000./p).^0.286;
% thetav = Tv.*(100000./p).^0.286;

% %recalc these for debugging
e_4ml = p.*(mr_4ml./0.622)./(1 + (mr_4ml./0.622)); %mb
eeee_4ml = e_4ml/100;
alph_4ml = (log(eeee_4ml./6.112))/17.67;
TDc_4ml = alph_4ml.*243.5./(1-alph_4ml);
% 
% esat = 6.112*exp( (17.67*(T-273.15))./(243.5+(T-273.15)) )*100; % pa
% RH = e./esat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% done with PDEPTH mean rewrite of sonde

clear T0 mr0
mr0 = mr_4ml(Zo);
T0 = T_4ml(Zo);

% % I think this is moot with John's code:
% %entropy "slope"
% sl =    ( cpd  +  mr0*cw+  (alv.*alv.*mrsat)./ (Rv*T.*T) )./T;  %for reversible
% % slp = ( cpd + mrsat*cw + (alv.*alv.*mrsat)./(Rv*T.*T) ) ./T ; % check


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% find LCL for height level = Zo   %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%from Emanuel
TSAT = 55 + 2840/( 3.5*log(T0)-log(e(Zo)/100)-4.805 );
LCL = h(Zo) + (cpd/g)*(T0-TSAT)*(1+mr0*(cpv/cpd))/(1+mr0); %LCL height (mASL???  - actually, I think AGL?)
ch = T0/(1669 - 122*RH(Zo) - T0);

%aboveLCL = find(h >= LCL);
%LCLi = aboveLCL(1);
%LCLh = h(LCLi);        % meters (agl) - I think?
%LCLp = p(LCLi);        % pa

%drylayer = find(h < LCL);
% mlLCLi = drylayer(end);
% mlLCL = h(mlLCLi);
% mlLCLp = p(mlLCLi)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define the Tv profile of the lifted parcel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[revT_rho_lif_ml,IsSat_JNM] = lift_parcel_adiabatic_jpeters_ZAwrf2(T_4ml, pres*100, sphum_4ml, Zo, 0, 0, h, 3, 2);
%find LCL from johns lifting saturation
mlLCLi = find(IsSat_JNM==1);  mlLCLi = mlLCLi(1);
mlLCL = h(mlLCLi);
%mlLCLp = p(mlLCLi)  /100

[mlCAPE_moi, mlCAPE_dry, ml_CIN, mlLFC, mlEL, ml_CAPEACBL, mlELtooshort] = calcCAPE_MoistAndDry_ZAwrf2(Tv_4ml,revT_rho_lif_ml,h,Zo,mlLCLi,IsSat_JNM);   %chokes here
%note, this LFC is the reversible "moist LFC", which you will use henceforth

%     lfci = find(mlLFC >= h) ; lfci = lfci(end); mlLFCp = p(lfci)/100
%     eli = find(mlEL >= h) ; eli = eli(end); mlELp = p(eli)/100 

% moist,warm psuedoadiabatic, just for the cape for now:
[pseT_rho_lif_ml,IsSat_JNM] = lift_parcel_adiabatic_jpeters_ZAwrf2(T_4ml, pres*100, sphum_4ml, Zo, 0, 1 /  max( gradient(h) ), h, 3, 2);   %(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
[mlCAPEpse_moi, mlCAPEpse_dry, ml_CINpse, mlLFCpse, mlELpse, ml_CAPEACBLpse, mlELtooshortpse] = calcCAPE_MoistAndDry_ZAwrf2(Tv_4ml,pseT_rho_lif_ml,h,Zo,mlLCLi,IsSat_JNM);
clear blah2 blah3 blah4 blah5 blah6



%%%%%%%%%%
%%%%%%%%%%  DONE LIFTING PARCELS
%%%%%%%%%%










%{

% % % leave this commented out unless you REALLY want skew-T plots. They slow
% % % the code down a lot to generate. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% plot in skewT form  %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rd      = 287;              % Specific gas constant for dry air
    rv      = 461;              % Specific gas constant for water vapor
    cp      = 1004;             % Specific heat capacity for dry air
    cpv     = 1850;             % Specific heat capacity for water vapor
    cw      = 4200;             % heat capacity of water
    lv0     = 2.5e6;            % latent heat of vaporization at t0
    t0      = 273.15;           % reference point for latent heat of vaporization
    epsilon = 0.622;            % Rd/Rv
    one     = ones(1,1001);     % dummy array of ones
    f1 = figure;
    hold on;


    % Values of process lines desired on the diagram (allows you to make the
    % diagram less cluttered).

    %isobars = Pprof;

    isobars   = [100000 92500 90000 87500 85000 82500 80000 78500 77500 76500 75000 73800 72800 72500 72300 72000 71800 71500 71200 71000 70700 70500 70300 70000 68500 67500 65000 62500 60000 57500 50000 45000 40000 35000 32500 30000 29000 28000 26500 25000 23500 21500 20000 17500 15000 12500 10000]; % isobars to be plotted
    isobars = [100000, 92500, 90000, 87500, 85000, 82500, 80000, 79000:-250:65000, 62500 60000 57500 50000 45000 40000 35000 32500 30000 29000 28000 26500 25000 23500 21500 20000 17500 15000 12500 10000 ]
    isotherms = [-100 -90 -80 -70 -60 -50 -40 -30 -20 -10 0 5 10 15 20 25 30 40 50];          % isoterms to be plotted (degree C)
   % isotherms = [-60:1.5:25];          % isoterms to be plotted (degree C)
    theta     = [233. 253. 273. 283. 293. 197. 300. 301. 302. 303. 304. 305. 306. 307. 308. 309. 310. 311. 312. 313. 314. 315. 323. 333. 353. 373. 393. 413];                   % potential temperatures to be plotted (Kelvin)
   % theta     = [250:1.5:350];
    mixing    = [.2e-3 1.e-3 2.e-3 3.e-3 4.e-3 5.e-3 6.e-3 7.e-3 7.5e-3 8.e-3 8.25e-3 8.5e-3 8.75e-3 9.e-3 9.25e-3 9.5e-3 9.75e-3 10.e-3 11.e-3 12.e-3  14.e-3 16.e-3 18.e-3];                       % mixing ratio values to be plotted
   % mixing    = [0.2e-3 18.e-3];
    thetae    = [280. 300. 320. 340. 360. 380. 400.];
  %  thetae    = [280 400]; % moist adiabats to be plotted (Kelvin)
    %thetae    = [320.:5:360.]

    % define the yaxis from the highest and lowest pressures desired on the graph
    % for energy conservation purposes we define y = -Rd * ln(p)

    ptop    = 10000;                            %pressure level at top of diagram
    pbot    = 105000;                           %pressure level at bottom of
    ytop    = -rd*log(ptop);
    ybot    = -rd*log(pbot);
    dy      = (ytop-ybot)/1000;                 % I arbitrarily define 1000 points in both x and y (to get smooth lines)
    yaxis   = ybot:dy:ytop;


    % define the xaxis from the lowest (top left) and highest (bottom right) temperatures desired
    % the x coordinate is x = T + skew/Rd y

    twarm   = 41.5;                              % warmest temperature on the graph
    tcold   = -110.;                             % coldest temperature on the graph
    skew    = (tcold-twarm)/(2*log(pbot/ptop));  % skew the temperature line to 45 degree angle
    tmin    = tcold+273 - skew/rd*ytop;
    tmax    = twarm+273 - skew/rd*ybot;
    dt      = (tmax-tmin)/1000;                  % a thousand points again
    xaxis   = tmin:dt:tmax;
    axis([tmin tmax ybot ytop]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   END USER INPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % draw the isobars

    for q=1:length(isobars)
        pp      = -rd*log(isobars(q));              % calculate the y coordinate of isobar
        plot(xaxis,pp.*one,'k');                   % plot isobar
    end

    % label isobars

    set(gca,'TickLength',[0.0 0.0]);
    set(gca,'YTick',-rd*log(isobars));
    set(gca,'YTickLabel',{isobars/100});

    % draw the isotherms

    for q=1:length(isotherms)
        tt      = isotherms(q)+273 - skew/rd*yaxis;     % determine the x coordinate of the isotherms using the y coordinates
        indices = find(tt >= tmin & tt <= tmax);   % find the array indices that are actually on the plot
        plot(tt(indices),yaxis(indices),'k');
    end

    % label isotherms

    tlow    = tmin + skew/rd*ybot - 273;            % determine lowest temperature at 1050 mb
    indices = find(isotherms >= tlow & isotherms <= twarm);
    set(gca,'XTick',isotherms(indices)+273 - skew/rd*ybot);
    set(gca,'XTickLabel',{isotherms(indices)});

    % draw the potential temperature lines

    b0      = -rd*log(100000);            % the y coordinate for 1000 mb
    for q=1:length(theta)
        th      = theta(q)*exp((b0-yaxis)/cp) - skew/rd*yaxis;  % find the x coordinate for all the potential temperature lines
        indices = find(th >= tmin & th <= tmax);    % find the array indices that are actually on the plot
        plot(th(indices),yaxis(indices),'k');
    end

    % label potential temperature lines at the 115 mb level

    b115    = -rd*log(11500);
    th      = theta*exp((b0-b115)/cp)- skew/rd*b115;  % coordinates of the theta lines at 115 mb
    indices = find(th>= tmin & th <= tmax);
    for q=1:length(indices)
        text(th(indices(q)),b115,num2str(theta(indices(q))-273));
    end

    % draw the mixing ratio lines

    tt      = 213.0:0.1:313;                               % define temperature grid on which to plot the mixing ratio lines
    es      = 610.8*exp(53.486-6810./tt-5.09*log(tt));     % equilibrium vapor pressure from temperature ala Bohren page 192
    b205    = -rd*log(20500);                            % presuure coordinate for mixing ratio labels
    for q=1:length(mixing)
        pp      = (1. + epsilon/mixing(q))*es;        % pressure along saturation mixing ratio line
        xmix    = tt + skew*log(pp);                  % find the x-coordinate at these pressure levels
        ymix    = -rd*log(pp);                        % and the y-coordinate
        indices = find( (pp >= 20000 & pp <= pbot) & (xmix >= tmin & xmix <= tmax));  % constrain points to fit on specified plot
        plot(xmix(indices),ymix(indices),'--k');

        %label mixing ratio lines at 205 mb, or along the y-axis where they run of
        if (ymix(indices(1)) > b205)
            text(xmix(indices(1)), b205, num2str(mixing(q)*1000));
            %     else
            %         xmix(indices(1));, ymix(indices(1));
            %         text(xmix(indices(1)), ymix(indices(1)), num2str(mixing(i)*1000));
        end
    end


    % draw the moist adiabats - the native code version
    tt      = (tcold+273):(twarm-tcold)/1000:(twarm+273);% use the full range of temperature values
    es      = 610.8*exp(53.486-6810./tt-5.09*log(tt));   % equilibrium vapor pressure from temperature ala Bohren page 192
    lv      = lv0 + (cpv-cw)*(tt-t0);                    % latent heat of vaporation as function of temp ala bohren page 197
    for q=1:length(thetae)
        th      = thetae(q);
        x0      = pbot;                                  % set the first guess for the solution

        % solve for the pressure value at thetae given the temperature

        % attempting pseudo moist adiab
        for k=1:length(tt);
            pp(k)    = fzero(@fthetae,x0,optimset('Display','off'),...
                tt(k),es(k),lv(k)*epsilon*es(k)/cp/tt(k),th);
            if isfinite(pp(k))
                x0   = pp(k);
            end
        end

        xth     = tt + skew*log(pp);                % find the x-coordinate at these pressure levels
        yth     = -rd*log(pp);                      % and the y-coordinate
        indices = find( (pp >= ptop & pp <= pbot) & ...
            (xth >= tmin & xth <= tmax));  % constrain points to fit on specified plot
        plot(xth(indices),yth(indices),'k');
    end





    % draw the moist adiabats - the native code version
    tt      = (tcold+273):(twarm-tcold)/1000:(twarm+273);% use the full range of temperature values
    es      = 610.8*exp(53.486-6810./tt-5.09*log(tt));   % equilibrium vapor pressure from temperature ala Bohren page 192
    lv      = lv0 + (cpv-cw)*(tt-t0);                    % latent heat of vaporation as function of temp ala bohren page 197

    Ts = (2840./(3.5*log(tt) - log(es) -4.805)) - 55;
    mrs = 0.622*(es/(pbot - es));

    for q=1:length(thetae)
        th      = thetae(q);
        x0      = pbot;                                  % set the first guess for the solution

        % solve for the pressure value at thetae given the temperature

        % attempting pseudo moist adiab
        for k=1:length(tt);
            pp(k)    = fzero(@fthetaeTemp,x0,optimset('Display','off'),...
                tt(k),Ts(k),mrs,th);  %T,Ts,mr,th
            if isfinite(pp(k))
                x0   = pp(k);
            end
        end

        xth     = tt + skew*log(pp);                % find the x-coordinate at these pressure levels
        yth     = -rd*log(pp);                      % and the y-coordinate
        indices = find( (pp >= ptop & pp <= pbot) & ...
            (xth >= tmin & xth <= tmax));  % constrain points to fit on specified plot
        plot(xth(indices),yth(indices),'k');
    end





   % convert mix ratio to dewpoint
    % C=(0.622+mr);
    % A=(mr.*p)./C;
    % B=log(A/611);
    % TDprof = 1./(0.003663 - 0.0001844.*B) ;

    % eee = e/100; %vap press [mb]
    % TDprof =273.15+  243.5./( (17.67./log(eee./6.112)) -1);    %  deg K      K.Emanuel


    pp=pres*100;
    ttT=T;
    T(10)
    %uu=u;
    %vv=v;
    hgt=h;


    pp=pres*100;
    tt=Tv;
    Tv(10)
    %uu=u;
    %vv=v;
    hgt=h;

%         %Env. T:
%     yt      = -rd*log(pp);
%     xt      = ttT  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
%         isfinite(xt));
%     plot(xt(indices),yt(indices),'c','linewidth',1)
% 
%             %Env. T:
%     yt      = -rd*log(pp);
%     xt      = tt  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
%         isfinite(xt));
%     plot(xt(indices),yt(indices),'y','linewidth',1)

%     Tv = Tv_preserve;
%     T = T_preserve;
%     TDc = TDc_preserve;
%     mr = mr_preserve;
%     thetav = thetav_preserve;
%     theta = theta_preserve;
%     e = e_preserve;
%     RH = RH_preserve;
%     esat = esat_preserve;

%    TDprof = TDc + 273.15;

%     pp=pres*100;
%     ttT=T;
%     T(10)
%     %uu=u;
%     %vv=v;
%     hgt=h;



%     %Env. T:
%     yt      = -rd*log(pp);
%     xt      = ttT  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
%         isfinite(xt));
%     plot(xt(indices),yt(indices),'k','linewidth',2)



%     %Env. Td:
%     ytd      = -rd*log(pp);
%     xtd     = (TDc_4ml + 273.15) + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtd >= tmin & xtd <= tmax) &...
%         isfinite(xtd));
%     plot(xtd(indices),ytd(indices),'co','linewidth',1)
% 
% 
%   %Env. Tt:
%     yt      = -rd*log(pp);
%     xt      = Tv_4ml  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
%         isfinite(xt));
%     plot(xt(indices),yt(indices),'go','linewidth',3)
% 
%    %Env. Tv:
%     yt      = -rd*log(pp);
%     xt      = T_4ml  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
%         isfinite(xt));
%     plot(xt(indices),yt(indices),'co','linewidth',3)





    Tv = Tv_preserve;
    T = T_preserve;
    TDc = TDc_preserve;
    mr = mr_preserve;
    thetav = thetav_preserve;
    theta = theta_preserve;
    e = e_preserve;
    RH = RH_preserve;
    esat = esat_preserve;


    %Env. Td:
    ytd      = -rd*log(pp);
    xtd     = (TDc + 273.15) + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xtd >= tmin & xtd <= tmax) &...
        isfinite(xtd));
    plot(xtd(indices),ytd(indices),'b','linewidth',1)


  %Env. Tt:
    yt      = -rd*log(pp);
    xt      = T  + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
        isfinite(xt));
    plot(xt(indices),yt(indices),'b','linewidth',3)

  %Env. Tt:
    yt      = -rd*log(pp);
    xt      = Tv  + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xt >= tmin & xt <= tmax) &...
        isfinite(xt));
    plot(xt(indices),yt(indices),'b','linewidth',3)



%     %Env. Td:
%     ytd      = -rd*log(pp);
%     xtd     = TDprof + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtd >= tmin & xtd <= tmax) &...
%         isfinite(xtd));
%     plot(xtd(indices),ytd(indices),':k','linewidth',3)
%     
%
%     %parcel VIRTUAL Temperature
%     ytp      = -rd*log(pp);
%     xtp      = parcelTprofR'  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp >= tmin & xtp <= tmax) &...
%         isfinite(xtp));
%     plot(xtp(indices),ytp(indices),'b','linewidth',0.5)
%         %parcel VIRTUAL Temperature
%     ytp      = -rd*log(pp);
%     xtp      = parcelTprofP'  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp >= tmin & xtp <= tmax) &...
%         isfinite(xtp));
%     plot(xtp(indices),ytp(indices),'m','linewidth',0.5)


%        %parcel VIRTUAL Temperature
%     ytp1      = -rd*log(pp);
%     xtp1      =     parcelTVprof_zos(1,:)'  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
%         isfinite(xtp1));
%     plot(xtp1(indices),ytp1(indices),'Color',[1 0.0 0],'linewidth',1)
% 
%        %parcel VIRTUAL Temperature
%     ytp1      = -rd*log(pp);
%     xtp1      =     parcelTVprof_zos(2,:)'  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
%         isfinite(xtp1));
%     plot(xtp1(indices),ytp1(indices),'Color',[1 0.4 0],'linewidth',1)
%     
%         %parcel VIRTUAL Temperature
%     ytp1      = -rd*log(pp);
%     xtp1      =  parcelTVprof_zos(3,:)'  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
%         isfinite(xtp1));
%     plot(xtp1(indices),ytp1(indices),'Color',[1 0.5 0],'linewidth',1)
    


       %parcel VIRTUAL Temperature
    ytp1      = -rd*log(pp);
    xtp1      = revT_rho_lif_sfc  + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
        isfinite(xtp1));
    plot(xtp1(indices),ytp1(indices),'r','linewidth',1)


    % % revT_rho_lif_mu(1:14) = NaN;

       %parcel VIRTUAL Temperature
    ytp1      = -rd*log(pp);
    xtp1      = revT_rho_lif_mu  + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
        isfinite(xtp1));
    plot(xtp1(indices),ytp1(indices),'b','linewidth',1)


       %parcel VIRTUAL Temperature
    ytp1      = -rd*log(pp);
    xtp1      = revT_rho_lif_ml  + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
        isfinite(xtp1));
    plot(xtp1(indices),ytp1(indices),'g','linewidth',1)
    
    
       %parcel VIRTUAL Temperature
    ytp1      = -rd*log(pp);
    xtp1      = pseT_rho_lif_mu  + skew*log(pp);
    indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
        isfinite(xtp1));
    plot(xtp1(indices),ytp1(indices),'g-.','linewidth',1)
    
    
    
%     
%            %parcel VIRTUAL Temperature
%     ytp1      = -rd*log(pp);
%     xtp1      = pseT_rho_lif_sfc  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
%         isfinite(xtp1));
%     plot(xtp1(indices),ytp1(indices),'r:','linewidth',1)


%        %parcel VIRTUAL Temperature
%     ytp1      = -rd*log(pp);
%     xtp1      = pseT_rho_lif_mu  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
%         isfinite(xtp1));
%     plot(xtp1(indices),ytp1(indices),'g:','linewidth',1)
% 
% 
%        %parcel VIRTUAL Temperature
%     ytp1      = -rd*log(pp);
%     xtp1      = pseT_rho_lif_ml  + skew*log(pp);
%     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
%         isfinite(xtp1));
%     plot(xtp1(indices),ytp1(indices),'c:','linewidth',1)
    
    

% 
% 
% %            %parcel VIRTUAL Temperature
% %     ytp1      = -rd*log(pp);
% %     xtp1      = parcelTprof_sfc'  + skew*log(pp);
% %     indices = find( (pp >= ptop & pp <= pbot) & (xtp1 >= tmin & xtp1 <= tmax) &...
% %         isfinite(xtp1));
% %     plot(xtp1(indices),ytp1(indices),'r','linewidth',0.5)
% % 
% % 
% %        %parcel VIRTUAL Temperature
% %     ytp2      = -rd*log(pp);
% %     xtp2      = parcelTprof_mu'  + skew*log(pp);
% %     indices = find( (pp >= ptop & pp <= pbot) & (xtp2 >= tmin & xtp2 <= tmax) &...
% %         isfinite(xtp2));
% %     plot(xtp2(indices),ytp2(indices),'b--','linewidth',0.5)
% % 
% % 
% %     %parcel VIRTUAL Temperature .
% %     ytp3      = -rd*log(pp);
% %     xtp3      = parcelTprof_ml  + skew*log(pp);
% %     indices = find( (pp >= ptop & pp <= pbot) & (xtp3 >= tmin & xtp3 <= tmax) &...
% %         isfinite(xtp3));
% %     plot(xtp3(indices),ytp3(indices),'g:','linewidth',0.5)
% 
% 
% 
% %     parcelTVprof_sfc
% %     parcelTVprof_mu
% %     parcelTVprof_mu
% 
% 
     axis([-133.0704398539077,-64.47088233178378,-3285.87077169578,-2684.877002193151])
%    % axis([-78.75194863893432,-71.4654639906638,-3270.24386930171,-3217.533902173995])
% 
% 
%     skewttitle = horzcat(num2str(YYMMDD),' ',num2str(TTT),' ',hhh,mmm,sss) 
%     title([ skewttitle])
% 
%     %xlabel(['Tv Profile: cape (LFC-EL) = ', num2str(round(CAPEtv)), ';  cin (Zo-LFC) = ', num2str(round(CINtv)),';  LFC(mb) =', num2str(round(LFCptv)),';  LCL =', num2str(round(LCLptv)), ', EL =', num2str(round(ELptv)) ])
% 
%    % outdir = '/Users/marq789/Documents/PROJECTS/ICLASS/CACTI-CI/interpsondes/testskewT/';
%     outdir = '/Users/marq789/Downloads/testskewT/';
%     %outdir = '/Users/marq789/Downloads/mergeskewT/';
%     mkdir(outdir)
% 
%     skewtlab = horzcat(num2str(YYMMDD),'_interpsonde_SkewT_',num2str(TTT),'_',hhh,mmm,sss,'.jpg')
%     %skewtlab = horzcat(num2str(YYMMDD),'_mergesonde_SkewT_',num2str(TTT),'_',hhh,mmm,sss,'.jpg');
%     
%     saveas(f1,horzcat(outdir,skewtlab))
% 
%     clear close
%     blah = get(groot,'CurrentFigure');
%     close(blah);

%}
