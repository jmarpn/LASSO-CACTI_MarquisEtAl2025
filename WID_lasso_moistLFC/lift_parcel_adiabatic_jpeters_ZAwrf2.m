function [T_rho_lif,IsSat_JNM]=lift_parcel_adiabatic_jpeters_ZAwrf2(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)
%function [T_lif,T_rho_lif,Tv_lif,Qv_lif,Qt_lif,B_lif]=lift_parcel_adiabatic_jpeters(T0,p0,q0,start_loc,fracent,prate,z0,T1,T2)


%JNMdebuggers
%   T0 = T;   p0 = pres*100;  q0 = mr;  start_loc = Zo;  fracent =0;  prate = 0;  z0 = h;  T1 = 3;  T2 = 2;





    %this function computes lifted parcel properties using the unsaturated
    %and saturated lapse rate formulas from (Peters et al. 2022)
    %https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 
    
    %input arguments
    %T0: sounding profile of temperature (in K)
    %p0: sounding profile of pressure (in Pa)
    %q0: sounding profile of water vapor mass fraction (in kg/kg)
    %start_loc: index of the parcel starting location (set to 1 for the
    %lowest: level in the sounding)
    %fracent: fractional entrainment rate (in m^-1)
    
    %output arguments
    %T_lif: lifted parcel temperature
    %Qv_lif: lifted parcel water vapor mass fraction
    %Qt_lif: lifted parcel total water mass fraction
    %B_lif: Lifted parcel buoyancy, computed using Eq. B6 in (Peters et al.
    %2022) (accounts for virtual temperature and loading effects)
    
    %IsSat_JNM       %JNM: creating vert log of saturated or unsat (0 = dry, 1 = sat)

    %prate: precipitation rate (in m^-1) large values make parcel more
    %pseudoadiabatic, small values make parcel more adiabatic.  This parameter is equal to the inverse of "L" used in the text, and its units are m^{-1}
    
    %z0: sounding profile of height above ground level (first height should
    %be 0 m)
    %T1 warmest mixed-phase temperature
    %T2 coldest mixed-phase temperature

    %CONSTANTS
    Rd=287.04; %dry gas constant
    Rv=461.5; %water vapor gas constant
    epsilon=Rd./Rv;   
    g=9.81; %gravitational acceleration

    % 
    %descriminator function between liquid and ice (i.e., omega defined in the
    %beginning of section 2e in Peters et al. 2022)
    omega_ = @(T) min(max((T - T1)./(T2-T1),0),1);
    %this is the height derivative of omega, which is also needed in the lapse
    %rate formula
    domega_ = @(T) (heaviside(T1-T) - heaviside(T2-T))./(T2-T1);
    
    T_lif=zeros(size(T0)); %temperature of the lifted parcel
    Qv_lif=zeros(size(T0)); %water vapor mass fraction of the lifted parcel
    Qt_lif=zeros(size(T0)); %total water mass fraction of the lifted parcel
    T_lif(1:1:start_loc)=T0(1:1:start_loc); %set initial values to that of the environment
    Qv_lif(1:1:start_loc)=q0(1:1:start_loc); %set initial values to that of the environment
    Qt_lif=Qv_lif; %set initial values to that of the environment

    IsSat_JNM = zeros(size(T0));    IsSat_JNM(:) = NaN;    %JNM: creating vert log of saturated or unsat (0 = dry, 1 = sat)
    Qs_lif_JNM = NaN;
     
    q_sat_prev=0;
    for iz=start_loc+1:length(z0)
        q_sat=(1-Qt_lif(iz-1)).*compute_rsat(T_lif(iz-1),p0(iz-1),1,T1,T2);
        
        Qs_lif_JNM = vertcat(Qs_lif_JNM,q_sat);

        if Qv_lif(iz-1)<q_sat %if we are unsaturated, go up at the unsaturated adiabatic lapse rate (eq. 19 in Peters et al. 2022)
            
            q0(iz-1);
            T_lif(iz) = T_lif(iz-1) + (z0(iz) - z0(iz-1))*drylift(T_lif(iz-1),Qv_lif(iz-1),T0(iz-1),q0(iz-1),fracent);
            Qv_lif(iz) = Qv_lif(iz-1) - (z0(iz) - z0(iz-1))*fracent*( Qv_lif(iz-1) - q0(iz-1) );
            Qt_lif(iz) = Qv_lif(iz);
            q_sat=(1-Qt_lif(iz)).*compute_rsat(T_lif(iz),p0(iz),1,T1,T2);

            if Qv_lif(iz)<q_sat
               IsSat_JNM(iz) = 0;
            end

            if Qv_lif(iz)>=q_sat %if we hit saturation, split the vertical step into two stages.  The first stage advances at the saturated lapse rate to the saturation point, and the second stage completes the grid step at the moist lapse rate
   
                IsSat_JNM(iz) = 1;
                
                satrat=(Qv_lif(iz)-q_sat_prev)./(q_sat-q_sat_prev);
                dz_dry=satrat*(z0(iz)-z0(iz-1));
                dz_wet=(1-satrat)*(z0(iz)-z0(iz-1));

                T_halfstep = T_lif(iz-1) + dz_dry*drylift(T_lif(iz-1),Qv_lif(iz-1),T0(iz-1),q0(iz-1),fracent);
                Qv_halfstep = Qv_lif(iz-1) - dz_dry*fracent*( Qv_lif(iz-1) - q0(iz-1) );
                Qt_halfstep = Qv_lif(iz);
                p_halfstep=p0(iz-1)*satrat + p0(iz)*(1-satrat);
                T0_halfstep=T0(iz-1)*satrat + T0(iz)*(1-satrat);
                Q0_halfstep=q0(iz-1)*satrat + q0(iz)*(1-satrat);

                T_lif(iz) = T_halfstep + dz_wet*moislif(T_halfstep,Qv_halfstep,...
                    (1-Qt_halfstep).*compute_rsat(T_halfstep,p_halfstep,0,T1,T2),(1-Qt_halfstep).*compute_rsat(T_halfstep,p_halfstep,2,T1,T2),p_halfstep,T0_halfstep,Q0_halfstep...
                    ,omega_(T_halfstep),domega_(T_halfstep),Qt_halfstep,fracent,prate);
                
                qent_ref=Qt_halfstep;
                
                Qt_lif(iz) = Qt_lif(iz-1) - (z0(iz) - z0(iz-1))*fracent*( qent_ref - Q0_halfstep );
                Qv_lif(iz) = (1-Qt_lif(iz)).*compute_rsat(T_lif(iz),p0(iz),1,T1,T2);

                if Qt_lif(iz)<Qv_lif(iz)
                    Qv_lif(iz)=Qt_lif(iz);
                end

                %
            end
            q_sat_prev=q_sat;
            %
        else %if we are already at saturation, just advance upward using the saturated lapse rate (eq. 24 in Peters et al. 2022)
            %
            IsSat_JNM(iz) = 1;

            prate_=prate; %this is switch really has no purpose, it's a relic of some old things i did with the code

            T_lif(iz) = T_lif(iz-1) + (z0(iz) - z0(iz-1))*moislif(T_lif(iz-1),Qv_lif(iz-1),...
                (1-Qt_lif(iz-1)).*compute_rsat(T_lif(iz-1),p0(iz-1),0,T1,T2),(1-Qt_lif(iz-1)).*compute_rsat(T_lif(iz-1),p0(iz-1),2,T1,T2),p0(iz-1),T0(iz-1),q0(iz-1)...
                ,omega_(T_lif(iz-1)),domega_(T_lif(iz-1)),Qt_lif(iz-1),fracent,prate_);
                     
            qent_ref=Qt_lif(iz-1); %this is switch really has no purpose, it's a relic of some old things i did with the code
                 
            Qt_lif(iz) = Qt_lif(iz-1) - (z0(iz) - z0(iz-1))*(fracent*( qent_ref - q0(iz-1) )  + prate_*( Qt_lif(iz-1)-Qv_lif(iz-1)) );
            Qv_lif(iz) = (1-Qt_lif(iz)).*compute_rsat(T_lif(iz),p0(iz),1,T1,T2);
            
            if Qt_lif(iz)<Qv_lif(iz)
                Qv_lif(iz)=Qt_lif(iz);
            end
        end
    end
        
    T_rho_lif=T_lif.*(1 - Qt_lif + Qv_lif)./( 1 + (epsilon - 1)./( ( epsilon.*(1 - Qt_lif)./Qv_lif - 1) ) );
    

    %QvQsRatio_lif = Qv_lif ./ Qs_lif_JNM ;



   % Tv_lif = T_lif.*( 1 +  Qv_lif*0.61);  %jnm added
    
    %this line is John's backwards way of calculating Trho. He said you can
    %replace it with T_0_lif = T0.*( 1 + (Rv/Rd - 1)*q0 ). Trho of env =
    %Tv):
%    T_0_lif=T0./(    1 + (epsilon - 1)./ ( ( epsilon.*(1 - q0)./q0 - 1) ) );  %


%     B_lif=g.*(T_rho_lif - T_0_lif)./T_0_lif;



end





%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

function gamma_m=moislif(T,qv,qvv,qvi,p0,T0,q0,omega,domega,qt,fracent,prate)
    
    %CONSTANTS
    Rd=287.04; %dry gas constant
    Rv=461.5; %water vapor gas constant
    epsilon=Rd./Rv;
    cp=1005; %specific heat of dry air at constant pressure
    g=9.81; %gravitational acceleration
    xlv=2501000; %reference latent heat of vaporization at the triple point temperature
    xls=2834000; %reference latent heat of sublimation at the triple point temperature
    cpv=1870; %specific heat of water vapor at constant pressure
    cpl=4190; %specific heat of liquid water
    cpi=2106; %specific heat of ice
    ttrip=273.15; %triple point temperature
 
    if qt<0
        zeroB=0;
        qt=-qt;
    else
        zeroB=1;
    end
    
    if qv<0
        noT=true;
        qv=-qv;
    else
        noT=false;
    end
    
    cpm = (1 - qt)*cp + qv*cpv + (1 - omega)*(qt-qv)*cpl + omega*(qt-qv)*cpi;
    Lv = xlv + (T - ttrip)*(cpv - cpl);
    Li = (xls-xlv) + (T - ttrip)*(cpl - cpi); 
    Rm0 = (1 - q0)*Rd + q0*Rv;
    
    if qt==qv
        pseudofac=(1-qv);
    else
        pseudofac=1;
    end
    

    
    T_rho=T*(1 - qt + qv/epsilon);
    T_rho0=T0.*( 1 - q0 + q0/epsilon );
    B = zeroB.*g*(T_rho - T_rho0)./(T_rho0);
    
    Qvsl = qvv/( epsilon - epsilon*qt + qv);
    Qvsi = qvi/( epsilon - epsilon*qt + qv);
    Q_M = (1 - omega)*qvv/(1 - Qvsl) + omega*qvi/(1 - Qvsi);
    L_M = Lv*(1 - omega)*qvv/(1 - Qvsl) + (Lv + Li)*omega*qvi/(1 - Qvsi);
    %
    
    eps_T = -fracent*(T - T0);
    eps_qv = -fracent*(qv - q0);
    eps_qt = -fracent*(qt - q0)-prate*(qt-qv);
    term1 = -B;
    if noT
        term2 = - pseudofac*Q_M*(Lv + Li*omega)*g/(Rd*T_rho);
    else
        term2 = - pseudofac*Q_M*(Lv + Li*omega)*g/(Rm0*T0);
    end
    term3 = -g;
    term4 = (cpm - Li*(qt - qv)*domega)*eps_T;
    term5 = (Lv + Li*omega)*(eps_qv + (qv/(1-qt))*eps_qt);
    %
    term6 = cpm;
    term7 = -Li*(qt - qv)*domega;
    term8 = pseudofac*(Lv + Li*omega)*(-domega*(qvv - qvi) + (1/(Rv*T^2))*(L_M));
    gamma_m =( term1 + term2 + term3 + term4 + term5)/(term6 + term7 + term8);
end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================

function gamma_d=drylift(T,qv,T0,qv0,fracent)
    %CONSTANTS
    Rd=287.04; %dry gas constant
    Rv=461.5; %water vapor gas constant
    cp=1005; %specific heat of dry air at constant pressure
    g=9.81; %gravitational acceleration
    cpv=1870; %specific heat of water vapor at constant pressure
    
    cpmv = (1 - qv)*cp + qv*cpv;
    B = g*( (T-T0)/T0 + (Rv/Rd - 1)*(qv - qv0) );
    eps = -fracent*(T - T0);
    gamma_d = - (g + B)/cpmv + eps;
end

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================



function [qsat] = compute_rsat(T,p,iceflag,T1,T2)

    %THIS FUNCTION COMPUTES THE SATURATION MIXING RATIO, USING THE INTEGRATED
    %CLAUSIUS CLAPEYRON EQUATION (eq. 7-12 in Peters et al. 2022).
    %https://doi-org.ezaccess.libraries.psu.edu/10.1175/JAS-D-21-0118.1 

    %input arguments
    %T temperature (in K)
    %p pressure (in Pa)
    %iceflag (give mixing ratio with respect to liquid (0), combo liquid and
    %ice (2), or ice (3)
    %T1 warmest mixed-phase temperature
    %T2 coldest mixed-phase temperature
    
    %NOTE: most of my scripts and functions that use this function need
    %saturation mass fraction qs, not saturation mixing ratio rs.  To get
    %qs from rs, use the formula qs = (1 - qt)*rs, where qt is the total
    %water mass fraction

    %CONSTANTS
    Rd=287.04; %dry gas constant
    Rv=461.5; %water vapor gas constant
    epsilon=Rd./Rv;
    cp=1005; %specific heat of dry air at constant pressure
    g=9.81; %gravitational acceleration
    xlv=2501000; %reference latent heat of vaporization at the triple point temperature
    xls=2834000; %reference latent heat of sublimation at the triple point temperature
    cpv=1870; %specific heat of water vapor at constant pressure
    cpl=4190; %specific heat of liquid water
    cpi=2106; %specific heat of ice
    ttrip=273.15; %triple point temperature
    eref=611.2; %reference pressure at the triple point temperature

    %descriminator function between liquid and ice (i.e., omega defined in the
    %beginning of section 2e in Peters et al. 2022)
    omeg = ((T - T1)./(T2-T1)).*heaviside((T - T1)./(T2-T1)).*heaviside((1 - (T - T1)./(T2-T1))) + heaviside(-(1 - (T - T1)./(T2-T1)));

    switch iceflag
        case 0 %only give mixing ratio with respect to liquid
             term1=(cpv-cpl)/Rv;
             term2=(xlv-ttrip*(cpv-cpl))/Rv;
             esl=exp((T-ttrip).*term2./(T.*ttrip)).*eref.*(T./ttrip).^(term1);
             esl = min( esl , p*0.5 );
             qsat=epsilon.*esl./(p-esl);         
        case 1 %give linear combination of mixing ratio with respect to liquid and ice (eq. 20 in Peters et al. 2022)
             term1=(cpv-cpl)/Rv;
             term2=(xlv-ttrip*(cpv-cpl))/Rv;
             esl_l=exp((T-ttrip).*term2./(T.*ttrip)).*eref.*(T./ttrip).^(term1);

             qsat_l=epsilon.*esl_l./(p-esl_l);

            term1=(cpv-cpi)/Rv;
             term2=( xls-ttrip*(cpv-cpi))/Rv;
             esl_i=exp((T-ttrip).*term2./(T.*ttrip)).*eref.*(T./ttrip).^(term1);

             qsat_i=epsilon.*esl_i./(p-esl_i);

             qsat=(1-omeg).*qsat_l + (omeg).*qsat_i;
        case 2 %only give mixing ratio with respect to ice
            term1=(cpv-cpi)/Rv;
             term2=( xls-ttrip*(cpv-cpi))/Rv;
             esl=exp((T-ttrip).*term2./(T.*ttrip)).*eref.*(T./ttrip).^(term1);
             esl = min( esl , p*0.5 );
             qsat=epsilon.*esl./(p-esl);
    end


end

