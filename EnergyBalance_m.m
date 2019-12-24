function [ustar,H, LE, G0, H_DL, H_WL, H_i,evap_fr,LEp_PM,LEp_PT09,LEp_PT126,LE0_fao56]= EnergyBalance_m(d0, z0m, z0h, fc_V,LAI_V, ..., 
                                                                            Rn, LST_V,...
                                                                            hpbl, Zref_V, Tref_K_V, Uref_V, Earef_V,qaref_V, Pref_V, Ps,G0)

                                                                        
% NOTES
% Rn=Rn(i); LST_K=LST_K(i);hpbl=Hi_PBL;
% Tref_K=Tref_K(i); Uref=Uref(i); Earef=Earef(i);qaref=qaref(i);Pref=Pref(i);G0=G(i);Ps=Pref;
%---------------------------------------------------------------------
% Description of parameters

% d0        -> Zero plane displacement height           (m)
% z0m       -> Roughness height for momentum tranfer    (m)
% z0h       -> Roughness height for heat tranfer        (m)
% fc        -> Fractional vegetaion cover               (-)
% LAI: canopy total leaf area index
% Rn:       -> net radiation                            (Watt/m^2)
% LST_K     -> Surface temperature                      (K)
% hpbl       -> Height of the PBL                        (m)
% Zref      -> Reference height                         (m)
% Uref      -> Wind speed at reference height           (m/s)
% Tref_K    -> Air temperature at reference height      (K)
% Pref      -> Pressure at reference height             (Pa)
% qaref     -> Specific humidity at reference height    (kg/kg)

% Ps        -> Surafce pressure                         (Pa)
% SWd       -> Downward Solar Radiation                 (Watt/m^2)
% LWd       -> Downward long wave radiation             (Watt/m^2)",
% albedo    -> Albedo                                   (-)

%---------------------------------------------------------------------
% Solving of 3 equations
% Nots: Here we start to solve the system of three equations
% i.e. the equation of momentum transfer, the equation of heat transfer, and the equation for the stability length.
% We use ASL functions of Brutsaert, 1999, if Zref < hst, the ASL height
% We use BAS functions of Brutsaert, 1999, if Zref > hst, Zref <= hpbl, the PBL height.

% Note: We will determine the Monin-Obukhov length using the definition
% given by Brutsaert (1982), P.65, eqn. 5.25.
% i.e. L = - ustar^3*rhoa/(k*g*((H/Ta*Cp)+0.61*E))
% By using the quantity, ef=Le*E/(Rn-G0), we can write
% H = (1-ef)*(Rn-G0) and E = ef/Le*(Rn-G0)
% So that L =-ustar^3*rhoam/(k*g*((1-ef)/(Ta*CP)+0.61*ef/Le)*(Rn-G0))
% From this eqn., it is obvious that L=f(ustar^3) and L=f(ef^-1)


% LIMITING CASES: A good idea could be to take L_DL and L_WL respectively as
% the max and min stability length.
% This has the advantage that the feedback effects between land
% surface and the hpbl are considered.
% For theoretical limits, H=Rn-G0, and E=Rn-G0 respectively.
% (Note: For wet surfaces in a dry climate, E may be bigger than
% Rn-G0 due to negative H,
% i.e. the Oasis effect. In this case, a small modification of
% L_WL can be expected. This is ignored for the time being)
% Previously, we interpolated mu_i between mu_i(0) and
% mu_i(-1500), by means of temperature difference
% Though other factors may aslo have influences as seen in the
% calculation of T0Ta_l and T0Ta_u,
% the uncertainties associated to those factors do not warrant
% their adequate applications for this case.
% This is consistant with the definition of resistences for the
% limiting cases, that is, that we ignore the stable case
% over the whole region.

%% 
constants
global P0 g k  Cpw Cpd L_e 
global Rd gamma

%% Meteorological Parameters
% Earef                               = 	Pref .* qaref * (Rv / Rd);                          % actual vapour pressure (based on Pressure at reference height)

Theta_a                             =   Tref_K_V .*((P0./Pref_V ).^0.286);                      % potential Air temperature    Eq. 2.23 P32 Brutsaert 2008
Theta_s                             =   LST_V  .*((P0./Ps   ).^0.286);                      % potential surface temperature

Theta_av                            =   Theta_a.* (1 + 0.61 * qaref_V);                       % virtual potential air temperature  P32 Brutsaert 2008 

% air densities (Eq 2.4,2.6 P25 Brutsaert 2008)
rhoa_m                              =   (Pref_V - 0.378 * Earef_V)./ (Rd * Tref_K_V);              % density of moist air [kg/m3]
rhoa_WL                             =   (Pref_V - 1.000 * Earef_V)./ (Rd * LST_V);               % density of dry air. [kg/m3]
% NOTE: rhoa_WL is only used for the wet-limit. To get a true upperlimit for the sensible heat
% the Landsurface Temperature is used as a proxy instead of air temperature.

% moist air density (kg / m3), this the same as rhoa_m
Cp                                  =   qaref_V* Cpw + (1-qaref_V)*Cpd;                         % Specific heat for constant pressure Page 29 of Brutsaert 2005
rhoa_m_Cp                           =   rhoa_m .* Cp;                                        % specific air heat capacity (J K-1 m-3)
rhoa_WL_Cp                          =   rhoa_WL .* Cp;                                       % specific air heat capacity (J K-1 m-3)

% Computing saturated vapor pressure for land surface at wet limit (WL) (this is for lower limit of PBL)
LST_C                               =   LST_V - 273.15;                                     % Landsurface temperature[C].
A                                   =   611;                                                % [Pa]
B                                   =   17.502;
C                                   =   240.97;                                             % [C]
esat_WL                             =   A * exp(B * LST_C./(LST_C + C));                    % Pa,(3.8),p.41 of Campbell & Norman, 1998
slope_WL                            =   B * C * esat_WL ./ ((C + LST_C).^2);                % Pa*0C-1,(3.9)

% NOTE: esat_WL is only used for the wet-limit. To get a true upperlimit for the sensible heat
% the Landsurface Temperature is used as a proxy instead of air temperature.
%% Rn and G0
% Rn                                  =   (1.0 - albedo) .* SWd + emissivity .* LWd - emissivity .* Sigma_SB .* LST_K.^4;

% global load_G0_ i ioverpass SCOPE_output
% if load_G0_
%     G0                              =   SCOPE_output.G0(ioverpass(i));
% else
%     Lambda_s                        =	0.315;                                              % bare soil (Kustas et al., 1989)
%     Lambda_c                        =	0.050;                                              % full vegetation canopy (Monteith, 1973)
%     G0                              =   Rn .* (Lambda_c + (1 - fc) * (Lambda_s - Lambda_c));
% end

%% ASL height
% hst= alfa*hpbl, with alfa=0.10~0.15 over moderately rough surfaces, or hst=beta*z0, with beta= 100~150.

alfa                                =   0.12;                                               % These are mid values as given by Brutsaert,1999
beta                                =   125;                                                % These are mid values as given by Brutsaert,1999
hst                                 =   max(alfa * hpbl, beta * z0m);                       % height of ASL


%% U* and L

CH                                  =   (Theta_s - Theta_a) .* k .* rhoa_m_Cp;              % from Eq. 2.55 P47 of Brutsaert 2008
CL                                  =   -rhoa_m_Cp .* Theta_av/ (k * g);                    % in this formula, air virtual potential temperature

z_d0                                =   Zref_V - d0;
ku                                  =   k * Uref_V;  
log_z_d0_z0m                        =   log(z_d0 ./ z0m);
log_z_d0_z0h                        =   log(z_d0 ./ z0h);

% Initial guess for u*, H and L assuming neutral stability
L                                   =   0;                                                  % initial L is zero for neutral condition
ustar                               =   ku ./ log_z_d0_z0m;                                 % U* in neutral condition when stability factors are zero
H                                   =   CH .* ustar ./ log_z_d0_z0h;                        % H  in neutral condition when stability factors are zero
errorH                              =   10;
H0                                  =   H;                                                  % H0 is H in neutral condition
steps                               =   0;
tic
if (Zref_V <= hst)
    while   (max(errorH(:)) > 0.1 && steps < 100)                                                  % MOS        
        L                           =   CL .* (ustar.^3) ./ H;
        psim1=PSIm(z_d0./L);
        psim2=PSIm(z0m./L);
        ustar                       =   ku ./ (log_z_d0_z0m - psim1 + psim2);
        psih1=PSIh(z_d0./L);
        psih2=PSIh(z0h./L);
        H                           =   CH .* ustar ./ (log_z_d0_z0h - psih1 + psih2);
        errorH                      =   abs(H0 - H);
        H0                          =   H;
        steps                       =   steps + 1  
        A=zeros(10,10)
    end
    C_i1                            =   PSIh(Zref_V./ L);
    C_i2                            =   PSIh(z0h ./ L);
else
   while    (max(errorH(:)) > 0.1 && steps < 100)                                                  % BAS
        L                           =   CL .* (ustar.^3) ./ H;                              % this is temporary L computed based on previous step Ustar and H
        ustar                       =   ku ./ (log_z_d0_z0m - Bw(hpbl, L, z0m, d0));        % Eq. 2.67 P 52 Brutsaert 2008
        H                           =   CH .* ustar./(log_z_d0_z0h - Cw(hpbl,L,z0m,z0h,d0));% Eq. 2.68 P 52 Brutsaert 2008
        errorH                      =   abs(H0 - H);
        H0                          =   H;
        steps                       =   steps + 1
        A=ones(10,10)
   end
   C_i1                             =   Cw(Zref_V, L, z0m, z0h, d0);
   C_i2                             =   0;
end
toc
%% Sensible heat Flux
I_1                                 =   log_z_d0_z0h + C_i2 >  C_i1;
I_2                                 =   log_z_d0_z0h + C_i2 <= C_i1;

re_i1                               =   (log_z_d0_z0h - C_i1 + C_i2)./(k * ustar);          % Actual resistance to heat transfer (s/m)
re_i2                               =   (log_z_d0_z0h      )./(k * ustar);
re_i                                =   I_1.*re_i1 + I_2.*re_i2;

H_i                                 =   rhoa_m_Cp.* (Theta_s-Theta_a)./re_i;                % Sensible heat flux
H_i                                 =   real (H_i);


%% Sensible heat flux at theoretical Dry limit
% Dry limit
%L_dry                              = 	0;
H_DL                                =   Rn - G0;
H_DL                                =   real (H_DL);
%% Sensible heat flux at theoretical wet limit
% Dry air is assumed.. eact=0, and we need to take the density of dry air
L_WL                                =   -(ustar.^3).* rhoa_WL./(k*g*(0.61* (Rn - G0)/ L_e)); 

% Bulk Stability Corrections
I_MOS                               =   (Zref_V < hst);
I_BAS                               =   (Zref_V >= hst);


C_WL_MOS                            =   PSIh(-Zref_V ./ L_WL);
C_WL_BAS                            =   Cw(Zref_V, L_WL, z0m, z0h,d0);
C_WL                                =   I_MOS.*C_WL_MOS + I_BAS.*C_WL_BAS;

% Calculating Resistances
I_1                                 =   log_z_d0_z0h>C_WL;
I_2                                 =   log_z_d0_z0h<=C_WL;

re_WL1                              =   (log_z_d0_z0h - C_WL)./(k * ustar);                 % Actual resistance to heat transfer at wet limit(s/m)
re_WL2                              =   (log_z_d0_z0h       )./(k * ustar);
re_WL                               =   I_1.*re_WL1 + I_2.*re_WL2;

H_WL                                =   ((Rn - G0) - (rhoa_WL_Cp./re_WL).*((esat_WL)/ gamma))./(1.0 + slope_WL/gamma);
H_WL                                =   real (H_WL);
%% Evaporative fraction 
H_i                                 =   min(H_i, H_DL);
H_i                                 =   max(H_i, H_WL);

% Relative evaporation
I_water                             =   (H_DL <= H_WL);
I_land                              =   (H_DL > H_WL);

evap_re_water                       =   1;                                                  % for water & wet surfaces
evap_re_land                        =   1 - (H_i - H_WL) ./ (H_DL - H_WL);
evap_re                             =   I_water.*evap_re_water + I_land.*evap_re_land;

% Evaporative fraction
I_1                                 =   ((Rn - G0) ~= 0);
I_2                                 =   ((Rn - G0) == 0);
evap_fr_1                           =   evap_re .* (Rn - G0 - H_WL) ./ (Rn - G0);
evap_fr_2                           =   1;                                                  % for negative available energy
evap_fr                             =   min(I_1.*evap_fr_1 + I_2.*evap_fr_2,1);

%% Latent heat
LE                                  =   evap_fr .* (Rn - G0);                                % Latent heat flux W m-2


%% Potential ET
%% Reference ETo. Formulación original de Xuelong, pero sin dar valores de ri_i

% LE0                                 = ((slope_WL .* re_i .* (Rn - G0) + rhoa_m_Cp .* (esat_WL - Earef)) ...
%                                         ./ (re_i.*(gamma + slope_WL) + gamma * ri_i));
%P-M
ri_i=0;              %al tener suelo y vegetaci?n es muy dif?cil establecer la resistencia de la cubierta, si pongo cero, tendr?a un l?mite te?rico

LEp_PM                                 = ((slope_WL .* re_i .* (Rn - G0) + rhoa_m_Cp .* (esat_WL - Earef_V)) ...
                                         ./ (re_i.*(gamma + slope_WL) + gamma * ri_i));
%P-T                                     
                                     
L_e 	   =    2.501-0.002361*(Tref_K_V-273.15); 	    % latent heat of vaporation (MJ kg-1)                                     


LEp_PT09                                 =(0.9/L_e).*(slope_WL./(gamma + slope_WL)).* (Rn - G0);

LEp_PT126                                 =(1.26/L_e).*(slope_WL./(gamma + slope_WL)).* (Rn - G0);
                                     
                                     
%P-M FAO56 reference ET for grass                                     
LE0_fao56                                 = ((slope_WL .* (208/Uref_V) .* (Rn - G0) + rhoa_m_Cp .* (esat_WL - Earef_V)) ...
                                         ./ ((208/Uref_V).*(gamma + slope_WL) + gamma * 70));
                                     
% ?21/?12/?2016
return


function [psim] = PSIm(zeta)
Y                                   =   -zeta;                                              %Brutsaert 2008, p50
% Integrated stability function
% 1a. Stability correction function for momentum, eq.(16)
% Y=-z/L, z is the height, L the Obukhov length, both specified in the calling statement.
% For stable conditions
% we use the expressions proposed by Beljaars and Holtslag (1991)and evaluated by Van den Hurk and
% Holtslag (1995) can be used.

a_s                                 =   1.0;                                                % constants, p. 122, van der Hurk & Holtslag (1997)
b_s                                 =   0.667;                                              % constants, p. 122, van der Hurk & Holtslag (1997)
c_s                                 =   5.0;                                                % constants, p. 122, van der Hurk & Holtslag (1997)
d_s                                 =   0.35;                                               % constants, p. 122, van der Hurk & Holtslag (1997)% QUESTION: In page 24, d_s=1 

a_u                                 =   0.33;                                               % constants, p. 49, Brutsaert(2008)
b_u                                 =   0.41;                                               % constants, p. 49, Brutsaert(2008)

I_s                                 =   (Y < 0.0);                                          % STABLE conditions (According to Beljaars & Holtslag, 1991, eq. 13)
I_u                                 =   (Y >= 0.0);                                         % UNSTABLE conditions(% According to Brutsaert 2008, p50)
I_u1                                =   (Y <= b_u^(-3));
I_u2                                =   (Y >  b_u^(-3));

y_s                                 =   -Y;                                                 %due to formulation of Beljaars and Holtslag 1991
y_u1                                =   Y;
y_u2                                =   b_u^(-3);

x_u1                                =   (y_u1/a_u).^(1/3);    
x_u2                                =   (y_u2/a_u).^(1/3);        

y_u                                 =   I_u1.*y_u1 + I_u2.*y_u2;
x_u                                 =   I_u1.*x_u1 + I_u2.*x_u2;
    
PSI0                                =   -log(a_u) + sqrt(3)*b_u*(a_u^(1/3))*pi/6; 

psim_s                              =   -(a_s * y_s + b_s*(y_s - c_s/d_s).*exp(-d_s*y_s)+b_s*c_s/d_s); 
psim_u                              =   log(a_u+y_u) - 3*b_u*(y_u.^(1/3)) +  b_u*a_u^(1/3)/2 * log((1+x_u).^2 ./ (1-x_u+x_u.^2)) +...
                                            sqrt(3) * b_u*a_u^(1/3) *   atan((2*x_u - 1)/sqrt(3)) + PSI0;
psim                                =   I_s .* psim_s + I_u .* psim_u;
return

function [psih] = PSIh(zeta) 
Y                                   =   -zeta;                                                %Brutsaert 2008, p50
% Integrated stability function 1b.
% Stability correction function for heat, eq.(17)
% Y=-z/L, z is the height, L the Obukhov length, both specified in the calling statement.
% For stable conditions
% we use the expressions proposed by Beljaars and Holtslag (1991)
% and evaluated by Van den Hurk and Holtslag (1995) can be used.

a_s                                 =   1.0;                                                % constants, p. 122, van der Hurk & Holtslag (1995)
b_s                                 =   0.667;                                              % constants, p. 122, van der Hurk & Holtslag (1995)
c_s                                 =   5.0;                                                % constants, p. 122, van der Hurk & Holtslag (1995)
d_s                                 =   0.35;                                               % constants, p. 122, van der Hurk & Holtslag (1995)% QUESTION: d_s=1 in page 34 of SEBS document of Su

c_u                                 =   0.33;                                               % constants, p. 443, Brutsaert, 2008
d_u                                 =   0.057;                                              % constants, p. 443, Brutsaert, 2008
n                                   =   0.78;                                               % constants, p. 443, Brutsaert, 2008

y_s                                 =   -Y;                                                 %due to formulation of Beljaars and Holtslag 1991
y_u                                 =   Y;

I_s                                 =   (Y < 0);
I_u                                 =   (Y >= 0);

psih_s                              =   -((1 + 2*a_s/3 * y_s).^1.5                      ... % STABLE conditions (According to Beljaars & Holtslag, 1991 eq. 13    )    
                                        +b_s * (y_s - c_s / d_s).*exp(-d_s * y_s)       ...
                                        +b_s * c_s / d_s - 1);
psih_u                              =   ((1 - d_u) / n) * log((c_u + (y_u.^ n)) / c_u);     % UNSTABLE conditions (According to Brutsaert 2008, p50    
psih                                =   I_s.*(psih_s) + I_u.*(psih_u);

return

function [bw] = Bw(hpbl, L, z0m, d0) 
% NOTES:
% Bulk Stability function for momentum, eq.(22), (26)
% PBL: Height of ABL or PBL
% L: The Obukhov length
% z0m: Surface roughnes height for momentum
% 
% The ASL height
% hst = alfa*PBL, with alfa=0.10~0.15,
% over moderately rough surfaces, or
% hst=beta*z0m, with beta= 100~150.
% 
% Typical values:
% The equations describe the Free convective conditions in the mixed layer,
% provided the top of the ABL satisfies the condition -hst > 2L.
%  
% For a surface with moderate roughness and PBL=1000m, alfa=0.12, this
% is -PBL/L > 17 and -L <60 (and similar values over very rough terrain).
% NOTE: The minus (-) sign infront of B1, B11, B22 are necessary,though not
% clearly specified by Brutsaert, 1999. This is consistent with the integral
% form given in (16)&(17) in which y, the variable is defined as -z/L.
%   (z0m LT (alfa/beta)*PBL): Bw = -ALOG(alfa) + PSIm(alfa*PBL/L) - PSIm(z0m/L) ;(22)
%   (z0m GE (alfa/beta)*PBL): Bw = ALOG(PBL/(beta*z0m)) + PSIm(beta*z0m/L)- PSIm(z0m/L) ;(26)
% B0 = (alfa / beta) * hpbl;
% B1 = -z0m / L;
% B11 = -alfa * hpbl / L;
% B21 = PBL / (beta * z0m);
% B22 = -beta * z0m / L;

alfa                                =   0.12;                                               % These are mid values as given by Brutsaert, 1999
beta                                =   125;                                                % These are mid values as given by Brutsaert, 1999

I_s                                 =   (-z0m ./ L) < 0;                                     % STABLE conditions (Brutsaert, 1982, Eq. 4.93, p.84)
bw_s                                =   -2.2 * log(1 + (-(z0m./L)));

I_u                                 =   (-z0m./L) >= 0;                                     % UNSTABLE conditions (Brutsaert 2008, p53)
I_mr                                =   (z0m<((alfa / beta) * hpbl));                       % for moderately rough terrain (Brutsaert 2008, p53)
I_vr                                =   (z0m>=((alfa / beta) * hpbl));                      % for moderately rough terrain (Brutsaert 2008, p53)

bw_umr                              =   PSIm(alfa*(hpbl-d0)./L) - PSIm(z0m./L) - log(alfa); 
bw_uvr                              =   PSIm(beta *(z0m    )./L) - PSIm(z0m./L) + log((hpbl-d0)./(beta* z0m));
bw_u                                =   max(I_mr.*bw_umr + I_vr.*bw_uvr,0);    

bw                                  =   I_s.*bw_s + I_u.*bw_u;

return

function [cw] = Cw(hpbl, L, z0m, z0h, d0) 
% Bulk Stability function for heat tranfer, eq.(23), (27)
% hpbl: Height of ABL or PBL
% L: The Obukhov length
% z0: Surface roughnes height for momentum
% z0h: Surface roughnes height for height transfer
% 
% The ASL height
% hst = alfa*hpbl, with alfa=0.10~0.15 over moderately rough surfaces, or hst=beta*z0, with beta= 100~150.

alfa                                =   0.12;                                               % These are mid values as given by Brutsaert, 1999
beta                                =   125 ;                                               % These are mid values as given by Brutsaert, 1999

I_s                                 =   ((-z0h ./ L) < 0);                                   % STABLE conditions (Brutsaert, 1982, Eq. 4.93,p.84)
cw_s                                =   -7.6 * log(1 + (-(-z0h ./ L)));      

I_u                                 =   ((-z0h ./ L) >=0);                                   % UNSTABLE conditions (Brutsaert 2008, p53)       
I_umr                               =   (z0m < ((alfa / beta) * hpbl));                     % for moderately rough terrain (Brutsaert 2008, p53)       
I_uvr                               =   (z0m >=((alfa / beta) * hpbl));                     % for very rough terrain (Brutsaert 2008, p53)
cw_umr                              =   PSIh(alfa * (hpbl-d0)./L) - PSIh(z0h./L) - log(alfa);
cw_uvr                              =   PSIh(beta * (z0m    )./L) - PSIh(z0h./L) + log((hpbl-d0)./(beta*z0m));
cw_u                                =   max(I_umr.*cw_umr + I_uvr.*cw_uvr ,0);


cw                                  =   I_s.*cw_s + I_u.*cw_u;

return



