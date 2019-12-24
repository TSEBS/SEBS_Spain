function [z0m, d0, z0h]=  kb_1(fc_V, LAI_V, hc_V, Zref_V, Uref_V, Pref_V, Tref_K_V,LST_V,qaref_V) 
% KB_1 by Massman, 1999, (24) This a surrogate for the full LNF model for describing the combined
% canopy and soil boundary layer effects on kBH_1. 
% Reference: Su (2001), HESS, Chen et al. 2012, JAMC, Chen et al.,2012,HESSD
% Input:
% fc : Fractional canopy cover (-)
% LAI: canopy total leaf area index
% hc: canopy height
% Zref: reference height of air temperature and humidity 
% Uref: wind speed (m/s)
% Pref: pressure,
% Tref_K: air temperature (K)
% LST_K: land surface temperature (K)
% output:
% d0        -> Zero plane displacement height           (m)
% z0m       -> Roughness height for momentum tranfer    (m)
% z0h       -> Roughness height for heat tranfer        (m)
%% Constants
Cd      =   0.2;        % Foliage drag coefficient
Ct      =   0.01;       % Heat transfer coefficient
Pr      =   0.71;       % Prandtl number
k       =   0.4;        % von Karman constant
T0      =   273.15;     % zero Kelvin [C]
P0      =   101325.;    % Standard pressure (Pa)
Lf      =   0.1;        % the size of the leaves
fs                  =   1 - fc_V; 
Nu                  =   1.327e-5 * (P0 ./ Pref_V) .* ((Tref_K_V / T0).^1.81);                   % Kinematic viscosity of air (Massman 1999b) (10 cm^2/s)
%% U*/U(h)
c1                  =   0.320;                                                              % model constants (Massman 1997)
c2                  =   0.264;                                                              % model constants (Massman 1997)
c3                  =   15.1;                                                               % model constants (Massman 1997)
ust2u_h             =   c1 - c2 * exp(-c3 * Cd * LAI_V);                                      % Ratio of ustar and u(h) (Su. 2001 Eq. 8)
Cd_bs               =   2*ust2u_h.^2;                                                       % Bulk surface drag cofficient (Su 2001)

%% within-canopy wind speed profile extinction coefficient
n_ec                =   Cd * LAI_V ./ (Cd_bs);                                                % windspeed profile extinction coefficient (Su. 2002 Eq 7.)
d2h                 =   1 - 1./(2*n_ec) .* (1 - exp(-2 * n_ec));                            % Ratio of displacement height and canopy height (derived from Su 2002, eq 9)

%% Displacement height and Roughness length for momentum
d0                  =   d2h .* hc_V;                                                          % displacement height

I_1                 =   (fc_V <= 0);
I_2                 =   (fc_V >  0);
z0m                 =   zeros(size(fc_V));
z0m(I_1)            =   0.004;  
% roughness height for momentum lower limit
z0m(I_2)            =   hc_V(I_2).*(1 - d2h(I_2)) .* exp(-k * ust2u_h(I_2).^(-1));            % roughness height for momentum (Eq 10 Su. 2001)

%% KB-1 for canopy only
I_1                 =   (n_ec ~= 0);
I_2                 =   (n_ec == 0);

kB1_c1              =   k * Cd ./(4 * Ct .* ust2u_h .* (1 - exp(-n_ec/2)));                 % KB-1 for Full canopy only (Choudhury and Monteith, 1988)
kB1_c2              =   0.0;                                                                % KB-1 for Full canopy only lower limit
kB1_c               =   I_1.*kB1_c1 + I_2.*kB1_c2;                                          % KB-1 for Full canopy only 

I_1                 =   (Nu ~= 0);
I_2                 =   (Nu == 0);

u_h                =    max(0,Uref_V.* log((hc_V-d0)./z0m)./log((Zref_V-d0)./z0m));               % (within canopy) Horizontal windspeed at the top of the canopy
Ustar_m             =   ust2u_h .* u_h;                                                     % friction velocity 
Re_m1               =   Ustar_m .* hc_V ./ Nu;                                                % roughness Reynolds number for 
Re_m2               =   0;                                                                  % roughness Reynolds number for 
Re_m                =   I_1.*Re_m1 + I_2.*Re_m2;                                            % roughness Reynolds number, = hu*/v    

%% KB-1 for mixed canopy and soil 
I_1                 =   (Nu ~= 0);
I_2                 =   (Nu == 0);

u_h                =    max(0,Uref_V.* log((hc_V-d0)./z0m)./log((Zref_V-d0)./z0m));               % (within canopy) Horizontal windspeed at the top of the canopy
Ustar_m             =   ust2u_h .* u_h;                                                     % friction velocity 
Re_m1               =   Ustar_m .* hc_V ./ Nu;                                                % roughness Reynolds number for 
Re_m2               =   0;                                                                  % roughness Reynolds number for 
Re_m                =   I_1.*Re_m1 + I_2.*Re_m2;                                            % roughness Reynolds number, = hu*/v    

Ct_m                =   (Pr^(-2/3))*(Re_m.^(-1/2));
kB1_m               =   (k * ust2u_h) .* (z0m ./ hc_V) ./ Ct_m;                               % KB-1 for Mixed canopy (Choudhury and Monteith, 1988)

%% KB-1 for soil, Changed by Chen et al. 2013 JAMC
hs=zeros(size(fc_V));
hs=hs+0.004;                                                                                % hs(1:length(LST_K),1) =   0.004; momentum roughness parameter (0.009 ~ 0.024)(Su et al., 1997, IJRS, p.2105-2124.)
Ustar_s             =   Uref_V * k ./ log(Zref_V ./ hs);                                        % Friction velocity in case of soil only. (Brutsaert 2008, Eq 2.41 P43 )[m/2]

I_1                 =   (Nu ~= 0 );
I_2                 =   (Nu == 0 );

Re_s1               =   Ustar_s .* hs ./ Nu;                                                % roughness Reynolds number, = hu*/v    
Re_s2               =   0.0;                                                                % roughness Reynolds number, = hu*/v    
Re_s                =   I_1.*Re_s1 + I_2.*Re_s2;                                            % roughness Reynolds number, = hu*/v   

[z0h ] = yang_kb_1(hs,Zref_V,Zref_V,Uref_V,LST_V,Tref_K_V,qaref_V, Pref_V);                            % Yang et al (2002)

kB1_s               =  log(hs./z0h); 
%% KB-1 Combined all, Chen et al. 2012, JAMC, Chen et al. 2012, HESS
kB_1                =   (fc_V.^2)      .* kB1_c    +...                                       % canopy only
                        2*(fc_V.*fs)   .* kB1_m    +...                                       % mixed canopy and soil
                        (fs.^ 2)     .* kB1_s;                                              % soil only 
                   
         
z0h                 =   z0m ./ exp(kB_1);                                                   % roughness height for heat (su 2002 eq 8)                   

z0h = min(Zref_V/10,max(z0h,1.0E-10));

return
