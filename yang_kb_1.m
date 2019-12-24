function [z0h]=  yang_kb_1(z0m,zm,zh,wspd,tsfc,tair,qair,psfc)
% zm=Zref;zh=Zref;wspd=Uref;tsfc=LST_K;tair=Tref_K;qair=qaref;psfc=Pref;
% OUT   strout    structure with output data, with
%		     z0h      heat roughness lenght (m)    
%
% IN    indata     structure with input data, at least with
%		     z0m      momentum roughness lenght (m)  
%            zm       air temperature observation height (m)
%		     zh       wind speed observation height (m)
%		     wspd     wind speed (m/s)
%		     tsfc     land surface temperature (K)
%		     tair     air temperature (K)
%            qair     specific humidity (kg/kg) 
%            psfc     pressure (pa)
%#######################################################################
kv     =  0.4;     % Von Karman constant
rd     =  287.0;   % Gas constant for dry air  (m**2/(s**2*K))
cp     =  1004.0;  % Specific heat of dry air at constant pressure(m^2/(s^2*K)).
rddcp  =  rd/cp;
g      =  9.8;     % Acceleration due to gravity at the earth surface.% (m/(s**2)) 
p0     =  1.0e5;   % Surface reference pressure, is 100000 Pascal. 
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%     Beginning of executable code...
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rhoair =  psfc./(rd*tair.*(1+0.61*qair));
ptair  =  tair .* (psfc./(psfc-rhoair*9.81.*zh)).^rddcp;
ptsfc  =  tsfc;       
pt1    =  ptair;      
c_u    =  kv ./log(zm./z0m);   %c_u 初始值
c_pt   =  kv ./log(zh./z0m);   %c_pt 初始值
tstr   =  c_pt.*(pt1 - ptsfc); %tstr 初始值
ustr   =  c_u .* wspd; 
% figure;plot(tstr,'.');
lmo    =  ptair.*ustr.^2./(kv*g*tstr);   %  M-O length      
nu     =  1.328e-5*(p0./psfc).* (pt1/273.15).^1.754;  % viscosity of air  
for i = 1:3   %此处循环用于收敛 z0h,而每次循环的的z0h值与，ustr,tstr有关，因此每次值与函数flxpar相关
    %use yang's model (2002,QJRMS)
    %CALL z0mz0h(zh,z0m,nu,ustr,tstr,z0h);
    [z0h]      =   z0mz0h(zh,nu,ustr,tstr);
    %CALL flxpar(zm,zh,z0m,z0h,wspd,ptsfc,pt1,lmo,c_u, c_pt);
    [c_u,c_pt] =   flxpar(zm,zh,z0m,z0h,wspd,ptsfc,pt1);  
    ustr       =   c_u.*wspd;  %ustr每次计算值与fortran都不同
    tstr       =   c_pt.*(pt1-ptsfc);
end            