
 
%Load input matrices:
%load('LST_V.mat');
%load('SWd_V.mat');
%load('LWd_V.mat');
%load('Pref_V.mat');
%load('Tref_K_V.mat');
%load('Earef_V.mat');
%load('Uref_V.mat');
%load('fc_V.mat');
%load('hc_V.mat');
%load('LAI_V.mat');
%load('albedo_V.mat');
%load('emissivity_V.mat');  
hpbl=1500; 
%??Rhref1 = convert_humidity (Pref, Tref_K, qaref, 'specific humidity', 'relative humidity','Bolton1980');

    
    %%Constants
    constants

    %%--------emissivity

    
    % preliminary processing
    %e0=VAPPREMATRIX(Tref_K-273.15);   %T ???????????? 

    LST_V(LST_V<240)=NaN;
    LWu=emissivity_V .* Sigma_SB .* LST_V.^4;
    
    Rn=SWd_V.*(1-albedo_V)+LWd_V-LWu;
    Lambda_s                        =	0.17;                                              % MPG: modified from0.315 to 0.17 bare soil (Kustas et al., 1989)
    Lambda_c                        =	0.05;                                              % full vegetation canopy (Monteith, 1973)
    G0                              =   Rn .* (Lambda_c + (1 - fc_V) * (Lambda_s - Lambda_c));
    
    %Zref         =    50;%above ground floor
    %******************************** z0m   
 
    
    %******************************** z0m     
    %ff=dir('.\Z0m_2010\GLASS01A01.V03.*.h25v05.2012253.Z_t.tif');
    %B=abs(doy1-doy);[x,index]=sort(B);
    %ind1=[index(1)];    
    %file=['.\Z0m_2010\',ff(ind1).name];
    %z0m= double(imread(file)); 
    %******************************** d 
    %ff=dir('.\d_2010\GLASS01A01.V03.A*.h25v05.2012253.d_t.tif');
    %B=abs(doy1-doy);[x,index]=sort(B);
    %ind1=[index(1)];    
    %file=['.\d_2010\',ff(ind1).name];
    %d0= double(imread(file)); 
    %[m,n]=size(NDVI);     
    
    %[z0m, d0, z0h]                                        =   kb_1(fc,NDVI, LAI,hc, Zref, Uref, Pref, Tref_K,LST_K,qaref,z0m,d0);
     [z0m, d0, z0h]                                        =   kb_1(fc_V, LAI_V, hc_V, Zref_V, Uref_V, Pref_V, Tref_K_V,LST_V,qaref_V);
                                          
    [ustar,H, LE, G0, H_DL, H_WL, H_i, evap_fr,LEp_PM,LEp_PT09,LEp_PT126,LE0_fao56]    =   EnergyBalance_m(d0, z0m, z0h, fc_V,LAI_V, ..., 
                                                                            Rn, LST_V,...
                                                                            hpbl, Zref_V, Tref_K_V, Uref_V, Earef_V,qaref_V, Pref_V, Pref_V,G0);

   
% H(H>1000|H<-300)=NaN;
% H_i(H_i>1000|H_i<-300)=NaN;
% LE(LE>1000|LE<-300)=NaN;
% Rn(Rn>1000|Rn<-300)=NaN;
% G0(G0>1000|G0<-300)=NaN;
      rho_w    =    single(0.998); % density of water [kg/(m2 mm)]                                                                  
      L_e 	   =    2.501-0.002361*(Tref_K_V-273.15); 	    % latent heat of vaporation (MJ kg-1)
      Rn_MJ    =    Rn *24*3600*30/10^6;    % Net Radiation (MJ m-2 month-1)MPG:chaged to 30 days on average. 
      ET_mm_rn =    evap_fr .* max(Rn_MJ,0) ./ (L_e*rho_w); 	% Monthly Evapotranspiration [mm month-1] derived only from Rn
      ET_mm    =    (LE *24*3600*30)./(10^6*L_e*rho_w);         % Monthly Evapotranspiration [mm month-1] derived from Rn and G0
      ESI_PM      =    LE./LEp_PM;
      ESI_PT09    =    LE./LEp_PT09;
      ESI_PT126   =    LE./LEp_PT126;
      ESI_ETo     =    LE./LE0_fao56;

      %      ett=ETdaily(132-3:132+3,239-3:239+3); % flux tower point, 132, 239

%      ETtower(doy) = mean(ett(:));
%-------------------------------------------------------------------------
%   fn=['.\','ET ',num2str(yy),num2str(mm,'%02d'),num2str(dd,'%02d'),'.mat'];
%   save(fn,'ETdaily','-mat'); 
%    end
% end
%end
%%
%save ETtower ETtower;
%%
%load('ET 20100713');
%imagesc(Tref_K);colorbar
%imagesc(ETdaily);colorbar
%ETdaily = double(ETdaily);
%save('try.txt', 'NDVI', '-ASCII');
%save('try.txt', 'ETdaily', '-ASCII');
%fid=fopen('ET 20100713', 'w');
%fwrite(fid, ETdaily);
%fclose(fid);
%imwrite(ETdaily,'ET 20100713.jpg','jpeg' )
%%
%fn=['.\','ET ',num2str(yy),num2str(mm,'%02d'),num2str(dd,'%02d'),'.tif'];
%imwrite(uint8(ETdaily),fn);

