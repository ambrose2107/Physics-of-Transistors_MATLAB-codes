clear all;
close all
%% close all - Pao-sah
h = 6.626e-34;
hcut= h/2/pi  ;
q = 1.6e-19;
mo = 9.1093837015e-31;
eps = 8.854e-12;
EOX = 3.9;
ESI = 11.7;
eox= eps*EOX;
esi = eps* ESI;
Eg_ox= 9;
Esi_offset = 3.1  ;
tox = 10e-9;
k=1.3806e-23;
T = 300
Cox = eox/tox
fait = k*T/q
tox= 10e-9
W=1e-6;  
L=1e-6; % in cm
ni = 1.5e16
Nd=3.36e23;
faif = fait *log(Nd/ni)
cox= eox/tox
Esi_offset = 3.1
Eg =  1.12   % metal workfunction
mueff = 100*1e-4
faim = 3.1
Vfb = faim + Eg  - (Esi_offset + Eg/2 -faif)
debye_lenght = sqrt(esi*fait/(q*Nd))

%% Pao-sha IDVG implement
 
vgg = []
Vgmax = -4
VDS1 = [ -0.1  -1.5 -3 -5]
%  VDS1 = [ -0.1  ]

dv = -0.01
delta_sai= -0.1
psi_s1 = 0
VG = Vgmax:0.25:0
% VG = -2
   SSAI = []
% for i=1 
for i=1:length(VDS1)

   
  Ids = []
   for Vg= Vgmax: 0.25:0
    IDV = []

  for V=0:dv:VDS1(i)
      p_sub = ni^2/Nd;n_sub = Nd;
%       Fun21=@(sai) ni^2/Na*(exp((  sai )/fait))
%       Fun1=@(sai) ni^2/Na*(exp((  sai - V)/fait))
    Fun1=@(sai) ni^2/Nd*(exp(( -sai +V)/fait))

%     Fun2=@(sai) sqrt(2 *k*T*Nd/esi).*sqrt( -q*sai/(k*T)-Fun1(sai)/Nd) ;
    Fun2=@(sai) sqrt(2*k*T/esi)*sqrt(abs(p_sub*(exp(V/fait)*(exp(- sai/fait)-1) + sai/fait)+ n_sub*((exp(sai/fait) -1)- sai/fait )));

    f1_by_f2 = @(sai) Fun1(sai)./Fun2(sai)
 
        VGS =Vg
         psi_s = -0.05;kbT = k*T
        eps_sub = esi;D_it  = 0
    %   for i1 = 1:1000
            for i1 = 1:1000  %Newton Rhapson method
            Vdd = V
            Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(Vdd/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
            a = Vfb - VGS - Qs_LF/Cox + psi_s;
            del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(Vdd/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*Vdd/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
            del_psi_s =  -a/del_a;
            psi_s = psi_s + del_psi_s;
            if abs(psi_s-psi_s1)<0.0001
                break
            else 
                psi_s1 = psi_s
            end
            end
            sai_s  =  psi_s 
      
    intergal_sai = integral(f1_by_f2,delta_sai, sai_s)
    IDV = [IDV q*mueff*W/L*intergal_sai ]
  end
  SSAI = [SSAI sai_s]
   vgg = [vgg Vg]
ID= sum(abs(IDV))*dv
 Ids = [Ids ID]

   end
 
 

figure(11)
semilogy( VG,abs(Ids)*1e6 ,'--'   , "linewidth",3)
hold on 
legend("IDVG @ Vd = -0.1 Pao-sah ","IDVG @ Vd = -1.5 Pao-sah","IDVG @ Vd = -3 Pao-sah","IDVG @ Vd = -5 Pao-sah")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('PMOS Id (uA/um) vs Vg (V): Pao-sah implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

figure(211)
plot(VG, abs(Ids)*1e6 ,'--'  , "linewidth",3)
hold on 
legend("IDVG @ Vd = -0.1 Pao-sah ","IDVG @ Vd = -1.5 Pao-sah","IDVG @ Vd = -3 Pao-sah","IDVG @ Vd = -5 Pao-sah")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Pao-sah implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

end

 
%% Pao-sha IDVD  
 
% Pao-sha IDVG implement
 
Vdmax = -5

VGS1 = [ -2 -3  -4 -5  ]
 
dv = -0.01
delta_sai= -0.01
VD = 0:dv:Vdmax

for i=1:length(VGS1)
% for i=1 
    Ids = []
    IDV = []
    SSAI = []
     
%      for V=0.0 
      for V=0.0 :dv:Vdmax

 
           for Vg=VGS1(i)
             Fun1=@(sai) ni^2/Nd*(exp(( -sai +V)/fait))

    Fun2=@(sai) sqrt(2 *k*T*Nd/esi).*sqrt( -q*sai/(k*T)+Fun1(sai)/Nd) ;

    f1_by_f2 = @(sai) Fun1(sai)./Fun2(sai)

%    fun3 = @(sai)  Vfb+ sai + (sqrt(2*esi*k*T*Nd)/Cox).*sqrt(q*sai./(k*T)  +Fun1(sai)/Nd)          
% 
%             VGS = fun3(Sai)
%             sai_s= interp1(real(VGS), real(Sai), Vg)

             VGS =Vg
        sai = -0.05;psi_s = -0.05;kbT = k*T
        eps_sub = esi; 
        psi_s1 = 0
    %   for i1 = 1:1000
        for i1 = 1:500  %Newton Rhapson method
         Vdd = V
        Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(Vdd/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
        a = Vfb - VGS - Qs_LF/Cox + psi_s;
        del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(Vdd/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*Vdd/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
        del_psi_s =  -a/del_a;
        psi_s = psi_s + del_psi_s;
        if abs(psi_s-psi_s1)<0.001
            break
        else 
                psi_s1 = psi_s
        end
        end
         
        sai_s  =  psi_s 
            
            SSAI = [SSAI sai_s]
            intergal_sai = integral(f1_by_f2,delta_sai, sai_s)
            IDV = [IDV q*mueff*W/L*intergal_sai ]
           end
    ID = sum(IDV) *dv
    Ids = [ Ids ID]
       end
figure(2)
plot(VD(1,1:end), Ids*1e6 ,'--'   , "linewidth",3)
hold on 
legend("IDVD @ VG = -2 Pao","IDVD @ VG = -3 Pao","IDVD @ VG = -4 Pao","IDVD @ VG = -5 Pao")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Pao-sah implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on
end

%% Brews clear all;
clear all;
 
h = 6.626e-34;
hcut= h/2/pi  ;
q = 1.6e-19;
 mo = 9.1093837015e-31;
eps = 8.854e-12;
EOX = 3.9;
ESI = 11.7;
eox= eps*EOX;
esi = eps* ESI;
Eg_ox= 9;
Esi_offset = 3.1  ;
tox = 10e-9;
k=1.3806e-23;
T = 300
Cox = eox/tox
fait = k*T/q
tox= 10e-9
W=1e-6;  
Na=3.335e23;
L=1e-6; % in cm
ni = 1.5e16  
faif = fait *log(Na/ni)
cox= eox/tox

Esi_offset = 3.1
Eg = 1.12
faim = 3.1   % metal workfunction
mueff = 100 *1e-4
Vfb =  + Eg/2 + faif 
debye_lenght = sqrt(esi*fait/(q*Na))

%%   Brews -IDVG --> DIrect VD
dv =  -0.01
Nd = Na
VDS1 = [-0.05 -1.5 -3 -5]
 psi_s1 = 0
Vgmax =  -4
nsub = Nd
psub = ni^2/Nd
VGs = 0:dv:Vgmax

for i=1:length(VDS1)
delta_sai= -0.001

Ids = []
SSAIo = []
Ids_bre = []
SSAIl = []
    for Vg=0:dv:Vgmax
    IDV = []
    for V = VDS1(i) 
    VGS =Vg;
    psi_s = -0.05;kbT = k*T
    eps_sub = esi;p_sub = ni^2/Na;n_sub = Na;D_it  = 0
     for i1 = 1:1000  %Newton Rhapson method
              VDS  = 0
            Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
            a = Vfb - VGS - Qs_LF/Cox + psi_s;
            del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*VDS/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
            del_psi_s =  -a/del_a;
            psi_s = psi_s + del_psi_s;
                if abs(psi_s-psi_s1)<0.0001
                    break
                else 
                    psi_s1 = psi_s
                end
     end
            sai_so =  psi_s 

            for i12 = 1:1000  %Newton Rhapson method
                VDS = V
            Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
            a = Vfb - VGS - Qs_LF/Cox + psi_s;
            del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*VDS/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
            del_psi_s =  -a/del_a;
            psi_s = psi_s + del_psi_s;
            if abs(psi_s-psi_s1)<0.0001
                break
            else 
                psi_s1 = psi_s
            end
            end
            sai_sl =  psi_s 
 
    SSAIo = [SSAIo sai_so]
    SSAIl = [SSAIl sai_sl]
     Id = @(sai)  Cox*(Vg-Vfb-  (sai)) +sqrt(2*esi*q*Na*abs( (sai)))  - 2*fait*(-Cox^2*(Vg-Vfb-  (sai))+esi*q*Na)/(-Cox*(Vg-Vfb- real(sai)) + sqrt(2*esi*q*Na*abs(real(sai))))     
         
     intergal_sai = integral(Id,sai_so, sai_sl)

    IDV = [IDV  mueff*W/L*intergal_sai ]
    end

    ID = sum(IDV)  ;
    Ids = [Ids ID];

    end
 
figure(113)
plot(VGs, SSAIo   ,'-. '   , "linewidth",3)
hold on 
plot(VGs, SSAIl   ,'-. '   , "linewidth",3)

figure(11)
semilogy(VGs,  abs(real(Ids))*1e6 ,'-. '   , "linewidth",3)
hold on 
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
legend("IDVG @ Vd = -0.1 Brew","IDVG @ Vd = -1.5 Brew","IDVG @ Vd = -3 Brew","IDVG @ Vd = -5 Brew")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

figure(211)
plot(VGs, abs(real(Ids))*1e6,'-.', "linewidth",3)
hold on 
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
legend("IDVG @ Vd = -0.1 Brew","IDVG @ Vd = -1.5 Brew","IDVG @ Vd = -3 Brew","IDVG @ Vd = -5 Brew")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

figure(2114)
plot( VGs, ( (Ids))*1e6,'-.', "linewidth",3)
hold on 
 
end
 
%%   Brews -IDVD --> DIrect VG
dv = -0.1
VGS1 = [ -2 -3 -4 -5 ]
Vdmax =  -5
psi_s1 = 0
 
% for i=1:length(VGS)
delta_sai= -0.01

Ids = []
SSAIo = []
Ids_bre = []
VDs = 0:dv:Vdmax

SSAIl = []
for p = 1:length(VGS1) 
    IDV = []
   for V=0:dv:Vdmax
    
        VGS = VGS1(p) 
        Vg  = VGS
%      for V=0:dv:Vds 
            for i1 = 1:1000  %Newton Rhapson method
              VDS  = 0
            Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
            a = Vfb - VGS - Qs_LF/Cox + psi_s;
            del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*VDS/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
            del_psi_s =  -a/del_a;
            psi_s = psi_s + del_psi_s;
            if abs(psi_s-psi_s1)<0.0001
                break
            else 
                psi_s1 = psi_s
            end
            end
            sai_so =  psi_s 

            for i1 = 1:1000  %Newton Rhapson method
                VDS = V
            Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
            a = Vfb - VGS - Qs_LF/Cox + psi_s;
            del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*VDS/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
            del_psi_s =  -a/del_a;
            psi_s = psi_s + del_psi_s;
            if abs(psi_s-psi_s1)<0.0001
                break
            else 
                psi_s1 = psi_s
            end
            end
            sai_sl =  psi_s 
 
    SSAIo = [SSAIo sai_so]
    SSAIl = [SSAIl sai_sl]
%      Id = @(sai)  Cox*(Vg-Vfb- real(sai)) +sqrt(2*esi*q*Na*abs(real(sai)))  + 2*fait*(-Cox^2*(Vg-Vfb- real(sai))+esi*q*Na)/(Cox*(Vg-Vfb- real(sai)) + sqrt(2*esi*q*Na*abs(real(sai))))     
          Id = @(sai)  Cox*(Vg-Vfb-  (sai)) +sqrt(2*esi*q*Na*abs( (sai)))  - 2*fait*(-Cox^2*(Vg-Vfb-  (sai))+esi*q*Na)/(-Cox*(Vg-Vfb- real(sai)) + sqrt(2*esi*q*Na*abs(real(sai))))     
    
     intergal_sai = integral(Id,sai_so, sai_sl)

    IDV = [IDV  mueff*W/L*intergal_sai ]
%   end
 
    ID =  (IDV) ;
 
   end
 
figure(2)
plot(VDs, IDV*1e6 ,'-'  , "linewidth",3)
hold on 
legend("IDVD @ VG = -2 Brew","IDVD @ VG = -3 Brew ","IDVD @ VG = -4 Brew","IDVD @ VG = -5 Brew")
xlabel({'VD (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs VD (V): Brews implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on
 end




  
%% clear all; Piece wise
 
h = 6.626e-34;
hcut= h/2/pi  ;
q = 1.6e-19;
mo = 9.1093837015e-31;
eps = 8.854e-12;
EOX = 3.9;
ESI = 11.9;
eox= eps*EOX;

esi = eps* ESI;
Eg_ox= 9;
Esi_offset = 3.1  ;
tox = 10e-9;
k=1.3806e-23;
T = 300
Cox = eox/tox
fait = k*T/q
tox= 10e-9
W=1e-6;  
L=1e-6; % in cm
ni = 1.5e16
cox= eox/tox
Esi_offset = 3.1
Eg = 1.12
faim = 3.1   % metal workfunction
mueff = 100 *1e-4
mu = mueff
%% iDVG 
 Sai = []
Nd=3.35e23;  
faif = fait *log(Nd/ni)

Vfb = faim + Eg - (Esi_offset + Eg/2 - faif)
debye_lenght = sqrt(esi*fait/(q*Nd))

sai = -2*faif
sai1=  ni^2/Nd*(exp((-sai )/fait))
VT =   Vfb + sai - (sqrt(2*esi*k*T*Nd)/Cox).*sqrt(-q*sai./(k*T)+sai1/Nd)          
% 
Vggg=  []
Vth = VT
 SAII = []
C_d = sqrt(q*Nd*esi/(2*abs(sai)));
m = 1+C_d/Cox
IDL=  []
Vgmax = 4
VGs = 0:-0.01:-Vgmax
 
syms SAI
VGs  = -0.01:-0.001:-Vgmax

VG  = -0.01:-0.01:-Vgmax
ID =[]
SAII =[]
x= 0
Na =Nd
VG1 = []
 ID1 = []
 
  VDS1  = [-0.05 -1.5 -3 -5]
%  VDS1  = [-0.5]

a= 0
b=0
 
for i=1:length(VDS1)
      VDS = VDS1(i)
      ID = []
 for i=1:length(VG)
     VGS =  VG(i)
     mu = mueff 
    
    if   VGS<Vth+m*VDS 

    %I_piecewise(j3,j1) = mu_eff*C_ox/L*((VGS-Vfb-2*phi_b-VDS/2)*VDS - 2/3*sqrt(2*eps_sub*q*p_sub)/C_ox*((VDS+2*phi_b)^(1.5)-(2*phi_b)^(1.5)));
    Id = mu*Cox*W/L*((VGS-Vth)*VDS-m*VDS^2/2);
    ID = [ID Id]
    elseif    (VGS) < (VT)-2*fait  && (VGS)> (VT+m*VDS) 
%    elseif VDS>(VGS-Vth)/m && VGS>Vth   +fait
    Id = mu*Cox*W/L*((VGS-Vth)^2/2/m); 
    ID = [ID Id]
    
    else
    sai = -0.45;psi_s = -0.45;kbT = k*T
    eps_sub = esi;p_sub = ni^2/Na;n_sub = Na;D_it  = 0
%   for i1 = 1:1000   
    for i1 = 1:1000  %Newton Rhapson method
    Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
    a = Vfb - VGS - Qs_LF/Cox + psi_s;
    del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*VDS/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
    del_psi_s =  -a/del_a;
    psi_s = psi_s + del_psi_s;
    end
   Vggg=  [Vggg VGS]
      Sai = [Sai psi_s];
    sai = psi_s
%     I_piecewise(j3,j1) = mu_eff/L*sqrt(eps_sub*q*p_sub/2/sai)*(kbT/q)^2*(n_sub/p_sub)*exp(q*sai/kbT)*(1-exp(-q*VDS/kbT));
    Id = mu*W/L*sqrt(esi*q*Na/2/abs(sai))*fait^2 *ni^2/Na^2*(exp( -(sai )/fait)) *(1-exp((VDS/fait)));    
    ID = [ID Id] 
   end
 end
 
 
ID1 = interp1(VG,ID, VGs) 
figure(211)
semilogy( VGs,  ID1*1e6, '. ' , "linewidth",3)
hold on
% legend("IDVG @ Vd = 0.1 piece","IDVG @ Vd = 1.5 piece","IDVG @ Vd = 3 piece","IDVG @ Vd = 5 piece")
legend("IDVG @ Vd = -0.05 piece","IDVG @ Vd = -1.5 piece","IDVG @ Vd = -3 piece","IDVG @ Vd = -5 piece")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Piece-wise ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

figure(11)
plot( VG, ID *1e6, '. ' , "linewidth",3)
 hold on
legend("IDVG @ Vd = -0.05 piece","IDVG @ Vd = -1.5 piece","IDVG @ Vd = -3 piece","IDVG @ Vd = -5 piece")
 
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Piece-wise ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on
  
end
%%    iDVD 

IDL=  []
Vdmax = -5
VD  = 0:-0.1:Vdmax
VGS =-2
syms SAI

VGG  = [-2 -3 -4 -5]
ID =[]
SAII =[]
x= 0

for j = 1:length(VGG)
      VGS =  VGG(j)
      ID =[]
for i=1:length(VD)
    VDS =  VD(i)
    if   VGS<Vth+m*VDS 

    %I_piecewise(j3,j1) = mu_eff*C_ox/L*((VGS-Vfb-2*phi_b-VDS/2)*VDS - 2/3*sqrt(2*eps_sub*q*p_sub)/C_ox*((VDS+2*phi_b)^(1.5)-(2*phi_b)^(1.5)));
    Id = mu*Cox*W/L*((VGS-Vth)*VDS-m*VDS^2/2);
    ID = [ID Id]
    elseif    (VGS) < (VT)-2*fait  && (VGS)> (VT+m*VDS) 
%    elseif VDS>(VGS-Vth)/m && VGS>Vth   +fait
    Id = mu*Cox*W/L*((VGS-Vth)^2/2/m); 
    ID = [ID Id]
    
    else
    sai = -0.45;psi_s = -0.45;kbT = k*T
    eps_sub = esi;p_sub = ni^2/Na;n_sub = Na;D_it  = 0
%   for i1 = 1:1000   
    for i1 = 1:1000  %Newton Rhapson method
    Qs_LF =  -sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
    a = Vfb - VGS - Qs_LF/Cox + psi_s;
    del_a = 1 -(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(VDS/fait)* (exp(-q*psi_s/kbT)-1) + q*psi_s/kbT )+ n_sub*((exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))) )*(p_sub*q/kbT*(1-exp(q*VDS/kbT)*exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
    del_psi_s =  -a/del_a;
    psi_s = psi_s + del_psi_s;
    end
   Vggg=  [Vggg VGS]
      Sai = [Sai psi_s];
    sai = psi_s
%     I_piecewise(j3,j1) = mu_eff/L*sqrt(eps_sub*q*p_sub/2/sai)*(kbT/q)^2*(n_sub/p_sub)*exp(q*sai/kbT)*(1-exp(-q*VDS/kbT));
    Id = mu*W/L*sqrt(esi*q*Na/2/abs(sai))*fait^2 *ni^2/Na^2*(exp( -(sai )/fait)) *(1-exp((VDS/fait)));    
    ID = [ID Id] 
   end
end
 
 
figure(2)
plot(VD, ID*1e6,'o')
hold on
legend("IDVD @ VG = -2 piece","IDVD @ VG= -3 piece","IDVD @ VG= -4 piece","IDVG @ Vd = -5 piece")
xlabel({'VD (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs VD (V): Piece-wise ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on


end
%%  Legends
figure(11 ) 
legend("IDVG @ Vd = -0.05 Pao-sah ","IDVG @ Vd = -1.5 Pao-sah","IDVG @ Vd = -3 Pao-sah","IDVG @ Vd = -5 Pao-sah","IDVG @ Vd = -0.05 Brew","IDVG @ Vd = -1.5 Brew","IDVG @ Vd = -3 Brew","IDVG @ Vd = -5 Brew","IDVG @ Vd = -0.05 piece","IDVG @ Vd = -1.5 piece","IDVG @ Vd = -3 piece","IDVG @ Vd = -5 piece")

figure(211 ) 
legend("IDVG @ Vd = -0.05 Pao-sah ","IDVG @ Vd = -1.5 Pao-sah","IDVG @ Vd = -3 Pao-sah","IDVG @ Vd = -5 Pao-sah","IDVG @ Vd = -0.05 Brew","IDVG @ Vd = -1.5 Brew","IDVG @ Vd = -3 Brew","IDVG @ Vd = -5 Brew","IDVG @ Vd = -0.05 piece","IDVG @ Vd = -1.5 piece","IDVG @ Vd = -3 piece","IDVG @ Vd = -5 piece")


figure(2 ) 
legend( "IDVD @ VG = -2 Pao","IDVD @ VG = -3 Pao","IDVD @ VG = -4 Pao","IDVD @ VG = -5 Pao","IDVD @ VG = -2 Brew","IDVD @ VG = -3 Brew ","IDVD @ VG = -4 Brew","IDVD @ VG = -5 Brew","IDVD @ VG = -2 piece","IDVD @ VG= -3 piece","IDVD @ VG= -4 piece","IDVG @ Vd = -5 piece")

% legend( "IDVD @ VG = -2 Brew","IDVD @ VG = -3 Brew ","IDVD @ VG = -4 Brew","IDVD @ VG = -5 Brew","IDVD @ VG = -2 piece","IDVD @ VG= -3 piece","IDVD @ VG= -4 piece","IDVG @ Vd = -5 piece","IDVD @ VG = -2 Pao","IDVD @ VG = -3 Pao","IDVD @ VG = -4 Pao","IDVD @ VG = -5 Pao")







 