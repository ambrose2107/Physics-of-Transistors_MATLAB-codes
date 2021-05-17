clear all;
% close all
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
Na=3.36e23;
faif = fait *log(Na/ni)
cox= eox/tox
Esi_offset = 3.1
Eg = 1.12
faim = 3.1   % metal workfunction
mueff = 200*1e-4
Vfb = faim - (Esi_offset + Eg/2+faif)
debye_lenght = sqrt(esi*fait/(q*Na))

%% Pao-sha IDVG implement
 
Vgmax = 5
VDS = [ 0.1  1.5 3 5 ]

dv = 0.01
delta_sai= 0.01


VG = 0:0.25:Vgmax
for i=1:length(VDS)
     vgg =  0 
  Ids = []
   for Vg=0:0.25:Vgmax
%  for Vg=2

    IDV = []
    SSAI = []
     vgg = vgg + 1
      Sai = -Vg-VDS(i)-abs(Vfb) :10e-3:Vgmax+VDS(i)
     
  for V=0:dv:VDS(i)
      Fun21=@(sai) ni^2/Na*(exp((sai )/fait))
Fun1=@(sai) ni^2/Na*(exp((sai - V)/fait))
Fun2=@(sai) sqrt(2 *k*T*Na/esi).*sqrt(q*sai/(k*T)+Fun1(sai)/Na) ;

f1_by_f2 = @(sai) Fun1(sai)./Fun2(sai)

% fun3 = @(sai)  Vfb+ sai + (sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)  +Fun1(sai)/Na)          
fun3 = @(sai)  Vfb+ sai + esi/Cox*Fun2(sai)
VGS = fun3(Sai)
sai_s= interp1(real(VGS), real(Sai), Vg)

SSAI = [SSAI sai_s]
intergal_sai = integral(f1_by_f2,delta_sai, sai_s)
 
IDV = [IDV q*mueff*W/L*intergal_sai ]
 end
 
ID= sum(IDV)*dv
 Ids = [Ids ID]
 end
figure(11)
semilogy(VG, Ids*1e6 ,'--'   , "linewidth",3)
hold on 
legend("IDVG @ Vd = 0.1 Pao-sah ","IDVG @ Vd = 1.5 Pao-sah","IDVG @ Vd = 3 Pao-sah","IDVG @ Vd = 5 Pao-sah")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Pao-sah implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

figure(211)
plot(VG, Ids*1e6 ,'--'  , "linewidth",3)
hold on 
legend("IDVG @ Vd = 0.1 Pao-sah ","IDVG @ Vd = 1.5 Pao-sah","IDVG @ Vd = 3 Pao-sah","IDVG @ Vd = 5 Pao-sah")
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
 
 
%% Pao-sha IDVD implement  
 
 
Vdmax = 5
VGS1 = [ 2 3 4 5]
dv = 0.01
delta_sai= 0.01
VD = 0:dv:Vdmax

for i=1:length(VGS1)
% for i=1 
    Ids = []
    IDV = []
    SSAI = []
     
%      for V=0.0 
      for V=0.0 :dv:Vdmax

       Sai = -VGS1(i)-V -abs(Vfb):0.01:Vdmax+VGS1(i)

           for Vg=VGS1(i)
            Fun1=@(sai) ni^2/Na*(exp((sai - V)/fait))
            Fun2=@(sai) sqrt(2 *k*T*Na/esi).*sqrt(q*sai/(k*T)+Fun1(sai)/Na) ;

            f1_by_f2 = @(sai) Fun1(sai)./Fun2(sai)

            fun3 = @(sai)  Vfb+ sai + (sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)  +Fun1(sai)/Na)          

            VGS = fun3(Sai)
            sai_s= interp1(real(VGS), real(Sai), Vg)

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
legend("IDVD @ VG = 2 Pao","IDVD @ VG = 3 Pao","IDVD @ VG = 4 Pao")
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

%% ---------------  Brews_____________________

clear all;
% close all
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
mueff = 200 *1e-4
Vfb = faim - (Esi_offset + Eg/2+faif)
debye_lenght = sqrt(esi*fait/(q*Na))

%%   Brews -IDVG --> DIrect VD
dv = 0.01
VDS = [ 0.1 1.5 3 5 ]
Vgmax = 5
 
for i=1:length(VDS)
delta_sai= 0.001

Ids = []
SSAIo = []
Ids_bre = []
VGS = 0:dv:Vgmax
SSAIl = []
for Vg=0:dv:Vgmax
%      for V=0:dv:Vds 
IDV = []
  for V = VDS(i) 

    Fun1so =@(sai) ni^2/Na*(exp((sai)/fait))

    Fun1=@(sai) ni^2/Na*(exp((sai - V)/fait))
    Fun2=@(sai) sqrt(2 *k*T*Na/esi).*sqrt(q*sai/(k*T)+Fun1(sai)/Na) ;

    fun_so = @(sai) Vg -(  Vfb+ sai + ((sqrt(2*esi*k*T*Na))/Cox)*sqrt(q*sai/(k*T)  +Fun1so(sai)/Na) )       
    fun_sl = @(sai) Vg -(  Vfb+ sai  + (sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)  +Fun1(sai)/Na) )         
    % VGS = fun3(Sai)
    % sai_s= interp1(VGS, Sai, Vg)
    sai_so= fsolve(fun_so,Vg); 
    sai_sl= fsolve(fun_sl,Vg); 
    SSAIo = [SSAIo sai_so]
   SSAIl = [SSAIl sai_sl]
     Id = @(sai)  Cox*(Vg-Vfb- sai) - sqrt(2*esi*q*Na*sai) + 2*fait*(Cox^2*(Vg-Vfb- sai)+esi*q*Na)/(Cox*(Vg-Vfb- sai) + sqrt(2*esi*q*Na*sai))     
     intergal_sai = integral(Id,sai_so, sai_sl)

    IDV = [IDV  mueff*W/L*intergal_sai ]
  end
 
ID = sum(IDV)  ;
Ids = [Ids ID];

end
 
figure(11)
semilogy(VGS,  Ids*1e6 ,'-. '   , "linewidth",3)
hold on 
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
legend("IDVG @ Vd = 0.1 Brew","IDVG @ Vd = 1.5 Brew","IDVG @ Vd = 3 Brew","IDVG @ Vd = 5 Brew")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

figure(211 )
plot(VGS, Ids*1e6,'-. '     , "linewidth",3)
hold on 
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
legend("IDVG @ Vd = 0.1 Brew","IDVG @ Vd = 1.5 Brew","IDVG @ Vd = 3 Brew","IDVG @ Vd = 5 Brew")
xlabel({'VG (V)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'Id (uA/um)'},'FontSize',20.5,'FontWeight','bold');
title('Id (uA/um) vs Vg (V): Brews implementation ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

 
end
 
%%   Brews -IDVD --> DIrect VG
dv = 0.1
VGS = [ 2 3 4 5 ]
Vdmax =  5
 
% for i=1:length(VGS)
delta_sai= 0.01

Ids = []
SSAIo = []
Ids_bre = []
VDS = 0:dv:Vdmax

SSAIl = []
for p = 1:length(VGS) 
    IDV = []
   for V=0:dv:Vdmax
    
        Vg = VGS(p) 
%      for V=0:dv:Vds 

    Fun1so =@(sai) ni^2/Na*(exp((sai)/fait))

    Fun1=@(sai) ni^2/Na*(exp((sai - V)/fait))
    Fun2=@(sai) sqrt(2 *k*T*Na/esi).*sqrt(q*sai/(k*T)+Fun1(sai)/Na) ;

    fun_so = @(sai) Vg -(  Vfb+ sai + ((sqrt(2*esi*k*T*Na))/Cox)*sqrt(q*sai/(k*T)  +Fun1so(sai)/Na)   )       
    fun_sl = @(sai) Vg -(  Vfb+ sai + (sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)  +Fun1(sai)/Na) )         
    % VGS = fun3(Sai)
    % sai_s= interp1(VGS, Sai, Vg)
    sai_so= fsolve(fun_so,Vg); 
    sai_sl= fsolve(fun_sl,Vg); 
    SSAIo = [SSAIo sai_so]
    SSAIl = [SSAIl sai_sl]
     Id = @(sai)  Cox*(Vg-Vfb- sai) - sqrt(2*esi*q*Na*sai) + 2*fait*(Cox^2*(Vg-Vfb- sai)+esi*q*Na)/(Cox*(Vg-Vfb- sai) + sqrt(2*esi*q*Na*sai))     
     
     intergal_sai = integral(Id,sai_so, sai_sl)

    IDV = [IDV  mueff*W/L*intergal_sai ]
%   end
 
    ID =  (IDV) ;
 
   end
 
figure(2)
plot(VDS, IDV*1e6 ,'-'  , "linewidth",3)
hold on 
legend("IDVD @ VG = 2 Brew","IDVD @ VG = 3 Brew ","IDVD @ VG = 4 Brew","IDVD @ VG = 5 Brew")
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


%%   ------------ Piece wise - --------

 clear all;
% close all
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
 tox=10e-9
k=1.3806e-23;
T = 300
Cox = eox/tox
fait = k*T/q
Na=3.335e23;  
W=1e-6;  
L=1e-6; % in cm
ni = 1.5e16
faif = fait *log(Na/ni)
cox= eox/tox
Esi_offset = 3.1
Eg = 1.12
faim = 3.1   % metal workfunction
mueff = 200 *1e-4
Vfb = faim - (Esi_offset + Eg/2+ faif)
debye_lenght = sqrt(esi*fait/(q*Na))

%% iDVG  piece wise
sai =2*faif  + 0.5*fait   
sai1=  ni^2/Na*(exp((sai )/fait))
VT =   Vfb+ sai + (sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)+sai1/Na)          
% VT = 0.8
 SAII = []
C_d = sqrt(q*Na*esi/(2*sai));
m = 1+C_d/Cox
 
mu = mueff
IDL=  []
Vgmax = 5
VGs = 0:0.1:Vgmax

syms SAI

VG  = 0.0:0.1:Vgmax
ID =[]
SAII =[]
vg= []
x= 0     
Sai = []
VGGG  =  []
Sai1 = []
mu_eff = mu
Vth= VT
VDS1  = [0.1 1.5 3 5]
for i=1:length(VDS1)
      VDS = VDS1(i)
      ID = []
 for i=1:length(VG)
    VGS =  VG(i)
   if VDS<(VGS-Vth)/m    && VGS>Vth  +fait
    %I_piecewise(j3,j1) = mu_eff*C_ox/L*((VGS-Vfb-2*phi_b-VDS/2)*VDS - 2/3*sqrt(2*eps_sub*q*p_sub)/C_ox*((VDS+2*phi_b)^(1.5)-(2*phi_b)^(1.5)));
    Id = mu*Cox*W/L*((VGS-Vth)*VDS-m*VDS^2/2);
    ID = [ID Id]
   elseif VDS>(VGS-Vth)/m && VGS>Vth   +fait
    Id = mu*Cox*W/L*((VGS-Vth)^2/2/m); 
    ID = [ID Id]
    
   else
    sai = 0.45
    psi_s = 0.45
    kbT = k*T
    eps_sub = esi
    p_sub = Na
    n_sub = ni^2/Na
    D_it  = 0
%   for i1 = 1:1000   
    psi_s =  0.5; %initial value of surface potential
    sai = 0.5
    for i1 = 1:3000  %Newton Rhapson method
    Qs_LF = q*D_it*(-psi_s)-sign(VGS-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(-q*psi_s/kbT) + q*psi_s/kbT -1)+ n_sub*(exp(-q*VDS/kbT)*(exp(q*psi_s/kbT) -1)- q*psi_s/kbT )));
    a = Vfb - VGS - Qs_LF/Cox + psi_s;
    del_a = 1 + q*D_it/Cox + sign(VGS-Vfb)*(1/Cox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(-q*psi_s/kbT) + q*psi_s/kbT -1)+ n_sub*(exp(-q*VDS/kbT)*(exp(q*psi_s/kbT) -1)- q*psi_s/kbT ))))*(p_sub*q/kbT*(1-exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(-q*VDS/kbT)*exp(q*psi_s/kbT)-1));
    del_psi_s = -a/del_a;
    psi_s = psi_s + del_psi_s;
    end
   
    Sai = [Sai psi_s];
      Sai1 = [Sai1 VGS]
   sai = psi_s
%     I_piecewise(j3,j1) = mu_eff/L*sqrt(eps_sub*q*p_sub/2/sai)*(kbT/q)^2*(n_sub/p_sub)*exp(q*sai/kbT)*(1-exp(-q*VDS/kbT));
    Id = mu*W/L*sqrt(esi*q*Na/2/sai)*fait^2 *ni^2/Na^2*(exp((sai )/fait)) *(1-exp((-VDS/fait)));    
    ID = [ID Id] 
   end
 end
figure(211)
semilogy( VG, ID*1e6, '-' , "linewidth",3)
hold on
% legend("IDVG @ Vd = 0.1 piece","IDVG @ Vd = 1.5 piece","IDVG @ Vd = 3 piece","IDVG @ Vd = 5 piece")
legend("IDVG @ Vd = 0.1 Pao-sah ","IDVG @ Vd = 1.5 Pao-sah","IDVG @ Vd = 3 Pao-sah","IDVG @ Vd = 5 Pao-sah","IDVG @ Vd = 0.1 Brew","IDVG @ Vd = 1.5 Brew","IDVG @ Vd = 3 Brew","IDVG @ Vd = 5 Brew","IDVG @ Vd = 0.1 piece","IDVG @ Vd = 1.5 piece","IDVG @ Vd = 3 piece","IDVG @ Vd = 5 piece")

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
plot( VG, ID *1e6, '-' , "linewidth",3)
 hold on
% legend("IDVG @ Vd = 0.1 piece","IDVG @ Vd = 1.5 piece","IDVG @ Vd = 3 piece","IDVG @ Vd = 5 piece","IDVG @ Vd = 0.1 Pao-sah ","IDVG @ Vd = 1.5 Pao-sah","IDVG @ Vd = 3 Pao-sah","IDVG @ Vd = 5 Pao-sah")
legend("IDVG @ Vd = 0.1 Pao-sah ","IDVG @ Vd = 1.5 Pao-sah","IDVG @ Vd = 3 Pao-sah","IDVG @ Vd = 5 Pao-sah","IDVG @ Vd = 0.1 Brew","IDVG @ Vd = 1.5 Brew","IDVG @ Vd = 3 Brew","IDVG @ Vd = 5 Brew","IDVG @ Vd = 0.1 piece","IDVG @ Vd = 1.5 piece","IDVG @ Vd = 3 piece","IDVG @ Vd = 5 piece") 
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
figure(11)
title('Id (uA/um) vs VD (V): Pao-sah, Brews and Piece-wise Comparison ','FontSize',17);
 
figure(211)
title('Id (uA/um) vs VD (V): Pao-sah, Brews and Piece-wise Comparison ','FontSize',17);
 
%% iDVD 
sai =2*faif  + 0.5* fait
sai1=  ni^2/Na*(exp((sai )/fait))
VT =   Vfb+ sai + (sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)+sai1/Na)          

 SAII = []
C_d = sqrt(q*Na*esi/(2*sai));
m = 1+C_d/Cox
mu = mueff
IDL=  []
Vdmax = 7
VD  = 0:0.1:Vdmax
VGS =2
syms SAI

VGG  = [ 2 3 4 5]
ID =[]
SAII =[]
x= 0

for j = 1:length(VGG)
      VGS =  VGG(j)
      ID =[]
for i=1:length(VD)
    VDS =  VD(i)
   
    if VGS> VT && VDS< (VGS-VT)/m 
     Id = mu*W/L*Cox*  (VGS-VT -m*VDS/2)*VDS
     ID = [ID Id]
     
    elseif VGS> VT && VDS> (VGS-VT)/m 
 
     Id = 0.5*mu*W/L*Cox*  (VGS-VT)^2/m    
     ID = [ID Id]
      
    else
     f3 =  @(SAI)  Vfb+ SAI +(sqrt(2*esi*k*T*Na)/Cox).*sqrt(q*sai./(k*T)+ni^2/Na*(exp((SAI-VDS )/fait))/Na)- VGS         
     Sai  =fsolve(f3,VDS)
         SAII = [SAII Sai]
     Id = mu*W/L*sqrt(esi*q*Na/2/Sai)*fait^2 *ni^2/Na^2*(exp((Sai )/fait)) *(1-exp((-VDS/fait)))   
     ID = [ID Id]
     
     
end
end
figure(2 )
plot(VD, ID*1e6,'*-')
hold on
% legend("IDVD @ VG = 2 piece","IDVD @ VG = 3 piece","IDVD @ VG = 4 piece","IDVD @ VG = 5 piece","IDVD @ VG = 2 Pao","IDVD @ VG = 3 Pao","IDVD @ VG = 4 Pao","IDVD @ VG = 2 Brew","IDVD @ VG = 3 Brew ","IDVD @ VG = 4 Brew","IDVD @ VG = 5 Brew")
legend("IDVD @ VG = 2 Pao","IDVD @ VG = 3 Pao","IDVD @ VG = 4 Pao","IDVD @ VG = 5 Pao","IDVD @ VG = 2 Brew","IDVD @ VG = 3 Brew ","IDVD @ VG = 4 Brew","IDVD @ VG = 5 Brew","IDVD @ VG = 2 piece  ","IDVD @ VG = 3 piece  ","IDVD @ VG = 4 piece  ","IDVD @ VG = 5 piece  ", "IDVD @ VG = 2 piece (2*phit)","IDVD @ VG = 3 piece (2*phit)","IDVD @ VG = 4 piece (2*phit)","IDVD @ VG = 5 piece (2*phit)")
xlim([0 5])
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
figure(2)
title('Id (uA/um) vs VD (V): Pao-sah, Brews and Piece-wise Comparison ','FontSize',17);


 
 



