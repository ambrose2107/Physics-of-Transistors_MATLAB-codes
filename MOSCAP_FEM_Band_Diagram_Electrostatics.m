clear all
clc
close all
%% Define dx and length
dx = 5E-8
x  = 0:dx:50E-7
N  = length(x)
VG = -1   % At what value to plot
NSUB  =  1E17; % Choose +vs for PMOS and -ve NMOS
Phi_m = 3.65  % Metal WF
tox   = 5E-7;  %oxide thickness
%% Variables
q = 1.6*10^-19;
epsilon = 8.854E-14;
e_si = 12;
e_ox = 4;
ni = 1E10;
tm = 15E-7  ; % Metal thickness
Eg_si = 1.1 ; %Si Eg
Esi_cboffset = 3.1;
Eg_ox = 9; 
fait = 0.026;
faif =  -1*sign(NSUB)*fait*log(abs(NSUB)/ni)
% Vfb = sign(NSUB)*Eg_si/2  - faif
tol = 1e-9  ;
T = 300     ;
k = 1.38e-23;
fait=k*T/q;
Vfb = Phi_m - (Esi_cboffset + Eg_si/2 + faif)   ;  %ideal W/o charges
% Vfb = 0
%% Inbuilt function for A
A = (-2*diag(ones(1,N),0)) + 1*diag(ones(1,N-1),1)+ ...
     1*diag(ones(1,N-1),-1)
A(1,:) = 0 ; A(1,1) = 1;     A(N,:) = 0 ; A(N,N) = 1
Sai = zeros(N,1); B=  zeros(N,1);
%% Sais vs VG -- DEFINE iterating parameters and Poisson solve
Cox=epsilon *e_ox/tox;
const = ni^2/(NSUB^2);
if sign(NSUB)== -1 % NMMOS
     sai_r = -0.2:0.01:1.2;  % NMOS range from Acc --> Dep --> Inv
elseif sign(NSUB) ==  1 % PMOS
     sai_r = -1.2:0.01:0.2;  % NMOS range from Acc --> Dep --> Inv
end   
% sai_r = -0.8
Vgg = []  % Extract the VG for given Sais
for s= 1:length(sai_r)
    sai = sai_r(s)
% Poissons
if sign(NSUB)== -1 % NMMOS
    F_sai1 = exp(-sai/fait)+(sai/fait)-1
    F_sai2 = const*(exp((sai)/fait) - (sai/fait)- 1)
elseif sign(NSUB) ==  1 % PMOS
    F_sai1 = const*(exp(-sai/fait)+(sai/fait)-1)
    F_sai2 = (exp((sai)/fait) - (sai/fait)- 1)
end
    F_sai = sqrt(F_sai1 + F_sai2);
    Qs = sqrt(2*epsilon*e_si*k*T*abs(NSUB))* F_sai;
if sign(NSUB)== -1 % NMMOS
    if sign(sai) == -1        % NMOS accumulation
       Vg = Vfb + sai - Qs/Cox ;
    else                      % NMOS - Dep + Inv
       Vg = Vfb + sai + Qs/Cox ;
    end
elseif sign(NSUB) ==  1 % PMOS
    if sign(sai)  == -1        % NMOS accumulation
       Vg = Vfb + sai - Qs/Cox ;
    else                      % NMOS --> Dep + Inv
       Vg = Vfb + sai + Qs/Cox ;
    end
end
 Vgg= [Vgg Vg] 
end
figure(5)
plot(Vgg, sai_r,'LineWidth',3)
ylabel({'Sai (V)'},'FontSize',14.5,'FontWeight','bold');
xlabel({ 'VG (V)'},'FontSize',14.5,'FontWeight','bold');
hold on 
set(gca, 'xlim', [-1.5 3])
title('Sai (V)  vs Vg ','FontSize',17);ax = gca;
c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k';ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; ax.Color = 'white'; 
grid on

%% % AT what point to claculate
SAI_surf = interp1(Vgg, sai_r, VG  ,'spline') 
% SAI_surf = 0
%% Choose NMOS or PMOS
if sign(NSUB)== -1 % NMMOS
    psub0 = abs(NSUB)
    nsub0 = ni^2/abs(NSUB)
elseif sign(NSUB)== 1 % PMOS
    nsub0 = abs(NSUB)
    psub0 = ni^2/abs(NSUB)
    else   % Default NMOS
           nsub0 = ni;
           psub0 = ni;
end
%% NR loop
for i=1:150
 p = psub0 * exp((-1*Sai)/fait)
 n = nsub0* exp(Sai/fait)
 rho =(p - n + (NSUB))
 B =  -q * dx^2*(rho)/(epsilon*e_si);  
 B(1,1) = SAI_surf; B(N,1) = 0;
 f= A*Sai-B;
% Jacobian
 delp= -1*p/fait
 deln= 1*n/fait
 delrho = (delp - deln)
 delb = -q *dx^2.*delrho/(epsilon*e_si);  
 delb(1) = 0; delb(N) = 0; % Since constants
 J = A-diag(delb)
 dv = -J\f
 if max(abs(dv))< tol
     break
 end
 Sai = Sai+dv
end
%% Electric field in Silicon and Oxide- Ex-V
E_Si = -gradient(Sai,dx)
E_S = (Sai(1)-Sai(2))/dx
E_ox=  e_si* E_S/e_ox
Ef_metal_Pos =   -VG
xox =  -tox:dx:0
 Sai_ox = Sai(1)- E_ox*xox
  Vox = E_ox*tox   % Direct
 VG = Vox + Vfb +Sai(1)
 Vox = VG - Vfb -Sai(1)
E_oxx = Vox/tox ; 
Sai_oxx= Sai(1) - E_oxx*xox

% figure(6)
% plot([xox x], [Sai_oxx' ;Sai ],[-tox -tox-tm],[VG-Vfb ; VG-Vfb],'LineWidth',3 ) % '--','LineWidth',3
% hold on
% ylabel({'Sai (V)'},'FontSize',20.5,'FontWeight','bold');
% xlabel({ 'x (cm)'},'FontSize',20.5,'FontWeight','bold');
% 
% figure(6)
% plot([xox x], [Sai_ox' ;Sai ],[-tox -tox-tm],[VG ; VG],'LineWidth',3 ) % '--','LineWidth',3
% hold on
% ylabel({'Sai (V)'},'FontSize',20.5,'FontWeight','bold');
% xlabel({ 'x (cm)'},'FontSize',20.5,'FontWeight','bold');
% 
% % ELetric field
% figure(9)
% 
% plot([xox(1:end-1) x(1:end -1)],[-diff(Sai_oxx)'/dx;-diff(Sai)/dx]  , 'LineWidth',3) 
% title('Electric field vs x ','FontSize',17);
%% Carrier conc. as f(x)
p = psub0* exp(-1*Sai/fait)
n = nsub0* exp(Sai/fait)
%% Band diagram
xm = -tm:dx:-tox
Ei = -Sai
EC = -Sai+Eg_si/2 
EV = -Sai-Eg_si/2
Ecox = -1*(Sai_ox)+Esi_cboffset + Eg_si/2  %Viv

E_ox_R=  EC(1) + Esi_cboffset
E_ox_L = E_ox_R - sign(Vox)*Vox
E_ox_CB = E_ox_R - E_oxx* abs(xox)
E_ox_VB = E_ox_CB - Eg_ox
% Ecox = -1*(Sai_ox)+Esi_cboffset + Eg_si/2
% Evox = Ecox - Eg_ox

Ef = (Ei(end)- faif)*ones(size(Ei))
Em = (Ef(end)-0)*ones(size(xm))-VG
Em1 = (E_ox_CB(1)-Phi_m)*ones(size(xm))-VG

figure(88)
 plot([xm xox x] ,[Em1'; E_ox_CB' ;EC],'LineWidth',3) % Oxide left side to metal
hold on
plot([xm xox x] ,[Em1'; E_ox_VB' ;EV],'LineWidth',3) % metla WF
plot([x] ,[Ef],'--','LineWidth',3,'color','k') % metla WF
plot([x] ,[Ei],'-.','LineWidth',3) % metla WF
xlabel({'x (cm)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'E (eV)'},'FontSize',20.5,'FontWeight','bold');
title('E (eV) vs x ','FontSize',17);ax = gca;
c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k';ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; ax.Color = 'white'; 
grid on
set(gca, 'xlim', [-2E-6  3E-6])

%%
figure(55)
subplot(2,2,1)
semilogy(x,p,x,n,'LineWidth',3)
hold on
semilogy([0 x(end)],[ni ni],'--','LineWidth',3 )
xlabel({'x (cm)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'n & p (cm^{-3})'},'FontSize',20.5,'FontWeight','bold');
title('n & p (cm^{-3}) vs x ','FontSize',17);
legend('p','n','ni  (cm^{-3})','FontSize',17)
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k'; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 1; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; ax.Color = 'white'; 
grid on
% Band diagram
subplot(2,2,2)
plot([xm xox x] ,[Em'; E_ox_CB' ;EC],'LineWidth',3) % Oxide left side to metal
hold on
plot([xm xox x] ,[Em'; E_ox_VB' ;EV],'LineWidth',3) % metla WF
plot([x] ,[Ef],'--','LineWidth',3,'color','k') % metla WF
plot([x] ,[Ei],'LineWidth',3) % metla WF
xlabel({'x (cm)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'E (eV)'},'FontSize',20.5,'FontWeight','bold');
title('E (eV) vs x ','FontSize',17);ax = gca;
c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k';ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; ax.Color = 'white'; 
grid on
% set(gca, 'xlim', [-2E-6  5E-6])


% potential
subplot(2,2,3)
hold on
plot([xox x], [Sai_oxx' ;Sai ],[-tox -tox-tm],[VG - Vfb; VG-Vfb],'LineWidth',3 ) % '--','LineWidth',3
xlabel({'x (cm)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'V (V)'},'FontSize',20.5,'FontWeight','bold');
title('V (V) vs x','FontSize',17);
% legend('n & p (cm^{-3})','FontSize',17)
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

% Electric field
subplot(2,2,4)
hold on
plot([xox(1:end-1) x(1:end -1)],[-diff(Sai_oxx)'/dx;-diff(Sai)/dx]  , 'LineWidth',3) 
xlabel({'x (cm)'},'FontSize',20.5,'FontWeight','bold');
ylabel({'E (V/cm)'},'FontSize',20.5,'FontWeight','bold');
title('E (V/cm) vs x ','FontSize',17);
ax = gca; c = ax.Color; ax.Color = 'white'; % This sets background color to black.
ax.Color = 'k' ; ax.YColor = 'r'; darkGreen = [0, 0.6, 0]; % Make the x axis dark green.
ax.XColor = darkGreen; ax.GridColor = 'k'; % Make the grid color yellow.
ax.GridAlpha = 0.8; % Set's transparency of the grid.
ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.FontWeight = 'bold'; 
ax.Color = 'white'; grid on

