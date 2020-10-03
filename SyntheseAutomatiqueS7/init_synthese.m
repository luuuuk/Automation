
close all
clc

% Valeurs numériques
% ------------------
J     = 4.6E-5;
f     = 7E-5; 
rm    = 0.33;
lm    = 0.001;
imax  =100;
nm    = 960;
k     = 0.026; % k=Phi (sujet)

k1    = 7;
k2    = 28.65;
k3    = 0.0159;
ki    = 1.6; 
kb    = 6.1;
Te    = 0.2;

taudes = (5/4.75);

% Réglages de la boucle de courant
% --------------------------------
trdesi=1.5*1e-3;
taudesi=trdesi/3;
tii   = 0.5*lm/rm;% 0.003;
kpi   = rm*tii*ki/taudesi;%3.2;

% Réglages de l'asservissement de l'angle T
% -----------------------------------------
tff = 0.025;   % filtrage de la consigne d'angle
xides=0.7;
trdes=0.1;
wndes=0.43*2*pi/trdes; 
kp=nm*J*wndes^2/(k2*k*ki); %27
kv=(2*xides*J*wndes-f)/(kp*ki*k*k3);%kv = 0.0934;

% Correcteur PID filtré
% -----------------------------------------


%% Commande par retour d'état

% Définition répresentation d'état
A = [0 1
    0 0];
B = [0
    kb/k2];
C = [k1 0];
D = 0;

sys = ss(A,B,C,D);
sysd = c2d(sys,Te);
[Ad,Bd,Cd,Dd] = ssdata(sysd)

% Pôles observateur
Po = [0.3 
    0.3];

% Placement des pôles/Calcul Ld, Ae, Be
Ld = acker(Ad',Cd',Po)'

Ae = [Ad(1,1) Ad(1,2) 0
    Ad(2,1) Ad(2,2) 0
    -Cd(1) -Cd(2) 1]
Be = [Bd
    0]

Pcom = [exp(-Te/taudes) exp(-Te/taudes) 0.05*exp(-Te/taudes)]

% Calcul State Feedback
Fe = acker(Ae, Be, Pcom)

dist = 1*pi/180;