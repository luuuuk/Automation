%% Automatique 2 - Taravil de synthèse
%% Commande et observation quadratiques

close all;
clear all;


%% Definition du système

% Definition des constantes
L = 1000;           % Longeur
I = 2.4 * 10^(-3);  % Pente
Cf = 0.001;         % Coefficient de frottement
g = 9.8;            % Accélération de la gravité

% Profil d'équilibre
he = 2;
qe = sqrt((I*he^3)/Cf);

alpha = g*he - (qe^2/he^2);
beta = (2*qe)/he;
gamma = g*I + 2*g*Cf* (qe^2/he^3);
delta = 2*g*Cf* (qe/he^2);
Delta = L/2;

% Definition des matrices du système
A = 1/Delta * [0 0 1/2 0; 
               0 0 0 1; 
               -alpha 0 beta-Delta*delta -beta; 
               gamma -alpha/2 beta/2 -Delta*delta];
           
B = 1/Delta * [0; 
               0; 
               alpha+Delta*gamma; 
               alpha/2];
           
E = 1/Delta * [-1/2; 
               -1; 
               0; 
               -beta/2];
           
C = [0 1 0 0];
D = 0;

% State space
syst = ss(A,[B E],C,D,'Statename',{'h_1','h_2','q_0','q_1'},'Inputname',{'h_0','q_2'},'Outputname',{'h_2'});


%% Q1 - Commandabilité et observabilité du système

rank(ctrb(A,B))
rank(obsv(A,C))

%%
% Le système est commandable puisque le rang de la matrice de
% commandabilité est égal à l'ordre du système (4). De même, le système
% est observable, comme la matrice d'obserabilité est du rang 4, ce qui est
% égal à l'ordre du système. Un approche quadratique est donc faisable.

eig(A)

%%
% De plus, le système est asymptotiquement stable, puisque tout les pôles
% sont à partie réelle négative. Vu que les pôles sont complémentaires, 
% nous attendons que le système donne une réponse à un échelon oscillante
% avec des oscillations décroissantes.

tzero(syst)

%%
% Le système n'a pas des zéros, ce qui assure qu'il n'y a pas de
% compensation de pôles.


%% Q2 - LQ retour d'état

%%
% Un approche LQ est faisable vu que le système est commandable

% Matrices de pondération
Q = C'*C;   

%%
% Ce choix de Q et justifie puisqu'il mets une pondération sur la sortie

R1 = 2;     % Grand R donne une petite sensibilité au bruit, mais lente
R2 = 0.02;  % Petit R donne une grande rapidité, mais sensible au bruit

%%
% La matrice R devoir être définie positive

% Gains de la commande pour R1 = 2
F1 = lqr(A,B,Q,R1)
G1 = inv(C*inv(-A+B*F1)*B)

% Gains de la commande pour R2 = 0.02
F2 = lqr(A,B,Q,R2)
G2 = inv(C*inv(-A+B*F2)*B)


%% Q3 - Comportements de h2 en BF à condition initiale non nulle

% Condition initiale non nulle
x0 = [0; 0; 0.5; 0];

% Pour R1
A_BF1 = A - B*F1;
B_BF1 = B*G1;

systInit1 = ss(A_BF1, B_BF1, C, D);
systInit1u = ss(A_BF1, B_BF1, -F1, G1);

figure();
subplot(2,1,1);
initial(systInit1, x0); grid on;
subplot(2,1,2);
initial(systInit1u, x0); grid on;

% Pour R2
A_BF2 = A - B*F2;
B_BF2 = B*G2;

systInit2 = ss(A_BF2, B_BF2, C, D);
systInit2u = ss(A_BF2, B_BF2, -F2, G2);

figure();
subplot(2,1,1);
initial(systInit2, x0); grid on;
subplot(2,1,2);
initial(systInit2u, x0); grid on;

% Plot des deux commandes
figure();
subplot(2,1,1);
initial(systInit1u, x0); title("Commande pour R1=2"); grid on;
subplot(2,1,2);
initial(systInit2u, x0); title("Commande pour R2=0.02"); grid on;

%%
% La commande la moins coûteuse est la commande avec R1 = 2. Néanmois,
% le temps de réponse du système est beaucoup plus lente que cel du système
% avec la commande qui utilise R2 = 0.02.


%% Q4 - Simulation

% consigne appliqué à t = 0s
r = 0.1;

% perturbation appliqué à t = 8000s
w = 0.1;

%%
% La simulation donne le même temps de réponse qu'on a trouvé avec
% initial() avant. La perturbation en entrée n'est pas rejetée.


%% Q5 - Observateur d'état

% Matrice de pondération W associé à la perturbation
W = E*E';

% Matrice de pondération V associé au bruit
V1 = 1;

%%
% Le choix des matrices de pondération est acceptable, vu que V>>W, 
% l'observateur va être plutôt lent mais par contre pas cher et pas très 
% sensible au bruit.

% Gain de l'observateur
L1 = lqr(A', C', W, V1)'


%% Q6 - Comportements fréquentiels

% Matrice de pondération V associé au bruit
V2 = 0.01;

%%
% V plus petit, alors nous attendons que le sytème soit moins sensible à la
% perturbation mais en même temps plus rapide et donc plus sensible au
% bruit.

% Gain de l'observateur
L2 = lqr(A', C', W, V2)'

% Analyse du comportement fréquentielle pour gain L1
AObs1 = (A - L1*C);
BObs1 = [L1 -E];
CObs1 = eye(4);
DObs1 = 0;

systObs1 = ss(AObs1, BObs1, CObs1, DObs1);


% Analyse du comportement fréquentielle pour gain L2
AObs2 = (A - L2*C);
BObs2 = [L2 -E];
CObs2 = eye(4);
DObs2 = 0;

systObs2 = ss(AObs2, BObs2, CObs2, DObs2);

figure();
bodemag(systObs1, 'r', systObs2, 'b'); legend(); grid on;

%%
% On observe que l'observateur 2 (pour V=0.01) est moins sensible à la perturbation.
% Néanmois, comme attendu, il est plus sensible au bruit puiqu'il est plus
% rapide.


%% Q7 - Comportement observateur conditions initiales differentes du système

% Construction observateur d'état
AObsFin = (A - L2*C);
BObsFin = [B L2];
CObsFin = eye(4);
DObsFin = 0;

systObs2 = ss(AObsFin, BObsFin, CObsFin, DObsFin);



%%
% On observe que l'erreur de l'estimation converge vers 0 en moins de 1300s.
% En même temps on voit bien le bruit.


%% Q8 - CDRE estimé: comportement en réponse à un échelon

% consigne toujours r=0.1 et perturbation toujours w=0.1

%%
% La simulation montre que la réponse à un échelon a un gain statique
% nul et que la perturbation en entrée n'est pas rejetée. Néanmois, onpeut
% observer le bruit de mésure.


%% Q9 - Bruit random comme perturbation

% variance varBruit=0.01, sample time Te=1
%%
% Perturbation comme bruit random: On observe que les perturbations ne sont pas rejetées.

% w=0 et erreur de paramètre de 20%sur gamma
gamma1 = 1.2*gamma;

% Système avec erreur de paramètre
AErrGamma = 1/Delta * [0 0 1/2 0; 
               0 0 0 1; 
               -alpha 0 beta-Delta*delta -beta; 
               gamma1 -alpha/2 beta/2 -Delta*delta];
           
BErrGamma = 1/Delta * [0; 
               0; 
               alpha+Delta*gamma1; 
               alpha/2];

CErrGamma = C;
DErrGamma = D;

%%
% Erreur de paramètre sur gamme: On observe une erreur d'estimation
% importante sur la première variable d'état, ce qui des due à l'influence 
% de gamma sur ce paramètre dans la matrice A. En plus, la commande devient
% très bruitée. Cela vient du fait que gamma acte sur la commande dans la
% matrice B.


%% Q10 - Action intégrale

% Matrices entendues
Ae = [A zeros(4,1); -C 0];
Be = [B; 0];
Qe = [Q zeros(4,1); 0 0 0 0 0.1]; % Poids moins important sur z que sur h2
Re = R1;

Fe = lqr(Ae,Be,Qe,Re)

Ki= Fe(1:4)
gi = Fe(5)

%%
% Simulation comportement en réponse à un échelon et à une perturbation:
% La perturbation à 8000s est rejetée.

%%
% Perturbation de type bruit random number

%%
% Erreur de paramètre sur gamma de 20%