%%
% Définition des paramètres du modèle 
a1 = -2;
a2 = 1;
b1 = (kb*k1*Te*Te)/(2*k2);
b2 = (kb*k1*Te*Te)/(2*k2);
%% 
% Définition des performances désirées 

p1 = -2*exp(-Te/taudes);
p2 = exp(-2*Te/taudes);
p3 = 0;
p4 = 0;
%%
% Résolution manuelle de l'équation de Bezout
ap1 = a1-1;
ap2 = a2-a1;
ap3 = -a2;

Pdes = [1 p1 p2 p3 p4];
x= inv([1 0 0 0 0; ap1 1 b1 0 0; ap2 ap1 b2 b1 0;ap3 ap2 0 b2 b1;...
   0 ap3 0 0 b2])*Pdes';

% Paramètres du correcteur
sp1 = x(2); % il s'agit de s'1 du cours.
r0  = x(3);
r1  = x(4);
r2  = x(5);