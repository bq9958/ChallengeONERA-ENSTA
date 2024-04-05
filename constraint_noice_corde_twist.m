function [Lp] = constraint_noice_corde_twist(params)

global Nb_disc

c = params(1:Nb_disc);
n = params(2*Nb_disc+1);
twist = params(Nb_disc+1:2*Nb_disc);
B = params(2*Nb_disc+2);
Lp_liste = [];

for THETA = 1:0.05*pi:pi/2
    COOR_POLAR = [0,THETA,1.5];
    Lp = constraint_noice_direction(c,n,twist,B,COOR_POLAR);
    Lp_liste = [Lp_liste,Lp];
end
Lp = max(Lp_liste);