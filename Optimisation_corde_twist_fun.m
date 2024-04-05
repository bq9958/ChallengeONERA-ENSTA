function [avLp,P_finale,T_finale,minparametre] = Optimisation_corde_twist_fun(Nb_pale,n_init)
%clear all
%close all
%clc

global Nb_disc

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex'); 

rtip=0.12;
c_init0 = [0.15 0.16 0.18 0.2 0.21 0.2 0.19 0.19 0.17 0.15 0.13 0.11 0.1 0.07 0.05]*rtip;
c_init = c_init0 + 0.1*rand(1,Nb_disc)*rtip;
twist_init0 = [35 40 38 36 34 32 29 25 23 20 18 16 14 13 12]*pi/180;
twist_init = twist_init0 + 5*rand(1,Nb_disc)*pi/180;

% n_init = 5000; %RPM
% Nb_pale = 3;

nonlin = @constraint_nonlin;

lb = [ones(1,Nb_disc)*0.08*rtip,ones(1,Nb_disc)*pi/10,3000,Nb_pale,Nb_disc];
ub = [[0.15 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.15 0.15 0.15 0.15 0.15 0.15 0.15]*rtip,ones(1,Nb_disc)*pi/5,8500,Nb_pale,Nb_disc];
x0 = [c_init,twist_init,n_init,Nb_pale,Nb_disc];

%twist_init = twist_init*180/pi;

T0 = constraint_aerody_corde_twist(x0);
Lp0 = constraint_noice_corde_twist(x0);
P0 = fun_aerody_corde_twist(x0);

A = [];
b = [];
Aeq = [];
beq = [];

options = optimset('MaxIter',600,'TolFun',1.e-2);
[minparametre,P_finale] = fmincon(@fun_aerody_corde_twist,x0,A,b,Aeq,beq,lb,ub,nonlin,options);
T_finale = constraint_aerody_corde_twist(minparametre);
Lp_finale = constraint_noice_corde_twist(minparametre);
c_finale = minparametre(1:Nb_disc)*100; %cm
twist_finale = minparametre(Nb_disc+1:2*Nb_disc)*180/pi;
n_finale = minparametre(2*Nb_disc);

% figure(1)
Lp_liste = constraint_noice_corde_twist_figure(minparametre);
angle = 0:0.02*pi:pi*4/9;
% plot(angle*180/pi,Lp_liste,'r','LineWidth',2)
% xlabel('$\theta$ (deg)')
% ylabel('$L_p$ (dB)')
% xlim([0,90])
avLp = max(Lp_liste);

function [cineq,ceq] = constraint_nonlin(params)
    ceq = [];
    cineq = zeros(1,2);
    a = constraint_aerody_corde_twist(params);
    %b = constraint_noice_corde_twist(params);
    cineq(1,:) = -a + 4.5;
    %cineq(2,:) = b - 45;
end
end

