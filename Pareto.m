clear;clc;close all;
n = 20; % nombre d'iteration
global Nb_disc
Nb_disc = 15;
taille_minparametre = 2*Nb_disc + 3;
Nb_pale = 4;
%n_init: RPM initial
n_init = 4000;

global dtau
dtau = 1e-5;

figure(1)
minparametre_liste = zeros(n,taille_minparametre+1);
avLp_liste = zeros(1,n);
P_liste = zeros(1,n);

for i = 1:n
    T_critere = false;
    [avLp,P,T_finale,minparametre] = Optimisation_corde_twist_fun(Nb_pale,n_init);
    if T_finale >= 4.5 && T_finale <= 4.8
        T_critere = true;
    end
    minparametre_liste(i,taille_minparametre+1) = T_critere;
    minparametre_liste(i,1:taille_minparametre) = minparametre;
    avLp_liste(i) = avLp;
    P_liste(i) = P;
    plot(avLp,P,'.','color','r','Markersize',43)
    text(avLp, P, sprintf('%d', i), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
    hold on;
    disp(['iteration',num2str(i)])
end

% for i = 1:n
%     leg_str{i} = ['data',num2str(i)];
% end
% legend(leg_str)
xlabel('$L_p$ (dB)','Interpreter','latex')
ylabel('P (W)','Interpreter','latex')
title(['Front de Pareto avec ' num2str(Nb_pale) ' pales'],"FontSize",20)

%%
ind_choisi = 6;
%minparametre_liste(ind_choisi,21) = true;
if minparametre_liste(ind_choisi,taille_minparametre+1) == false
    disp(['Error:Portance insuffisante!!!' 'La portance est ' num2str(constraint_aerody_corde_twist(minparametre_liste(ind_choisi,1:taille_minparametre))) ' N'])
else
    disp(['OK!!!' 'La portance est ' num2str(constraint_aerody_corde_twist(minparametre_liste(ind_choisi,1:taille_minparametre))) ' N'])

    minparametre = minparametre_liste(ind_choisi,1:taille_minparametre);
    filename = strcat('r1.5pale',num2str(Nb_pale),'n',num2str(n_init),'.mat');
    filepath = 'D:\Document\CoursENSTA\2A\PIE\Rapport final\minparametre\20240330\';
    save([filepath,filename],'minparametre')
    
    figure(2)
    Lp_liste = constraint_noice_corde_twist_figure(minparametre);
    angle = 0:0.02*pi:pi/2;
    plot(angle*180/pi,Lp_liste,'r','LineWidth',2)
    title('Bruit à différentes directions','FontSize',20)
    xlabel('$\theta$ (deg)','FontSize',18)
    ylabel('$L_p$ (dB)','FontSize',18)
    xlim([0,85])
    filename = strcat('noice-r1.5pale',num2str(Nb_pale),'n',num2str(n_init),'.jpg');
    filepath = 'D:\Document\CoursENSTA\2A\PIE\Rapport final\figures\';
    saveas(2,[filepath,filename])

    filename = strcat('r1.5pale',num2str(Nb_pale),'n',num2str(n_init),'.jpg');
    filepath = 'D:\Document\CoursENSTA\2A\PIE\Rapport final\figures\20240330\';
    saveas(1,[filepath,filename])
end
%%
retrace = true;
if retrace == true
    figure(1)
    for i = 1:n
        plot(avLp_liste(i),P_liste(i),'.','color','r','Markersize',43)
        text(avLp_liste(i), P_liste(i), sprintf('%d', i), 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center');
        hold on;
    end
end

% for i = 1:n
%     leg_str{i} = ['data',num2str(i)];
% end
% legend(leg_str)

xlabel('$L_p$ (dB)','Interpreter','latex')
ylabel('P (W)','Interpreter','latex')
title(['Front de Pareto avec ' num2str(Nb_pale) ' pales'],"FontSize",20)


