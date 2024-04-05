function [Cl,Cd,AoA_liste] = Viterna_Corrigan(c,rtip)

Cl = [];
Cd = [];
AoA_liste = [];
for alpha_deg = 15:0.5:45
    alpha = alpha_deg*pi/180;
    AR = rtip/mean(c);
    Cdmax = 1.11 + 0.018*AR;
    B1 = Cdmax;
    A1 = B1/2;
    Cls = 1.48; %NACA4412
    Cds = 0.07; %NACA4412
    alphas_deg = 15; %NACA4412
    alphas = alphas_deg*pi/180;
    A2 = (Cls - Cdmax*sin(alphas)*cos(alphas))*sin(alphas)/cos(alphas)^2;
    B2 = Cds - Cdmax*sin(alphas)^2/cos(alphas);
    
    Cl = [Cl,A1*sin(2*alpha) + A2*(cos(alpha)^2)/sin(alpha)];
    Cd = [Cd,B1*sin(alpha)^2 + B2*cos(alpha)];
    AoA_liste = [AoA_liste,alpha_deg];
end
end

% figure(1)
% plot(alpha_liste,Cl,'k-','LineWidth',2)
% hold on
% plot(alpha_liste,Cd,'r-','LineWidth',2)
% set(gca,'FontSize',15)
% grid on
% xlabel('$\alpha (^\circ)$')
% legend('$C_L$','$C_D$','Location','best')