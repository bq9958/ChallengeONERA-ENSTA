clear;
rtip = 0.12;
Nb_disc = 15;
load('minparametre.mat')
params = minparametre;
c = params(1:Nb_disc);
n = params(2*Nb_disc+1);
twist = params(Nb_disc+1:2*Nb_disc);


set(groot,'defaulttextinterpreter','latex');  
    set(groot, 'defaultAxesTickLabelInterpreter','latex');  
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    
    
    % 1) Fixer une valeur initiale pour l'angle phi entre le plan de rotation
    % et la direction de l'écoulement incident W (implique des valeurs de u et W)
    % Parametres constants
    %rtip = 0.15; % Rayon du rotor (m)
    B = 1;
    r_hub = 0.0125; % hub radius
    r = linspace(0.1,0.99,Nb_disc)*rtip;
    %r = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99]*rtip;
    %n = 5000; % RPM
    rho = 1.225; % Masse volumique air 15 degrees Celsius pressure = 1 atm
    
    omega = n*2*pi/60;
    Vx = 0;
    Vy = omega*r;
    %c = [0.17 0.22 0.23 0.22 0.2 0.17 0.15 0.1 0.07]*rtip; % Repartition de la corde en fonction de la irition en r, en portion de rtip
    sigma = B*c./(2*pi*r); % "Solidite" locale
    % F = 1; % tip loss factor
    fsave = 'polar_0to15deg_NACA4412';
    load(fsave,'AoA','CL','CD')
    %twist = [45 34 29 23 18 16 14 13 12]*pi/180; % twist angle in rad noted theta in Ning (2021)
    
    [Cl_liste,Cd_liste,AoA_liste] = Viterna_Corrigan(c,rtip);
    AoA = [AoA,AoA_liste];
    CL = [CL,Cl_liste];
    CD = [CD,Cd_liste];
    
    % figure(1)
    % plot(AoA,CL,'k-','LineWidth',2)
    % hold on
    % plot(AoA,CD,'r-','LineWidth',2)
    % set(gca,'FontSize',15)
    % grid on
    % xlabel('$\alpha (^\circ)$')
    % legend('$C_L$','$C_D$','Location','best')
    
    % On fixe phi, puis on deduit W.
    
    % phi_deg0 = 45; % arbitraire
    % phi0 = phi_deg0*pi/180;
    
    
    % W0 = Vy/cos(phi0);
    % u = W*sin(phi); ?
    
    % search in quadrant 1
    phi_test = linspace(1e-6,pi/2,91);
    % dphi = phi_test(2)-phi_test(1);
    for ir = 1:length(r)
        % ir = 8;
        
        % tip loss correction
        F_tip = 2/pi*acos( exp(-B/2*(rtip-r(ir))./(r(ir)*abs(sin(phi_test)))) );
        % hub loss correction
        F_hub = 2/pi*acos( exp(-B/2*(r(ir)-r_hub)./(r_hub*abs(sin(phi_test)))) );
        % tip and hub loss factor
        F = F_tip.*F_hub;
    
        
    %     phi_test = linspace(0,twist(ir),21);
        alpha_test = twist(ir) - phi_test; % angle d'attaque = phi - theta
    %     CL_test = interp1(AoA,CL,alpha_test*180/pi,'linear','extrap'); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
    %     CD_test = interp1(AoA,CD,alpha_test*180/pi,'linear','extrap'); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
        CL_test = interp1(AoA,CL,alpha_test*180/pi); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
        CD_test = interp1(AoA,CD,alpha_test*180/pi); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
        Cn_test = CL_test.*cos(phi_test)-CD_test.*sin(phi_test);   % coefficient normal (dans l'axe du rotor)
        u_test = (sigma(ir)*Cn_test*Vy(ir))./(4*F.*sin(phi_test).*cos(phi_test)); % Vitesse induite (tangentielle)
        residual_phi_test = sin(phi_test)./u_test - cos(phi_test)/Vy(ir);
        % phi_test_final = atan(u_test./Vy(ir)); % angle phi "a la fin" de l'iteration
        
        % find lower and upper bounds for phi (sign change)
    %     ind = find(abs(diff(sign((phi_test_final-phi_test)))) == 2);
        ind = find(abs(diff(sign(residual_phi_test))) == 2);
        if length(ind) == 0 % no sign change found
            phiL(ir) = NaN;
            phiU(ir) = NaN;
        else
            phiL(ir) = phi_test(ind(1));
            phiU(ir) = phi_test(ind(1)+1);
        end
        clear ind
    end
    
    % initialiser les variables
    niter_max = 20;
    % tolU_pc = 1; % tolerance on u in percent
    tolU_pc = 0.5; % tolerance on u in percent
    alpha_iter = NaN(length(r),niter_max);
    CL_iter = NaN(length(r),niter_max);
    CD_iter = NaN(length(r),niter_max);
    Cn_iter = NaN(length(r),niter_max);
    u_iter = NaN(length(r),niter_max);
    phi_iter = NaN(length(r),niter_max);
    residual_phi_iter = NaN(length(r),niter_max);
    CL_final = NaN(1,Nb_disc);
    CD_final = NaN(1,Nb_disc);
    u_final = NaN(1,Nb_disc);
    phi_final = NaN(1,Nb_disc);
    W_final = NaN(1,Nb_disc);
    
    stall=0;
    
    % trouver le resultat par methode de dichotomie
    ind_segments = find(isnan(phiL) == 0);
    % for ir = 1:length(r) % Boucle sur les segments en r
    for ir = ind_segments % Boucle sur les segments en r avec [phiL phiU] réels
        %disp(['Segment ',num2str(ir),': r = ',num2str(r(ir)*1000,'%.0f'),'mm'])
    
        ite=0;
        convergence = 0;
        
    %     phi_ite0 = phiL(ir); % angle phi initial = AoA of 5 degrees
    %     phi_ite0 = twist(ir) - 5*pi/180; % angle phi initial = AoA of 5 degrees
    %     phi_ite0=twist(ir); % angle phi initial = zero incidence
        
        % calcul de u et de phi a partir de phiL
        F_tip = 2/pi*acos( exp(-B/2*(rtip-r(ir))./(r(ir)*abs(sin(phiL(ir))))) );
        F_hub = 2/pi*acos( exp(-B/2*(r(ir)-r_hub)./(r_hub*abs(sin(phiL(ir))))) );
        F = F_tip.*F_hub;
        alpha_L = twist(ir) - phiL(ir); % angle d'attaque = phi - theta
        CL_L = interp1(AoA,CL,alpha_L*180/pi,'linear','extrap'); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
        CD_L = interp1(AoA,CD,alpha_L*180/pi,'linear','extrap'); % coefficient de trainee recalcule par interpolation pour l'incidence exacte consid�r�e
        Cn_L = CL_L*cos(phiL(ir))-CD_L*sin(phiL(ir));   % coefficient normal (dans l'axe du rotor)
        u_L = sign(phiL(ir))*(sigma(ir)*Cn_L*Vy(ir))/(4*F.*sin(phiL(ir))*cos(phiL(ir))); % Vitesse induite (tangentielle)   
        residual_phiL = sin(phiL(ir))/u_L - cos(phiL(ir))/Vy(ir);
        % phiL_iter(ir) = atan(uL_iter(ir)/Vy(ir)); % angle phi "a la fin" de l'iteration
        
        % calcul de u et de phi a partir de phiU
        F_tip = 2/pi*acos( exp(-B/2*(rtip-r(ir))./(r(ir)*abs(sin(phiU(ir))))) );
        F_hub = 2/pi*acos( exp(-B/2*(r(ir)-r_hub)./(r_hub*abs(sin(phiU(ir))))) );
        F = F_tip.*F_hub;
        alpha_U = twist(ir) - phiU(ir); % angle d'attaque = phi - theta
        CL_U = interp1(AoA,CL,alpha_U*180/pi,'linear','extrap'); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
        CD_U = interp1(AoA,CD,alpha_U*180/pi,'linear','extrap'); % coefficient de trainee recalcule par interpolation pour l'incidence exacte consid�r�e
        Cn_U = CL_U*cos(phiU(ir))-CD_U*sin(phiU(ir));   % coefficient normal (dans l'axe du rotor)
        u_U = sign(phiU(ir))*(sigma(ir)*Cn_U*Vy(ir))/(4*F.*sin(phiU(ir))*cos(phiU(ir))); % Vitesse induite (tangentielle)   
        residual_phiU = sin(phiU(ir))/u_U - cos(phiU(ir))/Vy(ir);
        % phiU_iter(ir) = atan(uU_iter(ir)/Vy(ir)); % angle phi "a la fin" de l'iteration
        
        % methode de dichotomie
        while ((convergence == 0) & (ite < niter_max))
            % while (abs((u_ite(ir,ite)-u_ite(ir,ite-1))/u_ite(ir,ite)) > 0.01) & (ite < niter_max)
    
            ite = ite+1;
    
            % calcul de u et de phi a partir de (phiL+phiU)/2
            phi_iter(ir,ite) = (phiL(ir) + phiU(ir))/2;
            F_tip = 2/pi*acos( exp(-B/2*(rtip-r(ir))./(r(ir)*abs(sin(phi_iter(ir))))) );
            F_hub = 2/pi*acos( exp(-B/2*(r(ir)-r_hub)./(r_hub*abs(sin(phi_iter(ir))))) );
            F = F_tip.*F_hub;
            alpha_iter(ir,ite) = twist(ir) - phi_iter(ir,ite); % angle d'attaque = phi - theta
            CL_iter(ir,ite) = interp1(AoA,CL, alpha_iter(ir,ite)*180/pi,'linear','extrap'); % coefficient de portance recalcule par interpolation pour l'incidence exacte consid�r�e
            CD_iter(ir,ite) = interp1(AoA,CD, alpha_iter(ir,ite)*180/pi,'linear','extrap'); % coefficient de trainee recalcule par interpolation pour l'incidence exacte consid�r�e
            Cn_iter(ir,ite) = CL_iter(ir,ite)*cos(phi_iter(ir,ite))-CD_iter(ir,ite)*sin(phi_iter(ir,ite));   % coefficient normal (dans l'axe du rotor)
            u_iter(ir,ite) = sign(phi_iter(ir,ite))*(sigma(ir)*Cn_iter(ir,ite)*Vy(ir))/(4*F*sin(phi_iter(ir,ite))*cos(phi_iter(ir,ite))); % Vitesse induite (tangentielle)
            residual_phi_iter(ir,ite) = sin(phi_iter(ir,ite))/u_iter(ir,ite) - cos(phi_iter(ir,ite))/Vy(ir);
            % phi_final(ir) = atan(u_iter/Vy(ir)); % angle phi "a la fin" de l'iteration
            
            if sign(residual_phi_iter(ir,ite)) == sign(residual_phiL) % solution between phi_test and phiU
                if (abs(u_iter(ir,ite)-u_U)/u_U*100 < tolU_pc) % convergence found
                    convergence = 1;
    %                 phi_final(ir) = atan(u_iter(ir,ite)/Vy(ir));
    %                 u_final(ir) = u_iter(ir,ite);
                    %disp(['Convergence found: : u = ',num2str(u_iter(ir,ite),'%.1f'),'m/s - phi = ',num2str(phi_iter(ir,ite)*180/pi,'%.1f'),' deg'])
    %                 disp(['Iteration ',num2str(ite),': u = ',num2str(u_ite(ir,ite),'%.1f'),'m/s - phi = ',num2str(phi_ite(ir,ite)*180/pi,'%.1f'),' deg'])
                else
                    phiL(ir) = phi_iter(ir,ite);
                    residual_phiL = residual_phi_iter(ir,ite);
                    u_L = u_iter(ir,ite);
                end
            else % solution between phiL and phi_test
                if (abs(u_iter(ir,ite)-u_L)/u_L*100 < tolU_pc) % convergence found
                    convergence = 1;
    %                 phi_final(ir) = atan(u_iter(ir,ite)/Vy(ir));
    %                 u_final(ir) = u_iter(ir,ite);
    %                 disp(['Convergence found: : u = ',num2str(u_final(ir),'%.1f'),'m/s - phi = ',num2str(phi_final(ir)*180/pi,'%.1f'),' deg'])
                    %disp(['Convergence found: : u = ',num2str(u_iter(ir,ite),'%.1f'),'m/s - phi = ',num2str(phi_iter(ir,ite)*180/pi,'%.1f'),' deg'])
                else
                    phiU(ir) = phi_iter(ir,ite);
                    residual_phiU = residual_phi_iter(ir,ite);
                    u_U = u_iter(ir,ite);
                end            
            end
        end
    
        if convergence == 0
            stall = ir;
        end
    
        CL_final(ir) = CL_iter(ir,ite);
        CD_final(ir) = CD_iter(ir,ite);
        u_final(ir) = u_iter(ir,ite);
        phi_final(ir) = phi_iter(ir,ite);
        W_final(ir) = Vy(ir)./cos(phi_final(ir));
        
    end
    
    % CLs = CL_final(stall+1);
    % CDs = CD_final(stall+1);
    % alphas = twist(stall) - phi_final(stall);
    % 
    % for i=1:stall
    %     AR = rtip/c(i);
    %     B1 = 1.11 + 0.018*AR;
    %     Cdmax = B1;
    %     A2 = CLs - Cdmax*sin(alphas)*cos(alphas);
    %     B2 = CDs - CDs*(sin(alphas)^2)/cos(alphas);
    %     alpha = twist(i) - phi
    % 
    %     CL(i) = A1*sin(2*) 
    % end
    
    alpha_final = twist - phi_final;
   
    
    %% Calcul du coefficient de poussee C_T et de moment C_Q
    Cn_final = CL_final.*cos(phi_final) - CD_final.*sin(phi_final);   
    Ct_final = CL_final.*sin(phi_final) + CD_final.*cos(phi_final);  
    qx = 0.5*rho*(W_final.^2); % dynamic pressure (Pa)
    dT = B*Cn_final.*qx.*c;
    dD = B*Ct_final.*qx.*c;
    dQ = B.*r.*Ct_final.*qx.*c;
    
    % r_hub = 0.0125; % hub radius
    % r = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99]*rtip;
    r_middle = (r(1:end-1)+r(2:end))/2; % segment centers
    dr = zeros(size(r));
    dr(1) = r_middle(1)-r_hub;
    dr(2:end-1) = r_middle(2:end)-r_middle(1:end-1);
    dr(end) = rtip-r_middle(end);