global dtau 
global Nb_disc
dtau = 1e-5;
Nb_disc = 15;

c = minparametre(1:15);
n = minparametre(31);
twist = minparametre(16:30);
B = minparametre(32);
COOR_POLAR = [0,pi/4,1.5];
constraint_noice_direction(c,n,twist,B,COOR_POLAR)

tau = 0:dtau:1;
tstart_liste = zeros(1,length(c));
tend_liste = zeros(1,length(c));
pL_liste = zeros(length(c),length(tau));
t_liste = zeros(length(c),length(tau));
    for Rind = 1:(length(c))
        [pL_ele,tstart,tend,t] = noice_ele_corde_twist(c,n,twist,Rind,COOR_POLAR);
        if (~any(isnan(pL_ele)))
            pL_liste(Rind,:) = pL_ele;
        end
        tstart_liste(Rind) = tstart;
        tend_liste(Rind) = tend ;
        t_liste(Rind,:) = t;
        %plot(t,pL_ele);
        %hold on
    end
%legend('1','2','3','4','5','6','7','8','9','Location','northeastoutside','box','off');
tob = max(tstart_liste):dtau:min(tend_liste);
pL = zeros(1,length(tob));
    for i = 1:length(c)
        pL = pL + interp1(t_liste(i,:),pL_liste(i,:),tob);
    end

fe=1/dtau;
BPF = B*n/60; % blade passing frenquency
ptot = pL;
for i = 1:(B-1)
    tc = 1/BPF*i;
    indtc = round(tc*fe);
    pL2 = pL(indtc:end);
    tob2 = tob(1:length(pL2)); 
    ptot = ptot(1:length(tob2)) + pL2;
    % plot(tob2,pL2(1:length(tob2)),'k','LineWidth',1)
    % hold on;
end
% plot(tob2,pL(1:length(tob2)),'r-','LineWidth',1)
% plot(tob2,ptot,'b-','LineWidth',2)
% xlim([0,5*tc])
pref = 2e-5;
Lp = 10*log10(rms(ptot).^2/pref^2);
