%% Parametres /////////////////////////////////////////////////////////////
clear; clc
ViewFormat;

nsimul      = 20; 
N           = round(logspace(1,4,nsimul));
trivial_    = true;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0;
p_          = 1.e0;
propMesh_   = ~trivial_;

N1_ = N; N2_ = N;
if(propMesh_)
    N2_ = MeshFactor*N;
end

filename2 = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

% Simulations
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
    disp('Exercice6 configuration.in');   
    system('Exercice6 configuration.in');
end

%% Analyse 1: plot de tous les graphes
% for i=1:nsimul
%     Analyse(filename2(i));
% end

%% Analyse: plot superposé
h1 = figure('Name','plot superposé phi');
intervalle = 4;
for i=1:intervalle:nsimul
    data    = load(filename2(i)+'_phi.out');
    r       = data(:,1);
    phi     = data(:,2);
    phi_ana = (R_^2 - r.^2)/4 + V0_; %solution analytique
    plot(r,phi,'-','Linewidth',lw);
    hold on
end
    Pana = plot(r,phi_ana,'r--','Linewidth',lw);
    xlabel('$r$ [m]'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = [string(steps),"$\phi_{ana}$"];
    leg = legend(legendStrings,'Location','southwest','NumColumns',3);
    title(leg, '$N$')
    xlim([r(1) r(end)]);
    ylim([min(phi_ana) max(phi_ana)]);
    MagInset(h1, -1, [0.0 0.03 220.0623 max(phi_ana)], [0.12 0.25 220.025 220.04], {'NW','NW';'NE','NE'});
    set(gca,'fontsize',fs-4);
SaveIMG('ComparaisonQualitativePhiTrivial');
    
h1 = figure('Name','plot superposé E & D');
    i = 4;
    data    = load(filename2(i)+'_E_D.out');
    r       = data(:,1);
    E       = data(:,2);
    D       = data(:,3);
        plot(r,E,'+','Linewidth',lw);
    hold on
        plot(r,D,'x','Linewidth',lw);
    xlabel('$r$ [m]'); ylabel('$E_r, D_r/\varepsilon_0$ [V/m]');
    grid on; hold on; set(gca,'fontsize',fs);
    leg = legend('$E_r$','$D/\varepsilon_0$','Location','southeast','NumColumns',1);
    title(leg,"$N=$ "+num2str(N(i)));
SaveIMG("ComparaisonQualitativeEDTrivial");

%% Analyse de convergence sur phi_0
phi_0 = (R_^2)/4 + V0_; %solution analytique
figure('Name',"Convergence phi_0")
Errphi = zeros(1,nsimul);
for i=1:nsimul
    data = load(filename2(i)+'_phi.out');
    r = data(:,1);
    phi = data(:,2);
    Errphi(1,i) = abs(phi(1)-phi_0)/phi_0 *100;
end
    loglog(N, Errphi ,'x','Linewidth',lw,'HandleVisibility','off');
    xlabel('$N$'); ylabel('Relative error on $\phi(0)$ [\%]');
    grid on; hold on; set(gca,'fontsize',fs);
    FitLOGLOG(N,Errphi,1);
    legend('Location','northeast','Numcolumns',2);
SaveIMG('ConvergencePhi0');