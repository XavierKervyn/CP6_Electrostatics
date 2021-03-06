%% Paramètres: cas non trivial et question d.ii ///////////////////////////
ViewFormat;

trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0.5;
p_          = 1.e0;
propMesh_   = trivial_;

N1_ = N; N2_ = N;
if (propMesh_)
    N2_ = MeshFactor_*N;
end

filename2  = "N1_"+ num2str(N1_(1)) + "N2_" + num2str(N2_(1)) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i))+"N2_"+num2str(N2_(i)) ;
end

% Simulations 2 (non trivial)
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
    disp('Exercice6 configuration.in');   
    system('Exercice6 configuration.in'); 
end
    
% Question d.ii
% figure('Name',"Plot Divergence")
% for i=1:nsimul
%     data   = load(filename2(i)+'_div_E_D.out');
%     r      = data(:,1);
%     divD   = data(:,4);
%     
%     plot(r,divD,'.-','Linewidth',lw);
%     xlabel('$r$ [m]'); ylabel('$\nabla \cdot D_r$ [V/m$^2$]');
%     grid minor; hold on; set(gca,'fontsize',fs);
% end
% SaveIMG('PlotDivergence');

figure('Name',"verification of Gauss' law")
intervalle = 2;
for i=1:intervalle:nsimul
    data   = load(filename2(i)+'_div_E_D.out');
    r      = data(:,1);
    rholib = data(:,2);
    divD   = data(:,4);
    
    plot(r, abs((rholib-divD)/a0_) * 100,'.-','Linewidth',lw);
    hold on
    xlabel('$r$ [m]'); ylabel('$\chi_d$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
end
    steps = N(1:intervalle:nsimul);
    legendStrings = [string(steps),"$N$"];
    leg = legend(legendStrings,'Location','northwest','NumColumns',2);
    title(leg, '$N$')
SaveIMG('VerificationGaussLaw');

%étude de convergence sur le max de la différence
figure('Name',"Convergence Gauss' law")
MaxDiff  = zeros(1,nsimul);
for i=1:nsimul
    data   = load(filename2(i)+'_div_E_D.out');
    r      = data(:,1);
    rholib = data(:,2);
    divD   = data(:,4);
    MaxDiff(i) = max(abs((rholib-divD)/a0_) * 100);
end
    plot(1./N,MaxDiff,'x','Linewidth',lw);
    xlabel('$1/N$'); ylabel('Max($\chi_d$) [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
    absc = (1./N);
    P = polyfit(absc,MaxDiff,1);
    hold on
    xrange = (0:1.e-3:0.1);
    plot(xrange, P(1)*xrange + P(2),'--','Linewidth',lw,'HandleVisibility','off');
SaveIMG('ConvergenceGaussLaw');

% calcul de la densité de charges de polarisation
h1 = figure('Name',"plot rho pol");
intervalle = 2;
for i=1:intervalle:nsimul
    data   = load(filename2(i)+'_div_E_D.out');
    r      = data(:,1);
    divE    = data(:,3);
    divD   = data(:,4);
    rho_pol = divE - divD;
    
    plot(r, rho_pol,'.-','Linewidth',lw);
    hold on
end
    xlabel('$r$ [m]'); ylabel('$\rho_{b} / \varepsilon_0 $ [V/m$^2$]');
    grid minor; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul);
    legendStrings = [string(steps),"$N$"];
    leg = legend(legendStrings,'Location','northeast','NumColumns',2);
    title(leg, '$N$')
    xlim([r(1) r(end)]);
    ylim([min(rho_pol) max(rho_pol)]);
    MagInset(h1, -1, [0.27 0.33 0 15e4], [0.05 0.25 1.e6 2.e6], {'NW','SW';'NE','SE'});
    set(gca,'fontsize',fs-4);
SaveIMG('plotRhoPol');
        
figure('Name',"plot rho lib")
    data   = load(filename2(nsimul)+'_div_E_D.out');
    r      = data(:,1);
    rholib = data(:,2);
    plot(r, rholib,'.','Linewidth',lw);
    xlabel('$r$ [m]'); ylabel('$\rho_{f} / \varepsilon_0 $ [V/m$^2$]');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG('plotRhoLib');

%% Calcul des charges de polarisation & étude de convergence
nsimul = 20;

trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MF = [2 1 0.5];
p_          = 1.e0;
propMesh_   = trivial_;

N1_ = logspace(1,3,nsimul);
N2_ = [MF(1)*N1_; MF(2)*N1_; MF(3)*N1_;];

filename2 = strings(3,10);
Qpol = zeros(1,nsimul);
figure('Name','ConvergenceQpol')
for i=length(MF):length(MF)
    % Simulations
    MeshFactor_ = MF(i);
    filename2(i,1)    = "N1_"+ num2str(N1_(1)) + "N2_" + num2str(N2_(i,1));
    for j = 2:nsimul
        filename2(i,j)= "N1_"+ num2str(N1_(j)) + "N2_" + num2str(N2_(i,j));
    end
    for j = 1:nsimul
        %simulation
        N1_loc = N1_(j);
        N2_loc = N2_(i,j);
        writeConfig;
        disp('Exercice6 configuration.in');   
        system('Exercice6 configuration.in'); 
        
        %calcul
        data      = load(filename2(i,j)+'_div_E_D.out');
        r         = data(:,1);
        divE      = data(:,3);
        divD      = data(:,4);
        rho_pol   = divE - divD;
        Somme   = 0;
        for k=1:length(r)-1
            Somme = Somme + (r(k+1)-r(k))*(r(k)*rho_pol(k) + r(k+1)*rho_pol(k+1)) ;
        end
        Qpol(j) = 8.85e-12*pi*Somme;
        QpolAna = -0.5*8.85e-12*a0_*b_^2;
    end
        loglog(N1_+N2_(i,:),abs(Qpol-QpolAna)/QpolAna * 100,'x','Linewidth',lw);
        
end
    hold on
    FitLOGLOG(N1_+N2_(i,:),abs(Qpol-QpolAna)/QpolAna * 100,1);
    xlabel('$N_{total}$'); ylabel('$\chi_q$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG("ConvergenceQpolMesh=05");