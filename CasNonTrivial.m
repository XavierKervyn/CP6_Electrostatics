%% Paramètres: cas non trivial avec facteur 2 /////////////////////////////
ViewFormat;

nsimul      = 20; 
N           = round(logspace(1,3,nsimul));
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = 1.e0;
propMesh_   = ~trivial_;

N1_ = N; N2_ = N;
if(propMesh_)
    N2_ = MeshFactor_*N;
end

filename2  = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

% Simulations 2 (non trivial)
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
    disp('Exercice6 configuration.in');   
    system('Exercice6 configuration.in'); 
end

%% Plot du epsilon_r et rho_lib
figure('Name','plot e_r et rho_lib')
    r1      = 0:1.e-3:0.3; r2 = 0.3:1.e-3:0.5; r  = [r1 r2];
    eps_r   = [ones(1,length(r1)) epsilon_r_*ones(1,length(r2))];
    rho_lib = [a0_*sin(3*pi*r1/b_) zeros(1,length(r2))];
    hold on
    yyaxis left
    plot(r, eps_r,'Linewidth',lw);
    xlabel('$r$ [m]');
    ylabel('$\varepsilon_r$');
    ylim([min(eps_r)-0.5 max(eps_r)+0.5]);

    yyaxis right
    plot(r, rho_lib,'Linewidth',lw);
    ylabel('$\rho_{f}/\varepsilon_0$ [V/m$^2$]');
    ylim([min(rho_lib)-1000 max(rho_lib)+1000]);
    hold off
    grid on; hold on; set(gca,'fontsize',fs);
SaveIMG("PlotEpsilonRrhoLib");
    
%% Analyse 2 (cas non trivial)
% for i=1:nsimul
%     Analyse(filename2(i));
% end

%% Plot phi et E pour différentes valeurs de N (cas non trivial)
figure('Name','plot phi non trivial')
    intervalle = 5;
    for i=1:intervalle:nsimul
        data = load(filename2(i)+'_phi.out');
        r = data(:,1);
        phi = data(:,2);
        plot(r,phi,'-','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = string(steps);
    leg = legend(legendStrings,'Location','southeast','NumColumns',2);
    title(leg, '$N$')
SaveIMG("PhiNonTrivialMesh"+num2str(MeshFactor_));
    
figure('Name','plot E non trivial')
    for i=1:intervalle:nsimul
        data = load(filename2(i)+'_E_D.out');
        r = data(:,1);
        E = data(:,2);
        plot(r,E,'-','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$E_r$ [V/m]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = string(steps);
    leg = legend(legendStrings,'Location','southeast','NumColumns',1);
    title(leg, '$N$')
SaveIMG("ENonTrivialMesh"+num2str(MeshFactor_));

%% Convergence de phi(r=b) pour différentes valeurs
%on refait une simulation plus détaillée pour une meilleure convergence
nsimul      = 20; 
N           = round(logspace(3,4,nsimul));
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = 1.e0;
propMesh_   = ~trivial_;

N1_ = N; N2_ = N;
if(propMesh_)
    N2_ = MeshFactor_*N;
end

filename2  = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

% Simulations 2 (non trivial)
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
    disp('Exercice6 configuration.in');   
    system('Exercice6 configuration.in'); 
end

figure('Name',"Convergence phi Ana r b")
phiAna = V0_ - (a0_*b_^2)/(12*pi) * log(b_/R_);
phirb  = zeros(1,nsimul); %valeur de phi en r=b
for i=1:nsimul
    data = load(filename2(i)+'_phi.out');
    r    = data(:,1);
    phirb(i)     = data(N1_(i)+1,2);
end
    err = abs((phirb - phiAna)/phiAna)*100;
    
    %plot 
    loglog(N, err,'x','Linewidth',lw,'HandleVisibility','off');
    hold on
    FitLOGLOG(N,err,1);
    xlabel('$N$; $N_2 = 2N_1$'); ylabel('$\chi_b$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG("ConvergencePhiAnaRBMesh"+num2str(MeshFactor_));

%% Plot de phi analytique
%on refait une simulation plus détaillée pour une meilleure convergence
nsimul      = 20; 
N           = round(logspace(1,3,nsimul));
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = 1.e0;
propMesh_   = ~trivial_;

N1_ = N; N2_ = N;
if(propMesh_)
    N2_ = MeshFactor_*N;
end

filename2  = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

% Simulations 2 (non trivial)
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
    disp('Exercice6 configuration.in');   
    system('Exercice6 configuration.in'); 
end

figure('Name','plot phi analytique')
    intervalle = 7;
    for i=1:intervalle:nsimul
        data = load(filename2(i)+'_phi.out');
        r = data(:,1);
        phi = data(:,2);
        fun = @(x)sin(3*pi*x/b_)/x;
        phiAna = zeros(size(r));
        for k=1:length(phiAna)
            if(r(k) < b_)
                phiAna(k) = V0_ - (a0_ * b_^2)/(9*pi^2)*( 3*pi/4*log(b_/R_) + integral(fun,b_,r(k),'ArrayValued',true) - sin((3*pi/b_).*r(k)));
            else
                phiAna(k) = V0_ - (a0_*b_^2)/(12*pi) * log(r(k)/R_);
            end
        end
        hold on
        plot(r,phi,'-','Linewidth',lw);
    end
    plot(r,phiAna,'r--','Linewidth',lw+1);
    xlabel('$r$ [m]'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = [string(steps) "$\phi_{ana}$"];
    leg = legend(legendStrings,'Location','southeast','NumColumns',2);
    title(leg, "$N_1$ ($N_2 =$"+num2str(MeshFactor_)+"$N_1$)")
SaveIMG("PhiAnalytiqueMesh"+num2str(MeshFactor_));

%erreur par rapport à phi analytique (on a seulement la dernière simulation
%numérique)
figure('Name','Erreur Phi Analytique')
    err = abs((phi-phiAna)./phiAna)*100;
    plot(r,err,'-','Linewidth',lw);
    xlabel('$r$ [m]'); ylabel('$\chi_a$ [\%]');
    grid on; hold on; set(gca,'fontsize',fs);
SaveIMG("ErreurPhiAnalytiqueMesh"+num2str(MeshFactor_)+"N="+num2str(N(nsimul)));