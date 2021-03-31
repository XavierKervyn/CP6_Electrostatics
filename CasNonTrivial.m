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

figure('Name',"phi r b en fonction de 1./N")
phirb  = zeros(1,nsimul); %valeur de phi en r = b
rind   = zeros(1,nsimul); %valeur de r à l'indice trouvé
ind    = zeros(1,nsimul); %indices relevés
rfin   = zeros(1,nsimul); %r final
phifin = zeros(1,nsimul); %phi final
for i=1:nsimul
    data = load(filename2(i)+'_phi.out');
    r    = data(:,1);
    [val,indice] = min(abs(r - b_));
    rind(i)      = r(indice+1);
    phirb(i)     = data(indice,2);
    ind(i)       = indice;
    rfin(i)      = r(end);
    phifin(i)    = data(end,2);
end
    %extrapolation
    P = polyfit(1./N, phirb,1); yintercept = P(2); 
    err = abs((phirb - yintercept)/yintercept)*100;
    
    %plot 
    plot(1./N, phirb,'x','Linewidth',lw,'HandleVisibility','off');
    hold on
    abscisse = 1./N; legendStrings = string(P(2)); xrange = (0:1.e-3:max(abscisse));
    plot(xrange, P(1)*xrange + P(2),'--','Linewidth',lw,'HandleVisibility','off');
    plot(0,P(2),'r +','Markersize',15,'Linewidth',lw+2); %marquer le y-intercept
    leg = legend(legendStrings,'Location','northwest','NumColumns',2);
    title(leg, '$\phi_{as}$ [V]')
    xlabel('$1/N$'); ylabel('$\phi(r=b)$ [V]');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG("PhiRBasymptotiqueMesh"+num2str(MeshFactor_));

figure('Name',"Convergence phi r b")
    FitLOGLOG(1./N,err,1);
    hold on
    loglog(1./N,err,'x','Linewidth',lw,'HandleVisibility','off');
    xlabel('$1/N$'); ylabel('Rel. error on $\phi_{as}(r=b)$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
    legend('Location','northwest','Numcolumns',1);
SaveIMG("ConvergencePhiRBasymptotiqueMesh"+num2str(MeshFactor_));