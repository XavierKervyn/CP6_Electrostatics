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
    disp('Exercice6 configuration_.in');   
    system('Exercice6 configuration_.in');
end

%% Analyse 1: plot de tous les graphes
% for i=nsimul:nsimul
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
    xlabel('$r [m]$'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = [string(steps),"$\phi_{ana}$"];
    leg = legend(legendStrings,'Location','southwest','NumColumns',3);
    title(leg, 'Steps $N$')
    xlim([r(1) r(end)]);
    ylim([min(phi_ana) max(phi_ana)]);
    %zone de zoom
    MagInset(h1, -1, [0.0 0.03 220.0623 max(phi_ana)], [0.062 0.23 220.025 220.045], {'NW','NW';'NE','NE'});
    %save les images
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
    xlabel('$r [m]$'); ylabel('$E_r, D_r/\varepsilon_0$ [V/m]');
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

%% Paramètres: cas non trivial avec facteur 2 /////////////////////////////
clear; clc
ViewFormat;

nsimul      = 20; 
N           = round(logspace(2,4,nsimul));
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 2;
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
    disp('Exercice6 configuration_.in');   
    system('Exercice6 configuration_.in'); 
end
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
    leg = legend(legendStrings,'Location','southwest','NumColumns',1);
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

%% Paramètres: cas non trivial avec facteur 0.5 ///////////////////////////
clear; clc
ViewFormat;

trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0.5;
p_          = 1.e0;
propMesh_   = ~trivial_;
nsimul      = 20; 
N           = round(logspace(2,4,nsimul));

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
    disp('Exercice6 configuration_.in');   
    system('Exercice6 configuration_.in'); 
end
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
    leg = legend(legendStrings,'Location','southwest','NumColumns',1);
    title(leg, '$N$')
SaveIMG("PhiNonTrivialMesh05");
    
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
SaveIMG("ENonTrivialMesh05");

%% Convergence de phi(r=b) pour différentes valeurs
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
    yas = plot(0,P(2),'r +','Markersize',15,'Linewidth',lw+2);
    leg = legend(legendStrings,'Location','northwest','NumColumns',2);
    title(leg, '$\phi_{as}$ [V]')
    xlabel('$1/N$'); ylabel('$\phi(r=b)$ [V]');
    grid minor; hold on; set(gca,'fontsize',fs);
    SaveIMG("PhiRBasymptotiqueMesh05");

figure('Name',"Convergence phi r b")
    FitLOGLOG(1./N,err,1);
    hold on
    loglog(1./N,err,'x','Linewidth',lw,'HandleVisibility','off');
    xlabel('$1/N$'); ylabel('Rel. error on $\phi_{as}(r=b)$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
    legend('Location','northwest','Numcolumns',1);
    SaveIMG("ConvergencePhiRBasymptotiqueMesh05");
    
%% Question d.ii
Questiond2.m;
    