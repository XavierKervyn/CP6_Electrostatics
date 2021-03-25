%% Parametres /////////////////////////////////////////////////////////////
clear all; clc
ViewFormat;

nsimul = 20; N = round(logspace(1,4,nsimul));
trivial_    = true;
if (trivial_)
    N1_ = N; N2_ = N;
else
    N1_ = N; N2_ = 2*N;
end

b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0.5;
p_          = 1.e0;
propMesh_   = false;

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
for i=nsimul:nsimul
    Analyse(filename2(i));
end

%% Analyse: plot superposé
figure('Name','plot superposé phi')
for i=1:4:nsimul
    data    = load([filename2(i)+'_phi.out']);
    r       = data(:,1);
    phi     = data(:,2);
    phi_ana = (R_^2 - r.^2)/4 + V0_; %solution analytique
    plot(r,phi,'.','Linewidth',lw);
    hold on
end
    Pana = plot(r,phi_ana,'r--','Linewidth',lw);
    xlabel('$r [m]$'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:4:nsimul); legendStrings = [string(steps),"$\phi_{ana}$"];
    leg = legend(legendStrings,'Location','southwest','NumColumns',2);
    title(leg, 'Steps $N$') 
    
    % create a new pair of axes inside current figure [.65 .175 .25 .25]
    axes('position',[0.65 0.68 0.30 0.30])
    box on % put box around new pair of axes
for i=[13 17]
    data    = load([filename2(i)+'_phi.out']);
    r       = data(:,1);
    indexOfInterest = (r < 0.07); % range of r near perturbation
    phi     = data(:,2);
    phi_ana = (R_^2 - r.^2)/4 + V0_; %solution analytique
    plot(r(indexOfInterest),phi(indexOfInterest),'.','Linewidth',lw);
    hold on
end
    %plot de la solution analytique, avec la range de la dernière simu
    plot(r(indexOfInterest),phi_ana(indexOfInterest),'r--','Linewidth',lw) % plot on new axes
    axis tight
    SaveIMG('ComparaisonQualitativePhiTrivial');
    
figure('Name','plot superposé E & D')
    i = 5;
    data    = load([filename2(i)+'_E_D.out']);
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
    SaveIMG('ComparaisonQualitativeEDTrivial');

%% Analyse de convergence sur phi_0

phi_0 = (R_^2)/4 + V0_; %solution analytique
figure('Name',"Convergence phi_0")
Errphi = zeros(1,nsimul);
for i=1:nsimul
    data = load([filename2(i)+'_phi.out']);
    r = data(:,1);
    phi = data(:,2);
    Errphi(1,i) = abs(phi(1)-phi_0)/phi_0 *100;
end
    loglog(N, Errphi ,'x','Linewidth',lw,'HandleVisibility','off');
    xlabel('$N$'); ylabel('Relative error on $\phi(0)$ [\%]');
    grid on; hold on; set(gca,'fontsize',fs);
    FitLOGLOG(N,Errphi,1);
    SaveIMG('ConvergencePhi0');

%% Paramètres: cas non trivial ////////////////////////////////////////////
clear all; clc
ViewFormat;

nsimul = 20; N = round(logspace(3,4,nsimul));
trivial_    = false;
if (trivial_)
    N1_ = N; N2_ = N;
else
    N1_ = N; N2_ = 2*N;
end

b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0.5;
p_          = 1.e0;
propMesh_   = true;

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
for i=1:nsimul
    Analyse(filename2(i));
end
%% Plot phi et E pour différentes valeurs de N (cas non trivial)
figure('Name','plot phi')
    intervalle = 5;
    for i=1:intervalle:nsimul
        data = load([filename2(i)+'_phi.out']);
        r = data(:,1);
        phi = data(:,2);
        plot(r,phi,'.','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = string(steps);
    leg = legend(legendStrings,'Location','southwest','NumColumns',1);
    title(leg, '$N$')
    SaveIMG('PhiNonTrivial');
    
figure('Name','plot E')
    for i=1:intervalle:nsimul
        data = load([filename2(i)+'_E_D.out']);
        r = data(:,1);
        E = data(:,2);
        plot(r,E,'.','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$E_r$ [V/m]');
    grid on; hold on; set(gca,'fontsize',fs);
    steps = N(1:intervalle:nsimul); legendStrings = string(steps);
    leg = legend(legendStrings,'Location','southeast','NumColumns',1);
    title(leg, '$N$')
    SaveIMG('ENonTrivial');

%% Convergence de phi(r=b) pour différentes valeurs
figure('Name',"phi r b en fonction de 1./N")
phirb  = zeros(1,nsimul); %valeur de phi en r = b
rind   = zeros(1,nsimul); %valeur de r à l'indice trouvé
ind    = zeros(1,nsimul); %indices relevés
rfin   = zeros(1,nsimul); %r final
phifin = zeros(1,nsimul); %phi final
for i=1:nsimul
    data = load([filename2(i)+'_phi.out']);
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
    zprime = polyval(P, 1./N);
    %plot(1./N, zprime,'--','Linewidth',lw);
    legendStrings = string(P(2));
    xrange = (0:0.05:1)*1.e-3;
    plot(xrange, P(1)*xrange + P(2),'--','Linewidth',lw,'HandleVisibility','off');
    yas = plot(0,P(2),'r +','Markersize',15,'Linewidth',lw+2);
    leg = legend(legendStrings,'Location','northwest','NumColumns',2);
    title(leg, '$\phi_{as}$ [V]')
    xlabel('$1/N$'); ylabel('$\phi(r=b)$ [V]');
    grid minor; hold on; set(gca,'fontsize',fs);
    SaveIMG('PhiRBasymptotique');

figure('Name',"Convergence phi r b")
    FitLOGLOG(1./N,err,1);
    hold on
    loglog(1./N,err,'x','Linewidth',lw,'HandleVisibility','off');
    xlabel('$1/N$'); ylabel('Rel. error on $\phi_{as}(r=b)$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
    SaveIMG('ConvergencePhiRBasymptotique');
    
%% Paramètres: cas non trivial et question d.ii ///////////////////////////
clear all; clc
ViewFormat;

nsimul = 1; N = round(logspace(3,4,nsimul));
trivial_    = false;
N1_ = N;
if (trivial_)
    N2_ = N;
else
    N2_ = 2*N;
end

b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
MeshFactor_ = 0.5;
p_          = 1.e0;
propMesh_   = true;

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
    
% Question d.ii
figure('Name',"verification of Gauss' law")
for i=1:nsimul
    data   = load([filename2(i)+'_div_E_D.out']);
    r      = data(:,1);
    rholib = data(:,2);
    divD   = data(:,4);
    
    plot(r, abs((rholib-divD)/a0_) * 100,'.');
    xlabel('$1/N$'); ylabel('Rel. err. on $A$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
end

figure('Name',"verification of Gauss' law avant b")
for i=1:nsimul
    data   = load([filename2(i)+'_div_E_D.out']);
    r      = data(:,1);
    rholib = data(:,2);
    divD   = data(:,4);
    
    err = abs((rholib-divD)/a0_) * 100;
    plot(r(1:N1_(i)-2), err(1:N1_(i)-2),'.');
    xlabel('$1/N$'); ylabel('Rel. err. on $A$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
end

figure('Name',"verification of Gauss' law après b")
for i=1:nsimul
    data   = load([filename2(i)+'_div_E_D.out']);
    r      = data(:,1);
    rholib = data(:,2);
    divD   = data(:,4);
    
    err = abs((rholib-divD)/a0_) * 100;
    plot(r(N1_(i)+3:end-1), err(N1_(i)+3:end-1),'.');
    xlabel('$1/N$'); ylabel('Rel. err. on $A$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
end

% calcul de la densité de charges de polarisation
divE    = data(:,3);
rho_pol = divE - divD;
% epsr    = zeros(length(r));
% for i=1:length(r)
%    if (r(i) < b_)
%        epsr(i) = 1;
%    else
%        epsr(i) = 4;
%    end
% end
% rho_pol = (1-epsr).*data(:,3);

figure('Name',"plot rho pol")
        plot(r, rho_pol,'.','Linewidth',lw);
    hold on
%         plot(r, -3*divE,'.','Linewidth',lw)
        xlabel('$r$ [m]'); ylabel('$\rho_{pol} / \varepsilon_0 $ [V/m$^2$]');
        grid minor; hold on; set(gca,'fontsize',fs);

figure('Name',"verification rho pol law avant b")
    plot(r(3:N1_(i)-2), rho_pol(3:N1_(i)-2),'.');
    xlabel('$r [m]$'); ylabel('$\rho_{pol}$ [V/m$^2$]');
    grid on; set(gca,'fontsize',fs);

figure('Name',"verification rho pol après b")
    plot(r(N1_(i)+3:end-1), rho_pol(N1_(i)+3:end-1),'.');
    xlabel('$r [m]$'); ylabel('$\rho_{pol}$ [V/m$^2$]');
    grid on; set(gca,'fontsize',fs);      
        
figure('Name',"plot rho lib")
        plot(r, rholib,'.','Linewidth',lw);
        xlabel('$r$ [m]'); ylabel('$\rho_{lib} / \varepsilon_0 $ [V/m$^2$]');
        grid minor; hold on; set(gca,'fontsize',fs);
    