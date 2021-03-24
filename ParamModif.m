%% Parametres
clear all; clc
format long;

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

lw=1; fs = 18; ms = 7;
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultAxesBox', 'on');

filename2 = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

% Simulations
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
%     disp('Exercice6_KervynLeMeur configuration_.in');   
%     system('Exercice6_KervynLeMeur configuration_.in'); 
    disp('Exercice6 configuration_.in');   
    system('Exercice6 configuration_.in');
end

% disp('Exercice6_KervynLeMeur configuration.in');   
% system('Exercice6_KervynLeMeur configuration.in');
%% Analyse 1: plot de tous les graphes
for i=1:nsimul
    Analyse(filename2(i));
end
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
    loglog(N, Errphi ,'x','Markersize',ms,'HandleVisibility','off');
    xlabel('$N$'); ylabel('Relative error on $\phi_0$ [\%]');
    grid on; hold on; set(gca,'fontsize',fs);
    P = polyfit(log(N), log(Errphi),1);
    z = polyval(P,log(N));
    loglog(N,exp(z),'--','Linewidth',lw);
    pentes = P(1); legendStrings = string(pentes);
    leg = legend(legendStrings,'Location','southwest','NumColumns',2);
    title(leg, 'Linear fit: slope')  

%% Paramètres: cas non trivial
clear all; clc
format long;

nsimul = 10; N = round(logspace(5,6,nsimul));
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

lw=1; fs = 18; ms = 7;
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultAxesBox', 'on');

filename2  = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

% Simulations 2 (non trivial)
for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
%     disp('Exercice6_KervynLeMeur configuration_.in');   
%     system('Exercice6_KervynLeMeur configuration_.in'); 
    disp('Exercice6 configuration_.in');   
    system('Exercice6 configuration_.in'); 
end

% disp('Exercice6_KervynLeMeur configuration.in');   
% system('Exercice6_KervynLeMeur configuration.in');
%% Analyse 2 (cas non trivial)

for i=1:nsimul
    Analyse(filename2(i));
end

%% Plot phi et E pour différentes valeurs de N (cas non trivial)
figure('Name','plot phi')
    for i=1:nsimul
        data = load([filename2(i)+'_phi.out']);
        r = data(:,1);
        phi = data(:,2);
        plot(r,phi,'-','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    
figure('Name','plot E')
    for i=1:nsimul
        data = load([filename2(i)+'_E_D.out']);
        r = data(:,1);
        E = data(:,2);
        plot(r,E,'-','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$E$ [V/m]');
    grid on; hold on; set(gca,'fontsize',fs);


%% Convergence de phi(r=b) pour différentes valeurs

figure('Name',"Convergence phi r b")
phirb = zeros(1,nsimul); %valeur de phi en r = b
for i=1:nsimul
    data = load([filename2(i)+'_phi.out']);
    r   = data(:,1);
    [val,indice] = min(abs(r - b_));
    phirb(i) = data(indice,2);
end
    loglog(N, phirb,'x','Markersize',ms,'HandleVisibility','off');
    xlabel('$N$'); ylabel('$\phi(r=b)$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
%     P = polyfit(log(N), log(phi),1);
%     z = polyval(P,log(N));
%     loglog(N,exp(z),'--','Linewidth',lw);
%     pentes = P(1); legendStrings = string(pentes);
%     leg = legend(legendStrings,'Location','southwest','NumColumns',2);
%     title(leg, 'Linear fit: slope')  
