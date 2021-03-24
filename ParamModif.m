%% Parametres

clear all; clc
repertoire = '';                 % Chemin d'acces au code compile
executable = 'Exercice6_KervynLeMeur';        % Nom de l'executable
input      = 'configuration_.in'; % Nom du fichier d'entree

format long;

nsimul = 20; % Nombre de simulations a faire

N       = round(logspace(1,4,nsimul));
trivial_ = true;

if (trivial_)
    N1_ = N; N2_ = N;
else
    N1_ = N; N2_ = 2*N;
end

b_  = 3.e-1 ;
R_  = 5.e-1;
a0_ = -3.e4;
epsilon_r_ = 4.e0;
V0_  = 2.2e2 ;
MeshFactor_   = 0.5;
p_  = 1.e0;
propMesh_  = false;

lw=1; fs = 18; ms = 7;
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultAxesBox', 'on');

%% Simulations

filename2  = "N1_"+ num2str(N1_(1),8) + "N2_" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename2(i)="N1_"+ num2str(N1_(i),8)+"N2_"+num2str(N2_(i),8) ;
end

for i = 1:nsimul
    N1_loc = N1_(i);
    N2_loc = N2_(i);
    writeConfig;
    disp('Exercice6_KervynLeMeur configuration_.in');   
    system('./Exercice6_KervynLeMeur configuration_.in'); 
end

% disp('Exercice6_KervynLeMeur configuration.in');   
% system('Exercice6_KervynLeMeur configuration.in');
%% Analyse
% Parcours des resultats de toutes les simulations
phip = zeros(1,nsimul);

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

%% Figures
figure
if strcmp(input,'trivial.in')
    R = 0.12;
    loglog(N,abs(phip-R^2/4)/(R^2/4),'k+')
    xlabel('N')
    ylabel('|\phi(0) - R^2/4| / (R^2/4) [V]')
elseif strcmp(input,'nontrivial.in')
    plot(1./N.^2,phip,'k+')
    xlabel('1/N^2 (N=N_1=N_2/2)')
    ylabel('\phi(b) [V]')
end
grid on

%% Cas non trivial

trivial = false;
if (trivial)
    N1 = N; N2 = N;
else
    N1 = N; N2 = 2*N;
end

paramstr  = 'N1';               N1_         = N1;  
paramstr0 = 'N2';               N2_         = N2;
paramstr1 = 'b';                b_          = 3.e-1 ;
paramstr2 = 'R';                R_          = 5.e-1;
paramstr3 = 'a0';               a0_         = -3.e4;
paramstr4 = 'epsilon_r';        epsilon_r_  = 4.e0;
paramstr5 = 'V0';               V0_         = 2.2e2 ;
paramstr6 = 'MeshFactor';       Mesh_       = 0.5;        %à changer
paramstr7 = 'p';                p_          = 1.e0;
paramstr8 = 'trivial';          trivial_    = trivial;
paramstr9 = 'proportionalMesh'; propMesh_   = true;

filename  = "N1="+ num2str(N1_(1),8) + "N2=" + num2str(N2_(1),8) ;
for i = 2:nsimul
    filename(i)="N1="+ num2str(N1_(i),8)+"N2="+num2str(N2_(i),8) ;
end

lw=1; fs = 18; ms = 7;
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultAxesBox', 'on');

%% Simulations
output = cell(4, nsimul);
% for i = 1:nsimul
%     output{1,i} = [paramstr, '=', num2str(N1_(i),8),paramstr0,'=',num2str(N2_(i),8), '_dataprof','.out'];
%     output{2,i} = [paramstr, '=', num2str(N1_(i),8),paramstr0,'=',num2str(N2_(i),8), '_phi','.out'];
%     output{3,i} = [paramstr, '=', num2str(N1_(i),8),paramstr0,'=',num2str(N2_(i),8), '_E_D','.out'];
%     output{4,i} = [paramstr, '=', num2str(N1_(i),8),paramstr0,'=',num2str(N2_(i),8), '_div_E_D','.out'];
%     
%     % Variant to scan N1 and N2 together:
%     cmd=sprintf('%s%s %s %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g %s=%.18g outputDataProfiles=%s outputPotential=%s outputElectricDisplacementFields=%s outputDivergences=%s', repertoire, executable, input, paramstr, N1_(i), paramstr0, N2_(i),paramstr1,b_,paramstr2,R_,paramstr3,a0_,paramstr4,epsilon_r_,paramstr5,V0_,paramstr6,Mesh_,paramstr7,p_,paramstr8,trivial_,paramstr9,propMesh_,output{1,i},output{2,i},output{3,i},output{4,i});
%     disp(cmd)   
%     system(cmd); 
% end


%% Analyse
% Parcours des resultats de toutes les simulations
phip = zeros(1,nsimul);

for i=1:nsimul
    Analyse(filename(i));
end

%% Plot phi et E pour différentes valeurs de N
figure('Name','plot phi')
    for i=1:nsimul
        data = load([filename(i)+'_phi.out']);
        r = data(:,1);
        phi = data(:,2);
        plot(r,phi,'-','Linewidth',lw);
        hold on
    end
    xlabel('$r$ [m]'); ylabel('$\phi$ [V]');
    grid on; hold on; set(gca,'fontsize',fs);
    
figure('Name','plot E')
    for i=1:nsimul
        data = load([filename(i)+'_E_D.out']);
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
    data = load([filename(i)+'_phi.out']);
    r   = data(:,1);
    [val,indice] = min(abs(r - b_));
    phi = data(:,2);
    phirb(i) = phi(indice);
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
