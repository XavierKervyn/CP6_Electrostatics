%% Parametres

clear all; clc
repertoire = '';                 % Chemin d'acces au code compile
executable = 'Exercice6_KervynLeMeur';        % Nom de l'executable
input      = 'configuration.in'; % Nom du fichier d'entree

format long;

nsimul = 20; % Nombre de simulations a faire

N       = round(logspace(1,4,nsimul));
trivial = true;

if (trivial)
    N1 = N; N2 = N;
else
    N1 = N; N2 = 2*N;
end

paramstr  = 'N1';               param     = N1;  
paramstr0 = 'N2';               param0    = N2;
paramstr1 = 'b';                param1    = 3.e-1 ;
paramstr2 = 'R';                param2    = 5.e-1;
paramstr3 = 'a0';               param3    = -3.e4;
paramstr4 = 'epsilon_r';        param4    = 4.e0;
paramstr5 = 'V0';               param5    = 2.2e2 ;
paramstr6 = 'MeshFactor';       param6    = 0.5;
paramstr7 = 'p';                param7    = 1.e0;
paramstr8 = 'trivial';          param8    = trivial;
paramstr9 = 'proportionalMesh'; param9    = false;

filename  = "N1=" + num2str(param(1),8);
for i = 2:nsimul
    filename(i)="N1=" + num2str(param(i),8);
end

lw=1; fs = 18; ms = 7;
set(groot, 'DefaultTextInterpreter', 'LaTeX');
set(groot, 'DefaultAxesTickLabelInterpreter', 'LaTeX');
set(groot, 'DefaultAxesFontName', 'LaTeX');
set(groot, 'DefaultLegendInterpreter', 'LaTeX');
set(groot, 'DefaultAxesBox', 'on');

%% Simulations
output = cell(4, nsimul);
%name=['_dataprof' '_phi' '_E_D' '_div_E_D'];
for i = 1:nsimul
    output{1,i} = [paramstr, '=', num2str(param(i),8), '_dataprof','.out'];
    output{2,i} = [paramstr, '=', num2str(param(i),8), '_phi','.out'];
    output{3,i} = [paramstr, '=', num2str(param(i),8), '_E_D','.out'];
    output{4,i} = [paramstr, '=', num2str(param(i),8), '_div_E_D','.out'];
    
    % Variant to scan N1 and N2 together:
    cmd=sprintf('%s%s %s %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g %s=%.15g outputDataProfiles=%s outputPotential=%s outputElectricDisplacementFields=%s outputDivergences=%s', repertoire, executable, input, paramstr, param(i), paramstr0, param0(i),paramstr1,param1,paramstr2,param2,paramstr3,param3,paramstr4,param4,paramstr5,param5,paramstr6,param6,paramstr7,param7,paramstr8,param8,paramstr9,param9,output{1,i},output{2,i},output{3,i},output{4,i});
    disp(cmd)   
    system(cmd); 
end
%% Analyse
% Parcours des resultats de toutes les simulations
phip = zeros(1,nsimul);

for i=1:nsimul
    Analyse(filename(i));
end


%% Analyse de convergence sur phi_0
phi_0 = (param2^2)/4 + param5; %solution analytique
figure('Name',"Convergence phi_0")
Errphi = zeros(1,nsimul);
for i=1:nsimul
    data = load([filename(i)+'_phi.out']);
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
paramstr6 = 'MeshFactor';       Mesh_       = 0.5;        %?? changer
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

%% Plot phi et E pour diff??rentes valeurs de N
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


%% Convergence de phi(r=b) pour diff??rentes valeurs

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
