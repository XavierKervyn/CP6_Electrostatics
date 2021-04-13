%% Cas Trivial
clear; clc
CasTrivial;

%% Cas Non trivial
clear; clc
MeshFactor_ = 2;
CasNonTrivial;

%% Etude N1 vs N2
clear; clc
nsimul  = 20; 
MF = 5;
printname = "N1=5N2";  MeshFactor_ = MF;   NonTrivialN1N2;
printname = "N1=02N2"; MeshFactor_ = 1/MF; NonTrivialN1N2;

%% Etude de convergence avec paramètre t
clear; clc
Ntotal = 1000;
tmin = 0.5; tmax = 0.99;
t = linspace(tmin,tmax,nsimul);
N1_ = t*Ntotal; N2_ = (1-t)*Ntotal;

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

figure('Name','plot param t')
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
    %plot 
    loglog(t, phirb,'x--','Linewidth',lw,'HandleVisibility','off');
    xlabel("$t$"); ylabel('$\phi(r=b)$ [V]');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG("PhiRBasymptotiqueMeshParamt");

%% Surface N1//N2
clear; clc;
nsimul = 20; DessineSurface(nsimul);

%% Question d.ii
clear all; clc
nsimul = 20;
N = round(logspace(1,3,nsimul));
Questiond2;

%% BONUS - variation de p.
clear; clc;
ViewFormat;
MeshFactor_ = 2;
% Convergence de phi(r=b) pour différentes valeurs
nsimul      = 10; np = 5;
N           = round(logspace(3,4,nsimul));
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = linspace(0,1,np);
propMesh_   = ~trivial_;

N1_ = N; N2_ = N;
if(propMesh_)
    N2_ = MeshFactor_*N;
end

filename2 = strings(nsimul,np);
for i = 1:nsimul
    for j = 1:np
    filename2(i,j)="N1_"+ num2str(N1_(i))+"N2_"+num2str(N2_(i)) + "p_" + num2str(p_(j));
    end
end

% Simulations
for i = 1:nsimul
    for j=1:np
        N1_loc = N1_(i);
        N2_loc = N2_(i);
        p_loc  = p_(j);
        writeConfigP;
        disp('Exercice6 configuration.in');   
        system('Exercice6 configuration.in');
    end
end

figure('Name',"Convergence phi RBP")
phiAna = V0_ - (a0_*b_^2)/(12*pi) * log(b_/R_);
phirb  = zeros(nsimul,np); %valeur de phi en r=b
for i=1:nsimul
    for j=1:np
        data = load(filename2(i,j)+'_phi.out');
        r    = data(:,1);
        phirb(i,j) = data(N1_(i)+1,2);
    end
end
    err = abs((phirb - phiAna)/phiAna)*100; %est une matrice
    
    %plot 
    pentes = zeros(1,np);
    for j=1:np
        loglog(N, err(:,j),'b x','Linewidth',lw,'HandleVisibility','off');
        hold on
        P = polyfit(log(N),log(err(:,j)),1); 
        z  = polyval(P, log(N));
        h = loglog(N,exp(z),'--','Linewidth',1);
        pentes(j) = P(1);
        %label(h,"p = " + num2str(p_(j)) + "; slope:" + num2str(P(1)),'location','right','slope','interpreter','latex','FontSize',fs-4)
    end
    legendStrings = "p =" + string(p_) + "; slope:" + string(pentes);
    leg = legend(legendStrings,'NumColumns',2,'Location','southeast','FontSize',fs-6);
    xlabel('$N$; $N_2 = 2N_1$'); ylabel('$\chi_b$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG("ConvergencePhiRBP"+num2str(MeshFactor_));

%% ordre en fonction de p
nsimul      = 10; np = 30;
N           = round(logspace(1,3,nsimul));
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = linspace(0,1,np);
propMesh_   = ~trivial_;

N1_ = N; N2_ = N;
if(propMesh_)
    N2_ = MeshFactor_*N;
end

filename2 = strings(nsimul,np);
for i = 1:nsimul
    for j = 1:np
    filename2(i,j)="N1_"+ num2str(N1_(i))+"N2_"+num2str(N2_(i)) + "p_" + num2str(p_(j));
    end
end

% Simulations
for i = 1:nsimul
    for j=1:np
        N1_loc = N1_(i);
        N2_loc = N2_(i);
        p_loc  = p_(j);
        writeConfigP;
        disp('Exercice6 configuration.in');   
        system('Exercice6 configuration.in');
    end
end


figure('Name',"Convergence phi RBP")
phiAna = V0_ - (a0_*b_^2)/(12*pi) * log(b_/R_);
phirb  = zeros(nsimul,np); %valeur de phi en r=b
for i=1:nsimul
    for j=1:np
        data = load(filename2(i,j)+'_phi.out');
        r    = data(:,1);
        phirb(i,j) = data(N1_(i)+1,2);
    end
end
    err = abs((phirb - phiAna)/phiAna)*100; %est une matrice
     
    ordres = zeros(1,np);
    for j=1:np
        P = polyfit(log(N),log(err(:,j)),1); 
        ordres(j) = -P(1);
    end
    plot(p_,ordres,'-','Linewidth',lw);
    xlabel('$p$'); ylabel('order of convergence');
    grid minor; hold on; set(gca,'fontsize',fs);
SaveIMG("OrdreConvergencePhiRBP"+num2str(MeshFactor_));