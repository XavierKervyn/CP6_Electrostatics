%% Cas Trivial
clear all; clc
CasTrivial;

%% Cas Non trivial
clear all; clc
MeshFactor_ = 2;
CasNonTrivial;

%% Etude N1 vs N2
clear all; clc
nsimul  = 20;
Ntotal  = round(logspace(2,3,nsimul)); 
MF = 5;
printname = "N1=5N2";  MeshFactor_ = MF;   NonTrivialN1N2;
printname = "N1=02N2"; MeshFactor_ = 1/MF; NonTrivialN1N2;

% Etude de convergence avec paramètre t
Ntotal = 10000;
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
    disp('Exercice6 configuration_.in');   
    system('Exercice6 configuration_.in'); 
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

%% Question d.ii
clear; clc
nsimul = 10;
N = round(logspace(1,3,nsimul));
Questiond2;
