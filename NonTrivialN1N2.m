%% Paramètres: cas non trivial avec facteur 2 /////////////////////////////
ViewFormat;
Ntotal  = round(logspace(1,2,nsimul));
N1_ = MeshFactor_*(MeshFactor_ +1)*Ntotal;
N2_ = 1/(MeshFactor_ +1)* Ntotal;
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = 1.e0;
propMesh_   = ~trivial_;

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
    steps = Ntotal(1:intervalle:nsimul); legendStrings = string(steps);
    leg = legend(legendStrings,'Location','southeast','NumColumns',2);
    title(leg, "$N_{tot}, N_1 =$" + num2str(MeshFactor_) +"$N_2$")
% SaveIMG("PhiNonTrivial"+printname);
    
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
    steps = Ntotal(1:intervalle:nsimul); legendStrings = string(steps);
    leg = legend(legendStrings,'Location','northeast','NumColumns',2);
    title(leg, "$N_{tot}$, $N_1 =$" + num2str(MeshFactor_) +"$N_2$")
% SaveIMG("ENonTrivial"+printname);

%% Convergence de phi(r=b) pour différentes valeurs
Ntotal  = round(logspace(3,4,nsimul));
N1_ = MeshFactor_*(MeshFactor_ +1)*Ntotal;
N2_ = 1/(MeshFactor_ +1)* Ntotal;
trivial_    = false;
b_          = 3.e-1 ;
R_          = 5.e-1;
a0_         = -3.e4;
epsilon_r_  = 4.e0;
V0_         = 2.2e2 ;
p_          = 1.e0;
propMesh_   = ~trivial_;

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
    P = polyfit(1./Ntotal, phirb,1); yintercept = P(2); 
    err = abs((phirb - yintercept)/yintercept)*100;
    
    %plot 
    plot(1./Ntotal, phirb,'x','Linewidth',lw,'HandleVisibility','off');
    hold on
    abscisse = 1./Ntotal; legendStrings = string(P(2)); xrange = (0:1.e-3:max(abscisse));
    plot(xrange, P(1)*xrange + P(2),'--','Linewidth',lw,'HandleVisibility','off');
    plot(0,P(2),'r +','Markersize',15,'Linewidth',lw+2); %marquer le y-intercept
    leg = legend(legendStrings,'Location','northwest','NumColumns',2);
    title(leg, '$\phi_{as}$ [V]')
    xlabel("$1/N_{tot}$, $N_1 =$" + num2str(MeshFactor_) +"$N_2$"); ylabel('$\phi(r=b)$ [V]');
    grid minor; hold on; set(gca,'fontsize',fs);
% SaveIMG("PhiRBasymptotiqueMesh"+printname);

figure('Name',"Convergence phi r b")
    FitLOGLOG(1./Ntotal,err,1);
    hold on
    loglog(1./Ntotal,err,'x','Linewidth',lw,'HandleVisibility','off');
    xlabel("$1/N_{tot}$, $N_1 =$" + num2str(MeshFactor_) +"$N_2$"); ylabel('Rel. error on $\phi_{as}(r=b)$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
    legend('Location','northwest','Numcolumns',1);
% SaveIMG("ConvergencePhiRBasymptotiqueMesh"+printname);


