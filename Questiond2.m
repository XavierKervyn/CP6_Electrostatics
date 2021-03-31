%% Paramètres: cas non trivial et question d.ii ///////////////////////////
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

N1_ = N; N2 = N;
if (propMesh_)
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
    
% Question d.ii
figure('Name',"verification of Gauss' law")
for i=1:nsimul
    data   = load(filename2(i)+'_div_E_D.out');
    r      = data(:,1);
    rholib = data(:,2);
    divD   = data(:,4);
    
    plot(r, abs((rholib-divD)/a0_) * 100,'.-','Linewidth',lw);
    xlabel('$1/N$'); ylabel('Rel. err. on $A$ [\%]');
    grid minor; hold on; set(gca,'fontsize',fs);
end

figure('Name',"verification of Gauss' law avant b")
for i=1:nsimul
    data   = load(filename2(i)+'_div_E_D.out');
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
    data   = load(filename2(i)+'_div_E_D.out');
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
        xlabel('$r$ [m]'); ylabel('$\rho_{f} / \varepsilon_0 $ [V/m$^2$]');
        grid minor; hold on; set(gca,'fontsize',fs);