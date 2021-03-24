%va écrire le configuration_.in
filename = 'configuration_.in';
fid = fopen(filename,'wt'); 

%Paramètres physiques:
fprintf(fid,['trivial=', num2str(trivial_),'\n']); % si vrai, le cas trivial est simule
fprintf(fid,['b=',num2str(b_),'\n']); 	   % rayon de la partie vide
fprintf(fid,['R=',num2str(R_),'\n']);        % rayon exterieur
fprintf(fid,['a0=',num2str(a0_),'\n']);  % valeur de reference du profil de densite de charges libres
fprintf(fid,['epsilon_r=',num2str(epsilon_r_),'\n']); % valeur de reference du profil de la constante dielectrique
fprintf(fid,['V0=',num2str(V0_),'\n']);   % potentiel au bord

%Paramètres numériques:
fprintf(fid,['proportionalMesh=',num2str(propMesh_),'\n']);% si vrai, N2=meshFactor*N1
fprintf(fid,['meshFactor=',num2str(MeshFactor_),'\n']);% mesh proportional factor 
fprintf(fid,['N1=',num2str(N1_loc),'\n']);% nombre d'intervalles domain 1
fprintf(fid,['N2=',num2str(N2_loc),'\n']);% nombre d'intervalles domain 2
fprintf(fid,['p=',num2str(p_),'\n']);% methode d'integration: p=1 trapeze p=0 midpoint

%écriture dans des fichiers:
fprintf(fid,['outputDataProfiles= N1_',num2str(N1_loc),'N2_',+num2str(N2_loc),'_dataprof.out\n']); % ficher output de  r_i epsilon_r(r_i) et rho_lib(r_i)/epsilon_0
fprintf(fid,['outputPotential = N1_',num2str(N1_loc),'N2_',num2str(N2_loc),'_phi.out\n']); % ficher output de r_i et phi(r_i)
fprintf(fid,['outputElectricDisplacementFields = N1_',num2str(N1_loc),'N2_',num2str(N2_loc),'_E_D.out\n']); %ficher d'output de r_mid,i, E_r(r_mid,i) et D_r(r_rmid,i)/epsilon_0
fprintf(fid,['outputDivergences = N1_',num2str(N1_loc),'N2_',num2str(N2_loc),'_div_E_D.out\n']); % ficher output de r_midmid,i, div E(r_midmid,i) et div D(r_midmid,i)/epsilon_0 

fclose(fid);