#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "ConfigFile.h"

using namespace std;
double const pi=3.14159265358979323846264338327950288419716939937510582097494459230e0;
// Resolution d'un systeme d'equations lineaires par elimination de Gauss-Jordan:
template <class T>
vector<T> solve(vector<T> const& diag,
                vector<T> const& lower,
                vector<T> const& upper,
                vector<T> const& rhs)
{
  vector<T> solution(diag.size());
  vector<T> new_diag(diag);
  vector<T> new_rhs(rhs);

  for(unsigned int i(1); i<diag.size(); ++i)
  {
    double pivot = lower[i-1]/new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  for(int i(diag.size()-2); i>=0; --i)
    solution[i] = (new_rhs[i] - upper[i]*solution[i+1]) / new_diag[i];

  return solution;
}

bool if_left(double const& r,double const& b)
{
	return r < b;
}

// TODO: Classe pour epsilon_r(r)
class Epsilonr {
public:
  Epsilonr(bool const& trivial_, double const& b_, double const& c_, \
  double const& relative_epsilon_)
    : trivial(trivial_), b(b_), R(c_),  relative_epsilon(relative_epsilon_) {};

  inline double operator()(double const& r, bool const& left) {
  // Le booleen "left" indique s'il faut prendre la limite a gauche ou a droite en cas de discontinuite
double eps(1e-12*b);
if(trivial or r<b-eps or (abs(r-b)<=eps and left))
  return 1.0;
else
  return relative_epsilon;
}

private:
  bool trivial=1;
  double b=0., R=0., relative_epsilon=0.;

};

// TODO: Classe pour rho_lib(r)/epsilon_0
class Rho_lib {
public:
  Rho_lib(bool const& trivial_, double const& b_, double const& a0_,double const& R_)
    : b(b_), a0(a0_), R(R_), trivial(trivial_) {};

  inline double operator()(double const& r) {
    if(trivial)
      return 1.0;
    else if(r<b and !trivial)
      return a0*sin(3*pi*r/b);
    else
      return 0;
  }
// pourquoi avoir R en attribut ?
private:
  double b, a0, R;
  bool trivial;
};

int main(int argc, char* argv[])
{

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice6 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice6 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Domaine :
  const double b(configFile.get<double>("b"));
  const double R(configFile.get<double>("R"));

  // Conditions aux bords :
  const double V0(configFile.get<double>("V0"));

  // Instanciation des objets :
  Epsilonr epsilonr(configFile.get<bool>("trivial"), b, R,configFile.get<double>("epsilon_r"));
  Rho_lib rho_lib(configFile.get<bool>("trivial"), b, configFile.get<double>("a0"),R);

  // Discretisation du domaine :
  int N1 = configFile.get<int>("N1");
  int N2;
  if(configFile.get<bool>("proportionalMesh")){
    // apply proportion
    N2 = (int)(configFile.get<double>("meshFactor")*((double)N1));
  } else {
    N2 = configFile.get<int>("N2");
  }
  //Show the actual value of N2
  //cout << " N2=" << N2 << endl;

  //
  // TODO: Intégration mixte : p*trapeze+(1-p)*mid-point, voir Eq.(3.48) des notes de cours.
  // Pour la règle des trapèzes vous aurez donc p=1.
  // L'implémentation du schéma d'intégration mixte est facultative.
  double p = configFile.get<double>("p");
  int ninters = N1 + N2;
  int npoints = N1 + N2 + 1;

  /*vector<double> r(npoints,0);
  vector<double> h(ninters,b/(N1));
  for(size_t k(N1); k<h.size();++k)
  {
	  h[k] = (R-b)/(N2);
  }
  for(size_t k(1); k<r.size(); ++k)
  {
	  //r[k] = r[k-1] + h[k-1];
    r[k] = r[k-1] + h[k];
  }*/

  double h1 = b/N1;
  double h2 = (R-b)/N2;

  //Positions
  vector<double> r(npoints);
  for(int i(0); i<N1; ++i){
    r[i] = i*h1;
  }
  for(int i(0); i<=N2; ++i){
    r[N1+i] = b + i*h2;
  }

  //Ecartement entre deux cases
  vector<double> h(ninters);
  h[0] = b/N1;
  for(int i(1); i<ninters; ++i){
    if(i<N1)
    h[i] = b/N1;
    else
    h[i] = (R-b)/N2;
  }
  /*for(int i(1); i<ninters; ++i){
    h[i] = r[i] - r[i-1];
  }*/

  vector<double> diag(npoints,0.);  // Diagonale
  vector<double> lower(ninters,0.); // Diagonale inferieure
  vector<double> upper(ninters,0.); // Diagonale superieure
  vector<double> rhs(npoints,0.);   // Membre de droite

  // TODO: Assemblage des elements de la matrice et du membre de droite
  // Boucle sur les intervalles (k)
  for(int k(0); k<ninters; ++k)
  {
  // trapezoidal
	double trap_kk = p*h[k]*(r[k]*epsilonr(r[k],r[k]<b)+r[k+1]*epsilonr(r[k+1],r[k+1]<b))/(2*pow(h[k],2));
	//double trap_kplus1kplus1 = p*h[k]*(r[k]*epsilonr(r[k],r[k]<b)+r[k+1]*epsilonr(r[k+1],r[k+1]<b))/(2*pow(h[k+1],2));
	double trap_kkplus1 = - p*h[k]*(r[k]*epsilonr(r[k],r[k]<b)+r[k+1]*epsilonr(r[k+1],r[k+1]<b))/(2*h[k]*h[k]);
    // mid-point
	double mid_kk = (1-p)*h[k]*(r[k]+r[k+1])*epsilonr((r[k]+r[k+1])/2,(r[k]+r[k+1])/2<b)/(2*pow(h[k],2));
	//double mid_kplus1kplus1 = (1-p)*h[k]*(r[k]+r[k+1])*epsilonr((r[k]+r[k+1])/2,(r[k]+r[k+1])/2<b)/(2*pow(h[k+1],2));
	double mid_kkplus1 = -(1-p)*h[k]*(r[k]+r[k+1])*epsilonr((r[k]+r[k+1])/2,(r[k]+r[k+1])/2<b)/(2*h[k]*h[k]);
    // Be careful, epsilonr has jumps at r=b, which is a grid point. Depending on

    diag[k] += trap_kk + mid_kk;
    diag[k+1] += trap_kk + mid_kk;
    //diag[k+1] += trap_kplus1kplus1 + mid_kplus1kplus1;
    lower[k] += trap_kkplus1 + mid_kkplus1;
    upper[k] += lower[k];

    rhs[k] += h[k]*(p*r[k]*rho_lib(r[k])/2+(1-p)*(r[k] + r[k+1])/4*rho_lib((r[k]+r[k+1])/2));
    rhs[k+1] += h[k]*(p*r[k]*rho_lib(r[k+1])/2+(1-p)*(r[k] + r[k+1])/4*rho_lib((r[k]+r[k+1])/2));

  /*  double trap = 0.5*p*(epsilonr(r[k],false)*r[k] + epsilonr(r[k+1],true)*r[k+1])/h[k];
    double mid  = 0.5*(1-p)*epsilonr(0.5*(r[k+1]+r[k]),false)*(r[k+1]+r[k])/h[k];

    // round-off errors we might fall on the 'wrong' side.
    // Call the See the boolean 'left'

    diag[k] += trap + mid;
    upper[k]-= trap + mid;
    lower[k]-= trap + mid;

    diag[k+1] += trap + mid;

    rhs[k]   += 0.5*p*rho_lib(r[k])*r[k]*h[k] + 0.25*(1-p)*rho_lib(0.5*(r[k]+r[k+1]))*(r[k]+r[k+1])*h[k];
  	rhs[k+1] += 0.5*p*rho_lib(r[k+1])*r[k+1]*h[k] + 0.25*(1-p)*rho_lib(0.5*(r[k]+r[k+1]))*(r[k]+r[k+1])*h[k];*/
  }

  // TODO: Condition au bord:
   diag.back()  = 1.0;
   rhs.back()   = V0;
   lower.back() = 0.0;


  // Resolution:
  vector<double> phi(solve(diag,lower,upper,rhs));

  // Export des resultats:
  // 1. epsilon_r, rho_lib
  ofstream ofs(configFile.get<string>("outputDataProfiles").c_str());
  ofs.precision(15);
  for(int i(0); i<npoints; ++i)
    ofs << r[i] << " " << epsilonr(r[i],false) << " " << rho_lib(r[i]) << endl;
  ofs.close();

  // 2. phi
  ofs.open(configFile.get<string>("outputPotential").c_str());
  for(int i(0); i<npoints; ++i)
    ofs << r[i] << " " << phi[i] << endl;
  ofs.close();

  // 3. E_r et D_r/epsilon_0
  // N.B.: On définira Dr comme étant D_r/epsilon_0
  vector<double> rmid(ninters);
    for(int k(0); k<ninters; ++k)
  {
	  rmid[k] = (r[k] + r[k+1])/2;
  }
  vector<double> Er(ninters);
  vector<double> Dr(ninters);
  for(int i(0); i<ninters; ++i)
  {
    // TODO: Calculer E_r et D_r/epsilon_0 aux milieux des intervalles
    // en utilisant la représentation en éléments finis
    Er[i] = -(phi[i+1] - phi[i])/(r[i+1]-r[i]);
    Dr[i] = Er[i]*(epsilonr(rmid[i],(rmid[i]<=b)));
  }
  ofs.open(\
  configFile.get<string>("outputElectricDisplacementFields").c_str());
  ofs.precision(15);
  for(int i(0); i<ninters; ++i)
    ofs << rmid[i] << " " << Er[i] << " " << Dr[i] << endl;
  ofs.close();

  // 4. rho_lib, div(E_r) et div(D_r)/epsilon_0
  vector<double> rmidmid(ninters-1);
  for(int k(0); k<ninters-1; ++k)
  {
	  rmidmid[k] = (rmid[k] + rmid[k+1])/2;
  }
  vector<double> div_Er(ninters-1);
  vector<double> div_Dr(ninters-1);
  for(int i(0); i < ninters-1; ++i) //à vérifier !!
  {
    // TODO: Calculer div(E_r) et div(D_r)/epsilon_0
    // en utilisant des différences finies centrées aux milieux des milieux des intervalles
	   div_Er[i] = 1.0 / rmidmid[i] * ( ( rmid[i+1] * Er[i+1] - rmid[i] * Er[i] ) / ( rmid[i+1] - rmid[i] ) );
     div_Dr[i] = 1.0 / rmidmid[i] * ( ( rmid[i+1] * Dr[i+1] - rmid[i] * Dr[i] ) / ( rmid[i+1] - rmid[i] ) );
  }
  ofs.open(configFile.get<string>("outputDivergences").c_str());
  ofs.precision(15);
  for(int i(0); i<ninters-1; ++i)
    ofs << rmidmid[i] << " " << rho_lib(rmidmid[i]) << " " << div_Er[i] << " " << div_Dr[i] << endl;
  ofs.close();

  return 0;
}
