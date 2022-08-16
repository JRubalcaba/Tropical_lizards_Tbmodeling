#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

////////////////////////// Transient heat model

// [[Rcpp::export]]
double Tb_t_cpp(double M, double a, int ground, double S, double Ta, double Tsk, double Tg, double v, double T0, double time){
  
  /// Derive morphological traits from allometry
  
  double A = 0.0314 * 3.14 * pow(M * 1e-3, 0.6666667); // Surface area (m2), O'Connor 1999
  double Ad = 0.6666667 * A; // dorsal area
  double Ag = 0.3333333 * A; // ventral area
  double L = 3.3 * pow(M * 1e-6, 0.33333); // Charactersitic dimension (m) -> Mitchell 1996
  double cp = 3.7; // Specific heat capacity body (J g-1 ?C-1)
  double C = M * cp;
  
  /// Convective heat transfer coefficient
  
  double nu = -1.1555e-14*pow(Ta+273, 3) + 9.5728e-11*pow(Ta+273, 2) + 3.7604e-08*(Ta+273) - 3.4484e-06; // Thermal conductivity of air (W m-1 K-1)   
  double k = 1.5207e-11*pow(Ta+273, 3) - 4.8574e-08*pow(Ta+273, 2) + 1.0184e-04*(Ta+273) - 3.9333e-04; // Kinematic viscosity of air (m2 s-1) 
  
  double Re = v * L / nu;
  double Nu = 0.1 * pow(Re, 0.74);
  double hc = Nu * k / L; // Convection heat transfer coef (W m-2 ?C-1)
  double k_skin = 0.5; // 0.027; // Conductivity of the skin (W m-1 ?C-1)
  double t_skin = 0.025 * pow(0.001 * M / (3.14 * 1000), 0.2); // thickness of skin (Stevenson 1985)
  double p_contact = 0.1; // Proportion of ventral skin in direct contact with substrate 
  double hg = p_contact * k_skin / t_skin; // Conduction heat transfer coef (W m-2 ?C-1)
  
  if(ground == 0) hg = 0;
  double epsilon = 0.98; // emissivity IR
  double beta = 5.67e-8; // Stefan Boltzmann constant (W m-2 K-4)
  double Ra = 4 * epsilon * beta * pow(Tsk+273, 3); // Radiative heat transfer coef (W m-2 ?C-1)
  double Rg = 4 * epsilon * beta * pow(Tg+273, 3);
  
  /// Body temperature
  double j = Ad / C * (a * S + Ra * Tsk + hc * Ta) + Ag / C * Tg * (Rg + hg);
  double theta = Ad / C * (Ra + hc) + Ag / C * (Rg + hg);
  
  double Tb = j/theta + (T0 - j/theta) * exp(-theta * time);
  return Tb;
}

////////////////////////// Operative temperature model

// [[Rcpp::export]]
double Te_model_cpp(double M, double a, int ground, double S, double Ta, double Tsk, double Tg, double v){
  
  /// Derive morphological traits from allometry
  
  double A = 0.0314 * 3.14 * pow(M * 1e-3, 0.6666667); // Surface area (m2), O'Connor 1999
  double Ad = 0.6666667 * A; // dorsal area
  double Ag = 0.3333333 * A; // ventral area
  double L = 3.3 * pow(M * 1e-6, 0.33333); // Charactersitic dimension (m) -> Mitchell 1996
  double cp = 3.7; // Specific heat capacity body (J g-1 ?C-1)
  double C = M * cp;
  
  /// Convective heat transfer coefficient
  
  double nu = -1.1555e-14*pow(Ta+273, 3) + 9.5728e-11*pow(Ta+273, 2) + 3.7604e-08*(Ta+273) - 3.4484e-06; // Thermal conductivity of air (W m-1 K-1)   
  double k = 1.5207e-11*pow(Ta+273, 3) - 4.8574e-08*pow(Ta+273, 2) + 1.0184e-04*(Ta+273) - 3.9333e-04; // Kinematic viscosity of air (m2 s-1) 
  
  double Re = v * L / nu;
  double Nu = 0.1 * pow(Re, 0.74);
  double hc = Nu * k / L; // Convection heat transfer coef (W m-2 ?C-1)
  double k_skin = 0.5; // Conductivity of the skin (W m-1 ?C-1)
  double t_skin = 0.025 * pow(0.001 * M / (3.14 * 1000), 0.2); // thickness of skin (Stevenson 1985)
  double p_contact = 0.1; // Proportion of ventral skin in direct contact with substrate 
  double hg = p_contact * k_skin / t_skin; // Conduction heat transfer coef (W m-2 ?C-1)
  
  if(ground == 0) hg = 0;
  double epsilon = 0.98; // emissivity IR
  double beta = 5.67e-8; // Stefan Boltzmann constant (W m-2 K-4)
  double Ra = 4 * epsilon * beta * pow(Tsk+273, 3); // Radiative heat transfer coef (W m-2 ?C-1)
  double Rg = 4 * epsilon * beta * pow(Tg+273, 3);
  
  /// Body temperature
  double j = Ad / C * (a * S + Ra * Tsk + hc * Ta) + Ag / C * Tg * (Rg + hg);
  double theta = Ad / C * (Ra + hc) + Ag / C * (Rg + hg);
  
  double Te = j/theta;
  return Te;
}

////////////////////////// Behavioral thermroegulation model

// [[Rcpp::export]]
NumericMatrix behav_therm(double M, double a, int ground, 
                          NumericVector S_sun, NumericVector Ta_sun, NumericVector Tsk_sun, NumericVector Tg_sun, NumericVector V_sun, 
                          NumericVector S_shade, NumericVector Ta_shade, NumericVector Tsk_shade, NumericVector Tg_shade, NumericVector V_shade,
                          double Tpref_mean, double Tpref_sd, double sigma_sun, double sigma_shade,
                          int max_iter){
  
  NumericMatrix Tb(1440, max_iter);
  for(int i=0; i<max_iter; i++){
    Tb(0,i) = 20; // initial Tb
  }
  NumericMatrix beta(1440, max_iter);
  for(int i=0; i<max_iter; i++){
    beta(0,i) = 0; //initial microenv (shade)
  }

  for(int iter = 0; iter < max_iter; iter++){
    
    for(int i=1; i < 1440; i++){
      
      double Tpref = R::rnorm(Tpref_mean, Tpref_sd);
      
      if(Tb(i-1,iter) < Tpref){
        if(sigma_sun > R::runif(0,1)){
          beta(i-1,iter) = 1; // if Tb < Tpref, move to the sun (beta 1) with prob sigma sun
        }
      }else{
        if(sigma_shade > R::runif(0,1)){
          beta(i-1,iter) = 0; // else, find the shade with prob sigma shade
        }
      }
      
      if(beta(i-1,iter) == 1){ // if beta 1 (sun)
        Tb(i,iter) = Tb_t_cpp(M, a, ground=ground, S_sun[i-1], Ta_sun[i-1], Tsk_sun[i-1], Tg_sun[i-1], V_sun[i-1], Tb(i-1,iter), 60);
      }else{
        Tb(i,iter) = Tb_t_cpp(M, a, ground=ground, S_shade[i-1], Ta_shade[i-1], Tsk_shade[i-1], Tg_shade[i-1], V_shade[i-1], Tb(i-1,iter), 60);
      }
      
      beta(i,iter) = beta(i-1,iter);
    }
  }

  NumericMatrix output(1440, 4); // column 0 is Tb; column 1 is beta
  for(int i=0; i < 1440; i++){
    output(i,0) = mean(Tb(i,_));
    output(i,1) = mean(beta(i,_));
    output(i,2) = sd(Tb(i,_));
    output(i,3) = sd(beta(i,_));
  }
  
  return output; 
}

