//
// Created by mikhail on 10/18/21.
//

#include <complex>
#include <iostream>
#include <random>
#include <vector>

#include "TProfile.h"

int GenerateMultiplicity( int mean_multiplicity );
double GenerateParticle( std::vector<double>  );
bool IsAccepted( double phi, const std::vector<double>& phi_all, double phi_resolution=0.05 );

int main(){
  using namespace std::complex_literals;
  int N_events = 100'000;
  std::complex<double> psi_rp = std::exp<double>(0.0);
  double phi_resolution=0.05;
  int N_generated = 0;
  int N_registered = 0;
  double rec_mean_cos=0.0;
  double rec_glob_sum_w=0.0;
  double tru_mean_cos=0.0;
  double tru_glob_sum_w=0.0;
  // ********************* Protons *********************
  auto v1 = 0.2;
  int proton_mean_multiplicity = 1;
  // ******************** All particles *********************
  double v1_all = 0.15;
  int all_mean_multiplicity = 15;
  // ******************** Loop through event *********************
  for( int i=0; i< N_events; ++i ) {
    // ******************** Generation of background particles *********************
    int all_multiplicity = GenerateMultiplicity( all_mean_multiplicity );
    std::vector<double> all_particles_phi;
    for( int n=0; n<all_multiplicity; ++n ){
      auto phi = GenerateParticle({v1_all});
      all_particles_phi.push_back(phi);
    }
    // ******************** Generation of particles of interest *********************
    int proton_multiplicity = GenerateMultiplicity( proton_mean_multiplicity );
    std::vector<double> protons_phi;
    std::complex<double> rec_q{0, 0};
    std::complex<double> tru_q{0, 0};

    double rec_sum_w=0.0;
    double tru_sum_w=0.0;

    N_generated+=proton_multiplicity;
    for( int n=0; n< proton_multiplicity; ++n ) {
      auto phi = GenerateParticle({v1});
      tru_q+=std::exp<double>(1i*phi );
      tru_sum_w+=1;
      if( !IsAccepted( phi, all_particles_phi, phi_resolution ) )
        continue;
      rec_q+=std::exp<double>(1i*phi );
      rec_sum_w+=1;
      N_registered++;
    }
    // ***********************************************************
    rec_q/=fabs(rec_sum_w) > std::numeric_limits<double>::min() ? rec_sum_w : 1.0;
    rec_mean_cos+=std::real(rec_q*conj(psi_rp));
    rec_glob_sum_w+=fabs(rec_sum_w) > std::numeric_limits<double>::min() ? 1.0 : 0.0;
    // ***********************************************************
    tru_q/=fabs(tru_sum_w) > std::numeric_limits<double>::min() ? tru_sum_w : 1.0;
    tru_mean_cos+=std::real(tru_q*conj(psi_rp));
    tru_glob_sum_w+=fabs(tru_sum_w) > std::numeric_limits<double>::min() ? 1.0 : 0.0;
  }
  // Directly calculated
  double v1_rec = rec_mean_cos/rec_glob_sum_w;
  double v1_tru = tru_mean_cos/tru_glob_sum_w;
  double efficiency = (double) N_registered / N_generated;
  // Directly calculated
  std::cout << "Calculated" << std::endl;
  std::cout << v1_tru << " " << v1_rec << " " << efficiency << std::endl;
  std::cout <<  std::endl;
  // Directly calculated
//  double expected_efficiency = 1.0 - all_mean_multiplicity * ( phi_resolution / (2*M_PI) + 1.0/M_PI * v1_all*v1*sin(phi_resolution) ) ;
  double expected_efficiency = 1.0 - all_mean_multiplicity * ( phi_resolution / (2*M_PI)  ) ;
  std::cout << "Expected" << std::endl;
  std::cout << expected_efficiency << std::endl;

  return 0;
}

int GenerateMultiplicity(int mean_multiplicity) {
  std::random_device random_device;
  std::default_random_engine generator(random_device());
  std::poisson_distribution<int> distribution(mean_multiplicity);
  int number=0;
  while( number ==0 )
    number = distribution(generator);
  return number;
}

double GenerateParticle( std::vector<double> harmonoics ){
  double phi;
  std::random_device random_device;
  std::default_random_engine random_engine(random_device());
  std::uniform_real_distribution<double> phi_dist(-M_PI, M_PI);
  std::uniform_real_distribution<double> p_dist(0, 1);
  phi = phi_dist(random_engine);
  double prob = 1.0 / (2.0*M_PI);
  int n=1;
  for( auto vn : harmonoics ) {
    prob += 2*vn *cos( n*phi ) / (2.0*M_PI) ;
    n++;
  }
  auto p = p_dist(random_engine);
  if( p < prob )
    return phi;
  else
    return GenerateParticle( harmonoics );
}

bool IsAccepted( double phi, const std::vector<double>& phi_all, double phi_resolution ){
  std::random_device random_device;
  std::default_random_engine random_engine(random_device());
  std::uniform_real_distribution<double> p_dist(0, 1);
  for( auto bg : phi_all ){
    if( fabs(bg-phi) < phi_resolution ) {
      auto p = p_dist(random_engine);
      if( p<0.5 )
        return false;
    }
  }
  return true;
}
