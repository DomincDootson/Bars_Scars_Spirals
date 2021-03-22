#ifndef BODIES
#define BODIES

#include <string>
#include <valarray>
#include <Eigen/Dense>

#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"
#include "../../DF_Class/Mestel.h"
#include "../../DF_Class/DFClass.h"

class Bodies {
public:
  Bodies(const int N) : m(N), xy(2*N), vxvy(2*N), n(N) {}

  template <class Tdf>
  Bodies(std::string fileName, const Tdf & df, int numbParticles);

  ~Bodies() {}
 
  std::valarray<double> m, xy, vxvy;
  const int n;
  
  void samplingDF();

  // some public functions that do the outputting of the coefficents
  std::valarray<double> radius(double softening2 = 0) const;
  std::valarray<double>  angle() const ;

  template <class Tbf>
  Eigen::VectorXcd responseCoefficents(const Tbf & bf, const double xi) const;

private:
  void convert2Cartesian();
};

template <class Tdf>
Bodies::Bodies(std::string fileName, const Tdf & df, int numbParticles) : m(numbParticles), xy(2*numbParticles), vxvy(2*numbParticles), n(numbParticles)
{
  std::random_device generator;
  std::uniform_real_distribution<double> uniform(0,1);
  std::vector<double> cumulativeDensity = df.readInCumulativeDensity(fileName);
  for (int i = 0; i < 2*numbParticles; i += 2){
    xy[i]     = df.radiusSampling(uniform(generator), cumulativeDensity);
    xy[i+1]   = 2 * M_PI * uniform(generator);

    vxvy[i]   = df.vRSampling();
    std::cout << "Sampling vPhi\n";
    vxvy[i+1] = df.vPhiSampling(xy[i], vxvy[i]);
  }
  convert2Cartesian();
}


void Bodies::convert2Cartesian() {

  for (int i =0; i < 2*n; i += 2){
    double r{xy[i]}, theta{xy[i+1]}, vr{vxvy[i]}, vPhi{vxvy[i+1]};
    xy[i]     = r * cos(theta); 
    xy[i+1]   = r * sin(theta);
    vxvy[i]   = cos(theta)*vr - sin(theta)*vPhi;
    vxvy[i+1] = sin(theta)*vr + cos(theta)*vPhi;
  }
}



std::valarray<double> Bodies::radius(double softening2) const 
{
  std::valarray<double> radius(n);
  
  for (int i = 0; i<n; i++)
  {
    radius[i] = sqrt(xy[2*i]*xy[2*i] +xy[2*i+1]*xy[2*i+1] +softening2);
  }

  return radius;
}




std::valarray<double>  Bodies::angle() const // Please for the love of god tidy this up
{

  std::valarray<double> angles(n);
  double x, y;

  for (int i = 0; i < n; ++i)
  {
    x = xy[2*i];
    y = xy[2*i+1];
    if ((x>0) && (y>0)){
      angles[i] = atan(y/x);
    }
    else if ((x<0) && (y>0)){
      angles[i] = M_PI - atan(abs(y/x));
    }
    else if ((x<0) && (y<0)){
      angles[i] = M_PI + atan(abs(y/x));
    }
    else if (x == 0 && y<0){
      angles[i] = 1.5*M_PI;
    }
    else if (x==0 && y>0){
      angles[i] = 0.5 * M_PI; 
    }
    else if (y==0 && x<0){
      angles[i] = M_PI;
    }
    else if (x==0 && y==0){
      angles[i] = 0;
    }
    else {
      angles[i] = 2 * M_PI - atan(abs(y/x));
    }
  }
  return angles;
}

template <class Tbf>
Eigen::VectorXcd Bodies::responseCoefficents(const Tbf & bf, const double xi) const{
  
  Eigen::VectorXcd coef = Eigen::VectorXcd::Zero(bf.maxRadialIndex()+1);
  std::valarray<double> rad{radius(0)}, ang{angle()};

  for (int i = 0; i < n; ++i){
    for (int j = 0; j < coef.size(); ++j){
      coef[j] += exp(-1i*(ang[i]*(bf.fourierHarmonic())))*bf.potential(rad[i], j);
    }
  }

  return - xi *(bf.scriptE()).inverse() * coef;
}

#endif 