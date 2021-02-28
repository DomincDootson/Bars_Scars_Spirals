#ifndef BODIES
#define BODIES

#include <string>
#include <valarray>
#include <Eigen/Dense>

#include "../../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

class Bodies {
public:
  Bodies(const int N) : m(N), xy(2*N), vxvy(2*N), n(N) {}
  Bodies(std::string fileName, const int numbParticles);

  ~Bodies() {}
 
  std::valarray<double> m, xy, vxvy;
  const int n;
  
  // some public functions that do the outputting of the coefficents
  std::valarray<double> radius(double softening2 = 0) const;
  std::valarray<double>  angle() const ;

  template <class Tbf>
  Eigen::VectorXcd responseCoefficents(const Tbf & bf, const double xi) const;
};


Bodies::Bodies(std::string fileName, int numbParticles) : m(numbParticles), xy(2*numbParticles), vxvy(2*numbParticles), n(numbParticles)
{
  /*particle2D * holdingParticle;
  holdingParticle = new particle2D[n];

  particleSampling(holdingParticle, n, fileName); 

  double weight{1/((double) n)};
  for (int i = 0; i < n; ++i)
  {
    
    m[i] = weight;
    xy[2*i]   = holdingParticle[i].x;
    xy[2*i+1] = holdingParticle[i].y;

    vxvy[2*i]   = holdingParticle[i].px;
    vxvy[2*i+1] = holdingParticle[i].py;
  }

  delete [] holdingParticle;*/ 
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