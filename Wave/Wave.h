#ifndef WAVE 
#define WAVE

#include <cmath>
#include <Eigen/Dense>

#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

class Wave
{
public:
	template <class Tbf>
	Wave(const Tbf & bf, double k, double centre, double sigma) : 
	m_k{k}, m_centre{centre}, m_sigma{sigma}, 
	m_responseCoef(1, bf.maxRadialIndex())
	{m_responseCoef(0) = densityResolving(bf, bf.maxRadius());} 


	~Wave() {}


	template <class Tbf>
	void density2dEvolution(const Tbf & bf, const std::string & filename, const int skip, const double rMax =10) const {m_responseCoef.write2dDensity2File(filename, bf, skip, rMax);}

	template <class Tbf>
	void density2dEvolution(int timeIndex, const Tbf & bf, const std::string & filename, const double rMax =10) const {m_responseCoef.write2dDensity2File(timeIndex, filename, bf, rMax);}
	void saveEvolutionCoeff(const std::string & filename) const {m_responseCoef.write2File(filename);}

	Eigen::VectorXcd operator()(int timeIndex) const {return m_responseCoef(timeIndex);}
	Eigen::VectorXcd& operator()(int timeIndex) {return m_responseCoef(timeIndex);}

	void resizeVector(const int nTimeStep) {m_responseCoef.resizeVector(nTimeStep);}
	double density(const double r, const double phi) const {return density2D(r, phi);}

	double density(double r) const {return 1/sqrt(2 * M_PI * m_sigma * m_sigma)*exp(-0.5 * pow((r - m_centre )/m_sigma,2));}

	template <class Tbf>
	void density2dFinal(const Tbf & bf, const std::string & outFilename, const int fromEnd = 1, const double rMax = 10, const double nStep = 201) const {m_responseCoef.write2dDensity2File(m_responseCoef.nTimeStep() - fromEnd, outFilename, bf, rMax, nStep);}

	void removeIC() {for (int time = 1; time < m_responseCoef.nTimeStep(); ++time) {m_responseCoef(time) -= m_responseCoef(0);}}

	void checkCoeff() const {for (int time = 0; time < m_responseCoef.nTimeStep(); ++time) {std::cout << m_responseCoef(time)(0) <<'\n';}}

private: 
	
	const double m_centre, m_sigma, m_k; 
	
	ExpansionCoeff m_responseCoef; 


	double density2D(const double r, const double phi, const double m = 2) const {return density(r) * cos(m*(phi+m_k*r));} 	
	Eigen::ArrayXXcd density2D(const int nGrid, const double rMax) const; 

	template <class Tbf>
	Eigen::VectorXcd densityResolving(const Tbf & bf, const double rMax, const int nGrid = 1000) const;
	
};


Eigen::ArrayXXcd Wave::density2D(const int nGrid, const double rMax) const {
	Eigen::ArrayXXcd density = Eigen::ArrayXXcd::Zero(nGrid, nGrid); 

	double x{}, y{}, theta{}, r{}, spacing{2*rMax/((double) nGrid-1)}, centre{0.5*(nGrid-1)};
	for (int i = 0; i < density.rows(); ++i)
	{for (int j = 0; j < density.cols(); ++j){
			x = spacing * (i - centre); y = spacing * (j - centre);
			r = sqrt(x*x+y*y); theta = atan2(y,x);
			if (r!= 0){density(i,j) = density2D(r, theta);}
		}
	}
	return density;
}

template <class Tbf>
Eigen::VectorXcd Wave::densityResolving(const Tbf & bf, const double rMax, const int nGrid) const {
	Eigen::ArrayXXcd density = density2D(nGrid, rMax);
	return bf.densityResolving(density, rMax);
}

#endif



