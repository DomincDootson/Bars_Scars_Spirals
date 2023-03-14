#ifndef SPIRAL2D
#define SPIRAL2D

#include <cmath>
#include <Eigen/Dense>

#include "../Potential_Density_Pair_Classes/KalnajsNBasis.h"
#include "../Volterra_Solver/ExpansionCoeff.h"

class Spiral2D
{
public:
	template <class Tbf>
	Spiral2D(const Tbf & bf, double k = 1, double scale = 5, double maxDensity = 1, int taper = 6) : 
	m_k{k}, m_scale{scale}, m_maxDensity{maxDensity}, m_taper{taper},
	m_responseCoef(1, bf.maxRadialIndex())
	{}//m_responseCoef(0) = densityResolving(bf, 15);} // This number is the upper limit of the basis function

	

	~Spiral2D() {}


	template <class Tbf>
	void density2dEvolution(const Tbf & bf, const std::string & filename, const int skip, const double rMax =10) const {m_responseCoef.write2dDensity2File(filename, bf, skip, rMax);}

	template <class Tbf>
	void density2dEvolution(int timeIndex, const Tbf & bf, const std::string & filename, const double rMax =10) const {m_responseCoef.write2dDensity2File(timeIndex, filename, bf, rMax);}
	void saveEvolutionCoeff(const std::string & filename) const {m_responseCoef.write2File(filename);}

	void susceptibilityEvolution(const std::string & filename) const;

	template <class Tbf>
	std::vector<double> densityPower(const Tbf & bf) const;


	Eigen::VectorXcd operator()(int timeIndex) const {return m_responseCoef(timeIndex);}
	Eigen::VectorXcd& operator()(int timeIndex) {return m_responseCoef(timeIndex);}

	void resizeVector(const int nTimeStep) {m_responseCoef.resizeVector(nTimeStep);}
	double density(const double r, const double phi) const {return density2D(r, phi);}
	double density(const double r) const {return (2 * m_maxDensity * pow(r/m_scale, 0.5*m_taper)) / (1+pow(r/m_scale, m_taper)); } 

	//double density(double r) const {return exp(-0.5 * pow((radius - )))}

	template <class Tbf>
	void density2dFinal(const Tbf & bf, const std::string & outFilename, const int fromEnd = 1, const double rMax = 10, const double nStep = 201) const {m_responseCoef.write2dDensity2File(m_responseCoef.nTimeStep() - fromEnd, outFilename, bf, rMax, nStep);}
	

	double maxDensity() const {return m_maxDensity;}
	void removeIC() {for (int time = 0; time < m_responseCoef.nTimeStep(); ++time) {m_responseCoef(time) -= m_responseCoef(0);}}

	template <class Tbf>
	void gaussianSpiral(const Tbf & bf, const double k, const double centre, const double width, const double rMax = 15); 

private: 
	const double m_k, m_scale, m_maxDensity;
	const int m_taper; 
	
	ExpansionCoeff m_responseCoef; 

	
	/* Functions defining Spiral Shape */ 
	/*-------------------------------- */ 

	//double density2D(const double r, const double phi, const double m = 2) const {return density(r) * cos(m*(phi-spiral(r)));} 

	double density2D(const double r, const double phi, const double m = 2) const {return density(r) * cos(m*(phi+m_k*r));} 
	double spiral(const double r, const double a = 1) const {return (1/m_k) * log(r/a);} 

	Eigen::ArrayXXcd density2D(const int nGrid, const double rMax) const; 

	template <class Tbf>
	Eigen::VectorXcd densityResolving(const Tbf & bf, const double rMax, const int nGrid = 500) const;
	
};


// We need to add a function to the Volterra Solve the does the evolution of the coefficents, i.e. pass it a spiral class
// A function that calculates the suseptibilty 


Eigen::ArrayXXcd Spiral2D::density2D(const int nGrid, const double rMax) const {
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
Eigen::VectorXcd Spiral2D::densityResolving(const Tbf & bf, const double rMax, const int nGrid) const {
	Eigen::ArrayXXcd density = density2D(nGrid, rMax);
	return bf.densityResolving(density, rMax);
}

template <class Tbf> 
std::vector<double> Spiral2D::densityPower(const Tbf & bf) const {
	int nStep{400}; double rMax{20};
	std::vector<double> power(m_responseCoef.nTimeStep()); 

	for (int i =0; i <power.size(); ++i) {
		if (i%20 ==0) {std::cout << "Outputting density power for time step: " << i << '\n';}	
		Eigen::ArrayXXcd den = bf.densityArray(m_responseCoef(i), nStep, rMax);
		power[i] = ((den * den.conjugate()).sum()).real(); // Don't include deltaR^2 as we normalise. 
	}
	double normalisation{power[0]};
	for (int i =0; i < power.size(); ++i) {power[i] = power[i]/normalisation;}

	return power; 
}

template <class Tbf>
void Spiral2D::gaussianSpiral(const Tbf & bf, const double k, const double centre, const double width, const double rMax ) {
	Eigen::VectorXcd coeff = Eigen::VectorXcd::Zero(bf.maxRadialIndex()+1);
	auto densityProfile = [centre, width] (double rad) {return 1/(width * sqrt(2*M_PI))*exp(-0.5 * pow((rad-centre)/width,2));};
	std::complex<double> unitComplex(0,1);
	double spacing =  0.01*width;
	for (int n = 0; n < bf.maxRadialIndex(); ++n) {
		for (double r = spacing; r < rMax; r+= spacing) {coeff(n) += r*r*bf.potential(r, n) * exp(unitComplex * r * k) * densityProfile(r);}
	}
	m_responseCoef(0) = -coeff*spacing; 
}


#endif
