#ifndef BAR2D
#define BAR2D 
//For sake of argument we assume that it is the positive harmoic

#include <Eigen/Dense>
#include <complex>
#include <iostream>

#include <../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h>

class Bar2D
{
public:
	Bar2D(const Eigen::VectorXcd expansionCoeff, const basisFunctionContainer &basisFunctions) :
	m_theta, m_omega, m_alpha{0},
	m_expansionCoeff{expansionCoeff};
	m_momentOfInertia{momentOfInertia(basisFunctions)}
	m_scriptE{basisFunction.getScriptE()};
	{}
	~Bar2D() {} 

	double torque(const Eigen::VectorXcd &diskCoeff) const;
	double momentOfInertia(const basisFunctionContainer &basisFunctions) const;

	void drift(const double timeStep);
	void kick(const double timeStep);

	void saveBarEvolution(ofstream &outPut);

	Eigen::MatrixXd barCoeff() const {return m_expansionCoeff;} // Maybe we could just overload = sign? 
private: 
	double m_theta, m_omega, m_alpha; // Note that alpha = Torque/MOI.
	Eigen::VectorXcd m_expansionCoeff;

	const double m_momentOfInertia;  
	const int m_fourierHarmonic;  
	const Eigen::MatrixXd m_scriptE; 

	

	Eigen::ArrayXXd rSquaredArray(const int nGrid, const double rMax) const;
};

void Bar2D::saveBarEvolution(ofstream &outPut) const {
	outPut << m_theta, ',' << m_omega << ',' << m_alpha << '\n';
}

double Bar2D::torque(const Eigen::VectorXcd &diskCoeff) const // Torque = il * (A.scriptE.conj(B) - b.scriptE.conj(A))
{
	std::complex<double>    firstTerm{ m_expansionCoeff.dot(m_scriptE *  diskCoeff.conjugate()) };
	std::complex<double> secondTerm{ m_diskCoeff.dot(m_scriptE *  m_expansionCoeff.conjugate()) }; 
	std::complex<double> unitComplex(0,1);
	return std::real(unitComplex * m_fourierHarmonic (firstTerm - secondTerm)); 

}

Eigen::ArrayXXd Bar2D::rSquaredArray(const int nGrid, const double rMax) const{ // Should all be fine, but a quick check won't hurt! 
	Eigen::ArrayXXd rSquaredArray{Eigen::ArrayXXd::Zero(nGrid, nGrid)};
	double x{}, y{}, rSquared{}, spacing{2*rMax/((double) nGrid-1)}, centre{0.5*(nGrid-1)};
	for (int i = 0; i < rSquaredArray.rows(); ++i){
		for (int j = 0; j < rSquaredArray.cols(); ++j){
			x = spacing * (i - centre); y = spacing * (j - centre); rSquared = x*x +y*y;
			rSquaredArray(i,j) = rSquared;
		}
	}
	return rSquaredArray;
}


double Bar2D::momentOfInertia(basisFunctionContainer &basisFunctions) const { // Ditto with this, a check wouldn't hurt 
	int nGrid{801}; // Parameters that set the size of the array. 
	double rMax{20};

	Eigen::ArrayXXd density{basisFunctions.densityArrayReal(m_expansionCoeff, nGrid, rMax)};
	Eigen::ArrayXXd rSquaredArray{rSquaredArray(nGrid, rMax)};
	return stepSize*stepSize*(density*rSquaredArray).sum();

}

void Bar2D::drift(const double timeStep){
	std::complex<double> unitComplex(0,1);
	double deltaTheta{m_omega*timeStep + 0.5 * m_alpha * timeStep * timeStep}; 
	
	m_expansionCoeff = exp(-unitComplex * deltaTheta * m_fourierHarmonic);
	m_theta += deltaTheta; 
}

void Bar2D::kick(const double timeStep, const Eigen::VectorXcd &diskCoeff){
	double oldAlpha = m_alpha;
	m_alpha = torque(diskCoeff)/m_momentOfInertia;
	return 0.5 * (oldAlpha + m_alpha) * timeStep; 
}


#endif