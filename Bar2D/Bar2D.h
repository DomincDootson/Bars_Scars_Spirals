#ifndef BAR2D
#define BAR2D 
//For sake of argument we assume that it is the positive harmoic

#include <Eigen/Dense>
#include <complex>
#include <iostream>


// Pass it the basis functions ?? 

class Bar2D
{
public:
	Bar2D() :
	m_theta, m_omega, m_alpha{0}, m_momentOfInertia{momentOfInertia()}
	m_scriptE{basisFunctionContainer.getScriptE()};
	{}
	~Bar2D() {} 

	double torque(const Eigen::VectorXcd &diskCoeff);
	double momentOfInertia(); // We porbably need to pass it the basis functions

	void drift(const double timeStep);
	void kick(const double timeStep);

	void saveBarEvolution(ofstream &outPut);

private: 
	double m_theta, m_omega, m_alpha, m_momentOfInertia;
	const int m_fourierHarmonic;  
	const Eigen::MatrixXd m_scriptE; 

	Eigen::VectorXcd m_expansionCoeff
};

void Bar2D::saveBarEvolution(ofstream &outPut){
	outPut << m_theta, ',' << m_omega << ',' << m_alpha << '\n';
}

double Bar2D::torque(const Eigen::VectorXcd &diskCoeff) // Torque = il * (A.scriptE.conj(B) - b.scriptE.conj(A))
{
	complex<double>    firstTerm{ m_expansionCoeff.dot(m_scriptE *  diskCoeff.conjugate()) };
	complex<double> secondTerm{ m_diskCoeff.dot(m_scriptE *  m_expansionCoeff.conjugate()) }; 
	complex<double> unitComplex(0,1);
	return unitComplex * m_fourierHarmonic (firstTerm - secondTerm); 

}

double momentOfInertia(){

}

void Bar2D::drift(const double timeStep){
	complex<double> unitComplex(0,1);
	double deltaTheta{m_omega*timeStep + 0.5 * m_alpha * timeStep * timeStep}; 
	
	m_expansionCoeff = exp(-unitComplex * deltaTheta * m_fourierHarmonic);
	m_theta += deltaTheta; 
}

void Bar2D::kick(const double timeStep);


#endif