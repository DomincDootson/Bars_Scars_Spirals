#ifndef BAR2D
#define BAR2D 
// Note that this assumes the bar is l = 2. 

#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <fstream>

#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"


class Bar2D
{
public: 
	Bar2D(const Eigen::VectorXcd expansionCoeff, const double omega0) :
	m_expansionCoeff{expansionCoeff},
	m_momentOfInertia{1},
	m_fourierHarmonic{2},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(0); v_omega.push_back(omega0); v_alpha.push_back(0);
		std::cout << "Please note the bar is not evolving (Might cause seg fault).\n\n";}

	Bar2D(const Eigen::VectorXcd expansionCoeff, const double omega0, const std::string growthFilename) :
	m_expansionCoeff{expansionCoeff},
	m_momentOfInertia{1},
	m_fourierHarmonic{2},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(0); v_omega.push_back(omega0); v_alpha.push_back(0); readInSize(growthFilename);}



	~Bar2D() {} 

	Eigen::VectorXcd barCoeff(const double time = -1) const {return barSize(time) * m_expansionCoeff;}
	double patternSpeed() const {return v_omega.back();}
	double angle() const {return v_theta.back();}
	
	double torque(const Eigen::VectorXcd &diskCoeff, const double time = -1) const;

	void drift(const double timeStep, const double freelyRotating);
	void kick(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating, const double time = -1);

	void saveBarEvolution(const std::string & evolutionFilename, const int skip = 1) const;
	
	
private: 
	Eigen::VectorXcd m_expansionCoeff;

	std::vector<double> v_theta, v_omega, v_alpha; // We use these to keep track of the evolution of the bar
	std::vector<double> v_times, v_size; // Use these vectors to set the relative size of the bar 


	const double m_momentOfInertia;  
	const int m_fourierHarmonic;  
	const Eigen::MatrixXd m_scriptE; 

	

	int findTimeIndex(const double time) const;
	double barSize(const double time) const;
	void readInSize(const std::string & filename);
};

void Bar2D::saveBarEvolution(const std::string & evolutionFilename, const int skip) const {
	std::ofstream out(evolutionFilename);
	for (int time = 0; time < v_theta.size(); time += skip){
		out << v_theta[time] << ',' << v_omega[time] << ',' << v_alpha[time] << ',' << m_momentOfInertia * v_alpha[time] << '\n';
	}
	out.close();
}


double Bar2D::torque(const Eigen::VectorXcd &diskCoeff, const double time) const // Torque = il * (A.scriptE.conj(B) - b.scriptE.conj(A))
{
	std::complex<double>    firstTerm{ m_expansionCoeff.dot(m_scriptE *  diskCoeff) }; // So the dot product of complex vectors
	std::complex<double> secondTerm{ diskCoeff.dot(m_scriptE *  m_expansionCoeff) };   // has the conjugation built in
	std::complex<double> unitComplex(0,1);
	double l{ (double) m_fourierHarmonic};
	return -barSize(time) * std::real(unitComplex * ((firstTerm - secondTerm) * l)); // The negative sign comes from force  = - grad(phi)

}

void Bar2D::drift(const double timeStep, const double freelyRotating){
	std::complex<double> unitComplex(0,1); double theta{v_theta.back()};
	double deltaTheta{v_omega.back()*timeStep + freelyRotating * 0.5 * v_alpha.back() * timeStep * timeStep}; 

	m_expansionCoeff *= exp(-unitComplex * (deltaTheta * m_fourierHarmonic)); 
	theta += deltaTheta; 
	v_theta.push_back(theta);
}

void Bar2D::kick(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating, const double time){
	double oldAlpha{v_alpha.back()}, omega{v_omega.back()};

	double alpha = torque(diskCoeff, time)/m_momentOfInertia; 
	v_alpha.push_back(alpha); 
	omega += freelyRotating * 0.5 * (oldAlpha + v_alpha.back()) * timeStep;  
	v_omega.push_back(omega);
}

int Bar2D::findTimeIndex(const double time) const {
	if (time >= v_times.back()) {return v_times.size()-1;}
	else if (time == 0) {return 0;}

	int lowerIndex{0};
	while (v_times[lowerIndex+1] < time) { lowerIndex += 1; }
	return lowerIndex; 
}

double Bar2D::barSize(const double time) const {
	if ((time >= v_times.back()) || (time == -1)) {return 1;}
	
	int lowerIndex{findTimeIndex(time)};
	return v_size[lowerIndex] + 
	(v_size[lowerIndex+1]-v_size[lowerIndex]) *(time - v_times[lowerIndex])/(v_times[lowerIndex+1] - v_times[lowerIndex]);
}

void Bar2D::readInSize(const std::string & filename){
	std::ifstream inFile; inFile.open(filename);
	if (!inFile.good()) {std::cout << "Bar size file doesn't exist.\n"; exit(0);}
	int n; inFile >> n;
	for (int i = 0; i < n; ++i){ 
		double time, size;
		inFile >> time >> size;
		v_times.push_back(time); v_size.push_back(size);
	}
	inFile.close();
}


#endif