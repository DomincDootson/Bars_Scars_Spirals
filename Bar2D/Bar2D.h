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
	Bar2D(const Eigen::VectorXcd expansionCoeff, const double omega0, const int harmonic = 2) : // Note that omega0 is the pattern speed! 
	m_expansionCoeff{expansionCoeff},
	m_momentOfInertia{1},
	m_fourierHarmonic{harmonic},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(0); v_omega.push_back(omega0); v_alpha.push_back(0);
		std::cout << "Please note the bar is not evolving (Might cause seg fault).\n\n";}

	Bar2D(const Eigen::VectorXcd expansionCoeff, const double omega0, const std::string growthFilename, const int harmonic = 2) :
	m_expansionCoeff{expansionCoeff},
	m_momentOfInertia{1},
	m_fourierHarmonic{harmonic},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(0); v_omega.push_back(omega0); v_alpha.push_back(0); readInSize(growthFilename);}

	Bar2D(const double omega0, const std::string growthFilename) : 
	m_momentOfInertia{1},
	m_fourierHarmonic{2},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(0); v_omega.push_back(omega0); v_alpha.push_back(0); readInSize(growthFilename);}

	Bar2D(const double omega0) : 
	m_momentOfInertia{1},
	m_fourierHarmonic{2},
	m_scriptE{Eigen::MatrixXd::Identity(m_expansionCoeff.size(), m_expansionCoeff.size())}
	{v_theta.push_back(0); v_omega.push_back(omega0); v_alpha.push_back(0);}


	~Bar2D() {} 

	Eigen::VectorXcd barCoeff(const double time = -1) const {return barSize(time) * m_expansionCoeff;}
	double patternSpeed() const {return v_omega.back();}
	double angle() const {return v_theta.back();}
	
	double torque(const Eigen::VectorXcd &diskCoeff, const double time = -1) const;

	void drift(const double timeStep, const double freelyRotating);
	void kick(const double timeStep, const Eigen::VectorXcd &diskCoeff, const double freelyRotating, const double time = -1);

	void saveBarEvolution(const std::string & evolutionFilename, const int skip = 1) const;
	void readInSize(const std::string & filename);

	template <class Tbf>
	void sormaniBar(const Tbf & bf, const double size = 1, const std::string & filename = "Bar2D/Bar_Potentials/Sormani_Medium.out"); 
	double sormaniPotential(const double rad, const double time) const;
	
	double barSize(const double time) const;
private: 
	Eigen::VectorXcd m_expansionCoeff;

	std::vector<double> v_theta, v_omega, v_alpha; // We use these to keep track of the evolution of the bar
	std::vector<double> v_times, v_size; // Use these vectors to set the relative size of the bar 

	std::vector<double> v_sormaniRadii, v_sormaniPotential; 


	const double m_momentOfInertia;  
	const int m_fourierHarmonic;  
	const Eigen::MatrixXd m_scriptE; 

	

	int findTimeIndex(const double time) const;
	
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
	std::complex<double>    firstTerm{ m_expansionCoeff.dot(diskCoeff) }; // So the dot product of complex vectors
	std::complex<double> secondTerm{ diskCoeff.dot(m_expansionCoeff) };   // has the conjugation built in
	std::complex<double> unitComplex(0,1);
	double l{ (double) m_fourierHarmonic};
	//std::cout << "torque: " << barSize(time) * std::real(unitComplex * ((firstTerm - secondTerm) * l))<<'\n';
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

/* Sormani Stuff */ 
/* ------------- */ 

template <class Tbf>
void Bar2D::sormaniBar(const Tbf & bf, const double size, const std::string & filename) {
	std::ifstream inFile; inFile.open(filename);
	if (!inFile.good()) {std::cout << "Bar size file doesn't exist.\n"; exit(0);}
	int n; double step, rad{}, pot{}; 
	inFile >> n >> step;

	Eigen::VectorXcd integration = Eigen::VectorXcd::Zero(bf.maxRadialIndex()+1);
	m_expansionCoeff *= 0;
	
	for (int r = 0; r < n; ++r) {
		inFile >> rad >> pot;
		v_sormaniRadii.push_back(rad);
		v_sormaniPotential.push_back(size*pot); 
		for (int i = 0; i < m_expansionCoeff.size(); ++i) {
			m_expansionCoeff(i) += - size*step * (2*M_PI*rad)*bf.density(rad, i) * pot; 
		}
	}  
	//std::cout << m_expansionCoeff() << '\n'; 
	inFile.close();
}


double Bar2D::sormaniPotential(const double rad, const double time) const {
	double step{v_sormaniRadii[1]-v_sormaniRadii[0]}; 
	double index = (rad - v_sormaniRadii.front())/step;  
	if ((index < 0) || (index >= v_sormaniPotential.size()-2)) {return 0;}
	return  (ceil(index) - index) * (v_sormaniPotential[(int) floor(index)])+(index - floor(index)) * (v_sormaniPotential[(int) ceil(index)]) * barSize(time);
}



#endif