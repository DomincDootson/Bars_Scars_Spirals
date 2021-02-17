#ifndef EXPANSIONCOEFF
#define EXPANSIONCOEFF 

#include <vector>
#include <cmath>

#include "../Potential_Density_Pair_Classes/PotentialDensityPairContainer.h"

class ExpansionCoeff
{
public:
	ExpansionCoeff(const int numbTimeStep, const int maxRadialIndex) : m_coeff(numbTimeStep) {
		for (auto i = m_coeff.begin(); i != m_coeff.end(); ++i){
			i -> setZero(maxRadialIndex+1);
		}
	}

		ExpansionCoeff(const Eigen::VectorXcd t0Coeff, const int numbTimeStep, const int maxRadialIndex) : m_coeff(numbTimeStep) {
		for (auto i = m_coeff.begin(); i != m_coeff.end(); ++i){			
			i -> setZero(maxRadialIndex+1);
		}
		m_coeff[0] = t0Coeff;
	}
	
	ExpansionCoeff(const std::string & filename, const int numbTimeStep) : m_coeff(numbTimeStep) 
	{coefficentReadIn(filename);}
	

	~ExpansionCoeff() {}

	
	
	Eigen::VectorXcd operator()(int timeIndex) const {return m_coeff[timeIndex];}
	Eigen::VectorXcd& operator()(int timeIndex) {return m_coeff[timeIndex];}
	
	void coefficentReadIn(const std::string &filename);

	void write2File(const std::string & filename, const int skip = 10) const;
	void writePerturbation2File(const std::string &filename) const;

	template <class Tbf>
	void writeDensity2File(const std::string & outFilename, const Tbf & bf, const int m_skip) const; 

private:
	
	std::vector<Eigen::VectorXcd> m_coeff;	
};
char sign(double number){
	if (number < 0){ return '-';}
	else {return '+';}
}


void ExpansionCoeff::coefficentReadIn(const std::string &filename)
{
	std::ifstream inFile(filename);
	int maxRadialIndex; inFile >> maxRadialIndex;
	assert(maxRadialIndex +1 == m_coeff[0].size() && "The perturbation does't have the same size as perturbation vector.");

	double real, imag;
	std::complex<double> unitComplex(0,1);
	for (auto it = m_coeff.begin(); it != m_coeff.end(); ++it){
	 	for (int n = 0; n < it -> size(); ++ n) {
	 		inFile >> real >> imag;
	 		(*it)(n) = real + unitComplex * imag;
	 	}
	 } 
	inFile.close();
}


void ExpansionCoeff::write2File(const std::string & filename, const int skip) const 
{
	std::ofstream out(filename);
	for (int time = 0; time < m_coeff.size(); time += skip)
	{
		for (int n = 0; n < m_coeff[time].size(); ++n)
		{
			if (n == m_coeff[time].size() -1){
				out << real(m_coeff[time](n)) << sign(imag(m_coeff[time](n)))  << abs(imag(m_coeff[time](n)))  << "j\n";
			}
			else{
			out << real(m_coeff[time](n)) << sign(imag(m_coeff[time](n)))  << abs(imag(m_coeff[time](n))) << "j,"; // Could we change the syntax t
			}
		}
	}
	std::cout << "Evolution saved to: " << filename << '\n';
	out.close();
}

void ExpansionCoeff::writePerturbation2File(const std::string &filename) const
{
	std::ofstream out(filename);
	out << m_coeff[0].size()-1 << '\n';
	for (int time = 0; time < m_coeff.size(); ++time)
	{
		for (int n = 0; n < m_coeff[time].size(); ++n)
		{
			if (n == m_coeff[time].size() -1){
				out << real(m_coeff[time](n)) << " " << imag(m_coeff[time](n))  << '\n';
			}
			else{
			out << real(m_coeff[time](n)) << " " << imag(m_coeff[time](n))  << " ";
			}
		}
	}
	std::cout << "Evolution saved to: " << filename << '\n';
	out.close();
}


void outputVector(std::ofstream & out, std::vector<double> vec){
	for (auto i = vec.begin(); i != vec.end() - 1; ++i){
		out << *i <<',';
	}
	out << vec.back() << '\n';
}

std::vector<double> radiiVector(double rMax, int nStep){
	double step{rMax/((double) nStep)};
	std::vector<double> radii;
	for (double radius = 0; radius <= rMax; radius += step){
		radii.push_back(radius);
	}
	return radii;
}

template <class Tbf>
void ExpansionCoeff::writeDensity2File(const std::string & outFilename, const Tbf & bf, const int skip) const 
{
	std::ofstream out(outFilename);
	std::vector<double> radii = radiiVector(10, 200);
	outputVector(out, radii);
	for (int time = 0; time < m_coeff.size(); time += skip){
		std::vector<double> densityOnLine = bf.oneDdensity(radii, m_coeff[time]);
		outputVector(out, densityOnLine);
	}
	out.close();
}


#endif