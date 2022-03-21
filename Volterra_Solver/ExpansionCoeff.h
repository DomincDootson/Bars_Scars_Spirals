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
	{coefficentReadInConstructor(filename);}
	

	~ExpansionCoeff() {}

	
	
	Eigen::VectorXcd operator()(int timeIndex) const {return m_coeff[timeIndex];}
	Eigen::VectorXcd& operator()(int timeIndex) {return m_coeff[timeIndex];}
	
	std::vector<double> energyEvolution(const Eigen::MatrixXd & scriptE);



	void coefficentReadIn(const std::string &filename);
	void coefficentReadInConstructor(const std::string &filename);

	void write2File(const std::string & filename, const int skip = 10) const;
	void writePerturbation2File(const std::string &filename) const;

	template <class Tbf>
	void writeDensity2File(const std::string & outFilename, const Tbf & bf, const int skip) const; 
	template <class Tbf>
	void write2dDensity2File(const std::string & outFilename, const Tbf & bf, const int skip, const double rMax=10, const int nStep=201) const; 
	template <class Tbf>
	void write2dDensity2File(int timeIndex, const std::string & outFilename, const Tbf & bf, const double rMax=10, const int nStep=201) const; 
	template <class Tbf>
	void write2dPotential2File(const std::string & outFilename, const Tbf & bf, const int skip) const;


	int nTimeStep() const {return m_coeff.size();}
	void resizeVector(const int nTimeStep); // Please note that this leaves the vecotor at index = 0 unchaged. 

private:
	
	std::vector<Eigen::VectorXcd> m_coeff;	
};

char sign(double number){
	if (number < 0){ return '-';}
	else {return '+';}
}

std::vector<double> ExpansionCoeff::energyEvolution(const Eigen::MatrixXd & scriptE)
{
	std::vector<double> energy;
	for (auto i = m_coeff.begin(); i != m_coeff.end(); ++i){
		double holding = -2 * real((*i).dot(scriptE * (*i)));
		energy.push_back(holding);
	}
	return energy;
}

void ExpansionCoeff::coefficentReadIn(const std::string &filename)
{
	std::ifstream inFile(filename);
	int maxRadialIndex; inFile >> maxRadialIndex;
	std::cout << "Reading perturbation in from: " << filename << '\n';
	std::cout << maxRadialIndex+1 << " " << m_coeff[0].size() << '\n';
	assert(maxRadialIndex +1 == m_coeff[0].size() && "The perturbation does't have the same size as perturbation vector."); // Note that the first row of the perturbation file must be nMax

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

void ExpansionCoeff::coefficentReadInConstructor(const std::string &filename)
{
	std::ifstream inFile(filename);
	int maxRadialIndex; inFile >> maxRadialIndex;
	std::cout << "Reading perturbation in from: " << filename << '\n';

	double real, imag;
	std::complex<double> unitComplex(0,1);
	int time{0};
	for (auto it = m_coeff.begin(); it != m_coeff.end(); ++it){
		it -> setZero(maxRadialIndex+1);
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
	for (int time = 0; time < m_coeff.size(); time += 1)
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

void outputArray(std::ofstream & out, Eigen::ArrayXXd array){
	for (int i = 0; i < array.cols(); ++i){
		for (int j =0; j <array.rows(); ++j){
			out << array(i,j) << ',';
		}
	}

	out << array(array.rows()-1, array.cols()-1) << '\n';
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

template <class Tbf>
void ExpansionCoeff::write2dDensity2File(const std::string & outFilename, const Tbf & bf, const int skip, const double rMax, const int nStep) const
{
	std::ofstream out(outFilename);
	for (int time = 0; time < m_coeff.size(); time += skip){ 
		if (time % (skip*10) == 0) {std::cout << "Outputting density for: " << time << '\n';}
		Eigen::ArrayXXd density = bf.densityArrayReal(m_coeff[time], nStep, rMax); // I suppose this could always be make quicker by saving the individual arrays for each bf
		outputArray(out, density);
	}
	std::cout << "Density evolution saved to: " << outFilename << '\n';
	out.close();
}

template <class Tbf>
void ExpansionCoeff::write2dDensity2File(int timeIndex, const std::string & outFilename, const Tbf & bf, const double rMax, const int nStep) const {
	std::ofstream out(outFilename);
	Eigen::ArrayXXd density = bf.densityArrayReal(m_coeff[timeIndex], nStep, rMax);;
	outputArray(out, density);
	out.close(); 

	std::cout << "Density saved to: " << outFilename << '\n'; 
}



template <class Tbf>
void ExpansionCoeff::write2dPotential2File(const std::string & outFilename, const Tbf & bf, const int skip) const
{
	std::ofstream out(outFilename);
	for (int time = 0; time < m_coeff.size(); time += skip){ 
		if (time % (skip*10) == 0) {std::cout << "Outputting potential for: " << time << '\n';}
		Eigen::ArrayXXd density = bf.potentialArrayReal(m_coeff[time], 201, 10); // I suppose this could always be make quicker by saving the individual arrays for each bf
		outputArray(out, density);
	}
	out.close();
}


void ExpansionCoeff::resizeVector(const int nTimeStep) {
	auto maxRadialIndex{m_coeff[0].size()};
	m_coeff.resize(nTimeStep); 
	for (int i = 1; i < m_coeff.size(); ++i) {m_coeff[i] = Eigen::VectorXcd::Zero(maxRadialIndex);}
}

#endif