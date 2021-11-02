#ifndef POTENTIALDENSITYPAIRCONTAINER
#define POTENTIALDENSITYPAIRCONTAINER

#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <typeinfo> // This allows us to find the type of a object
// We want to inlude some file management
#include <iostream>


template <class T>
class PotentialDensityPairContainer
{
public:
	PotentialDensityPairContainer(std::vector<double> params, int maxN, int l) :
	m_potentialDensityContainer{}, m_maxRadialIndex{maxN}, m_fourierHarmonic{l}, 
	m_scriptE{Eigen::MatrixXd::Zero(m_maxRadialIndex+1, m_maxRadialIndex+1)}
	{
		for (int i = 0; i <= m_maxRadialIndex; ++i){
			m_potentialDensityContainer.emplace_back(params, i, m_fourierHarmonic);
		}
		m_scriptE = calculateScriptE();
	}

	~PotentialDensityPairContainer() {}



	T operator()(int index) const{return m_potentialDensityContainer[index];}

	double potential(const double radius, const int n) const;
	double   density(const double radius, const int n) const;

	void scriptE(const std::string &filename) const;
	Eigen::MatrixXd calculateScriptE() const;

	Eigen::MatrixXd scriptE() const {return m_scriptE;}

	int maxRadialIndex() const {return m_maxRadialIndex;}
	int fourierHarmonic() const {return m_fourierHarmonic;}

	Eigen::ArrayXXcd   densityArray(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const; 
	Eigen::ArrayXXcd potentialArray(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const;

	std::vector<Eigen::ArrayXXcd> individualPotential(const int nGrid, const double rMax) const;

	Eigen::ArrayXXd   densityArrayReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const;
	Eigen::ArrayXXd potentialArrayReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const;



	std::vector<double> oneDdensity(const std::vector<double> & radii, const Eigen::VectorXcd &coeff) const;
	std::vector<double> oneDpotential(const std::vector<double> & radii, const Eigen::VectorXcd &coeff) const;


	Eigen::VectorXcd potentialResolving(const Eigen::ArrayXXcd &potentialArray, const double rMax) const; // Not great at resolving l =0
	Eigen::VectorXcd densityResolving(const Eigen::ArrayXXcd &densityArray, const double rMax) const;

	Eigen::VectorXcd potentialFitting(const Eigen::ArrayXXcd &potentialArrayTrue, const double rMax) const;  // This will calculate individual parts and solve matrix equation



private:	
	
	std::vector<T> m_potentialDensityContainer;
	int m_maxRadialIndex,  m_fourierHarmonic;
	Eigen::MatrixXd m_scriptE;
	
	// Function that calculates the matrix elements for fitting
	std::complex<double> potentialFittingElement(const double rMax, const int nGrid, const int i, const int j) const;
	Eigen::MatrixXcd     potentialFittingMatrix(const double rMax, const int nGrid) const; 
	std::complex<double>     potentialInhomogenity(const Eigen::ArrayXXcd &potentialArrayTrue, const double rMax, const int i) const;

	void saveParamters(std::ofstream & out) const;
	double scriptEelement(int k, int j) const;

};


// Our potential and density function

template <class T>
double PotentialDensityPairContainer<T>::potential(const double radius, const int n) const
{
	assert(n<=m_maxRadialIndex);
	return m_potentialDensityContainer[n].potential(radius);
}

template <class T>
double   PotentialDensityPairContainer<T>::density(const double radius, const int n) const
{
	assert(n<=m_maxRadialIndex);
	return m_potentialDensityContainer[n].density(radius);
}

// Grid generation Function

template <class T>
Eigen::ArrayXXcd   PotentialDensityPairContainer<T>::densityArray(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const
{
	Eigen::ArrayXXcd grid = Eigen::ArrayXXcd::Zero(nGrid, nGrid);
	for (int i = 0; i <= m_maxRadialIndex; ++i){
		grid += coefficents(i) * m_potentialDensityContainer[i].densityGrid(nGrid, rMax);
	}
	return grid;
}


template <class T>
Eigen::ArrayXXcd PotentialDensityPairContainer<T>::potentialArray(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const
{
	Eigen::ArrayXXcd grid = Eigen::ArrayXXcd::Zero(nGrid, nGrid);
	for (int i = 0; i <= m_maxRadialIndex; ++i){
		grid += coefficents(i) * m_potentialDensityContainer[i].potentialGrid(nGrid, rMax);
	}
	return grid;
}

template <class T>
std::vector<Eigen::ArrayXXcd> PotentialDensityPairContainer<T>::individualPotential(const int nGrid, const double rMax) const {
	std::vector<Eigen::ArrayXXcd> individualPot;
	for (int n = 0; n < m_maxRadialIndex; ++n){
		individualPot.emplace_back(nGrid, nGrid);
		individualPot[n] = m_potentialDensityContainer[n].potentialGrid(nGrid, rMax);
	}
	return individualPot;
}


template <class T>
Eigen::ArrayXXd   PotentialDensityPairContainer<T>::densityArrayReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const{
	return 2*(densityArray(coefficents, nGrid, rMax)).real();
}

template <class T>
Eigen::ArrayXXd PotentialDensityPairContainer<T>::potentialArrayReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const{
	return 2*(potentialArray(coefficents, nGrid, rMax)).real(); 
}

template <class T>
std::vector<double> PotentialDensityPairContainer<T>::oneDdensity(const std::vector<double> & radii, const Eigen::VectorXcd &coeff) const
{
	std::vector<double> densityVector;
	for (auto r = radii.begin(); r != radii.end(); ++r){
			std::complex<double> den{0};
			for (int n = 0; n <= m_maxRadialIndex; ++n){
				den += coeff[n] * density(*r, n);
			}
			densityVector.push_back(real(den));
		}
	return densityVector;	
}

template <class T>
std::vector<double> PotentialDensityPairContainer<T>::oneDpotential(const std::vector<double> & radii, const Eigen::VectorXcd &coeff) const
{
	std::vector<double> potentialVector;
	for (auto r = radii.begin(); r != radii.end(); ++r){
			std::complex<double> pot{0};
			for (int n = 0; n <= m_maxRadialIndex; ++n){
				pot += coeff[n] * potential(*r, n);
			}
			potentialVector.push_back(real(pot));
		}
	return potentialVector;		
}

// Resolving function
template <class T>
Eigen::VectorXcd PotentialDensityPairContainer<T>::potentialResolving(const Eigen::ArrayXXcd &potentialArray, const double rMax) const
{
	double spacing{2*rMax/((double) potentialArray.rows()-1)};
	Eigen::VectorXcd coefficents(m_maxRadialIndex+1);
	for (int i = 0; i <= m_maxRadialIndex; ++i)
	{
		Eigen::ArrayXXcd densityArray{((m_potentialDensityContainer[i].densityGrid(potentialArray.rows(), rMax)).conjugate())};
		coefficents(i) = spacing*spacing*(densityArray * potentialArray).sum();
	}
	return -(m_scriptE.inverse())*coefficents; 
}

template <class T>
Eigen::VectorXcd PotentialDensityPairContainer<T>::densityResolving(const Eigen::ArrayXXcd &densityArray, const double rMax) const
{
	double spacing{2*rMax/((double) densityArray.rows()-1)};
	Eigen::VectorXcd coefficents(m_maxRadialIndex+1);
	for (int i = 0; i <= m_maxRadialIndex; ++i)
	{
		Eigen::ArrayXXcd potentialArray{((m_potentialDensityContainer[i].potentialGrid(densityArray.rows(), rMax)).conjugate())};
		coefficents(i) = spacing*spacing*(potentialArray * densityArray).sum();
	}
	return -(m_scriptE.inverse())*coefficents; 
}


// ScriptE Generating functions

template <class T>
void PotentialDensityPairContainer<T>::saveParamters(std::ofstream & out) const
{
	out<< typeid(m_potentialDensityContainer[0]).name() << " " << m_maxRadialIndex << " " << m_fourierHarmonic << '\n'; // we also want to include a basis function tag
}

template <class T>
double PotentialDensityPairContainer<T>::scriptEelement(int k, int j) const // k - potential index, j - density index
{
	int steps{10000};
	double stepSize{25/((double) steps)}, scriptEelement{0}, r{};
	for (int i = 1; i < steps; ++i){
		r = i * stepSize;
		scriptEelement += (2*M_PI*r)*stepSize*potential(r, k)*density(r, j);
	}
	return -scriptEelement;
}

template <class T>
void PotentialDensityPairContainer<T>::scriptE(const std::string &filename) const
{
	std::ofstream out(filename);
	

	for (int i = 0; i<= m_maxRadialIndex; ++i)
	{
		for(int j = 0; j<= m_maxRadialIndex; ++j)
		{
			out << (m_scriptE(i,j)) << " ";
		}
		out << '\n';
	}
	out.close();
}

template <class T>
Eigen::MatrixXd PotentialDensityPairContainer<T>::calculateScriptE() const
{
	Eigen::MatrixXd scriptE = Eigen::MatrixXd::Zero(m_maxRadialIndex+1, m_maxRadialIndex+1);
	for (int i = 0; i<= m_maxRadialIndex; ++i)
	{
		for(int j = 0; j<= m_maxRadialIndex; ++j)
		{
			scriptE(i,j) = PotentialDensityPairContainer<T>::scriptEelement(i,j);
		}
	}
	return scriptE; 
}

template <class T> 
Eigen::VectorXcd PotentialDensityPairContainer<T>::potentialFitting(const Eigen::ArrayXXcd &potentialArrayTrue, const double rMax) const {
	Eigen::VectorXcd inhomo{m_maxRadialIndex+1}; Eigen::MatrixXcd potentialMatrix = potentialFittingMatrix(rMax, potentialArrayTrue.cols());
	for (int i = 0; i < m_maxRadialIndex+1; ++i) {inhomo(i) = potentialInhomogenity(potentialArrayTrue, rMax, i);} 

	
	return (potentialMatrix.inverse()) * inhomo;
}

template <class T> 
std::complex<double> PotentialDensityPairContainer<T>::potentialFittingElement(const double rMax, const int nGrid, const int i, const int j) const {
	double spacing{(2*rMax)/((double) nGrid)};
	Eigen::ArrayXXcd arrayI{m_potentialDensityContainer[i].potentialGrid(nGrid, rMax).conjugate()}, arrayJ{m_potentialDensityContainer[j].potentialGrid(nGrid, rMax)};
	return spacing*spacing*((arrayI*arrayJ).sum());
}

template <class T> 
Eigen::MatrixXcd PotentialDensityPairContainer<T>::potentialFittingMatrix(const double rMax, const int nGrid) const {
	Eigen::MatrixXcd matrix = Eigen::MatrixXcd::Zero(m_maxRadialIndex+1, m_maxRadialIndex+1);
	for (int i = 0; i < m_maxRadialIndex+1; ++i) {for (int j = i; j < m_maxRadialIndex+1; ++j) {
		matrix(i, j) = potentialFittingElement(rMax, nGrid, i, j); 
		matrix(j, i) = matrix(i, j);
	}}
	std::cout << "det is: " << matrix.determinant() << '\n'; 
	return matrix; 
}

void saveArray(const Eigen::ArrayXXcd grid) {
	std::ofstream out("Plotting/inhomo.csv");
	for (int i = 0; i < grid.rows(); ++i) {
		for (int j = 0; j < grid.rows()-1; ++j) { 
			out << grid(i,j).real() << ',';
		}
		out << grid(i, grid.cols()-1).real() << '\n';
	}
	out.close();
}


template <class T> 
std::complex<double> PotentialDensityPairContainer<T>::potentialInhomogenity(const Eigen::ArrayXXcd &potentialArrayTrue, const double rMax, const int i) const {
	int nGrid{(int) potentialArrayTrue.cols()}; double spacing{(2*rMax)/((double) nGrid)};
	Eigen::ArrayXXcd arrayI{(m_potentialDensityContainer[i].potentialGrid(nGrid, rMax)).conjugate()}; // 

	//saveArray(arrayI*potentialArrayTrue); 

	return spacing*spacing*(arrayI*potentialArrayTrue).sum();
}

#endif 
