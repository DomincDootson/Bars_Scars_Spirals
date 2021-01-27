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
	m_scriptE{Eigen::MatrixXcd::Zero(m_maxRadialIndex+1, m_maxRadialIndex+1)}
	{
		for (int i = 0; i <= m_maxRadialIndex ; ++i)
		{
			m_potentialDensityContainer.emplace_back(params, i, m_fourierHarmonic);
		}
		m_scriptE = scriptE();
	}

	~PotentialDensityPairContainer() {}



	T operator()(int index) const{return m_potentialDensityContainer[index];}

	double potential(const double radius, const int n) const;
	double   density(const double radius, const int n) const;

	void scriptE(const std::string &filename) const;
	Eigen::MatrixXcd scriptE() const;

	int maxRadialIndex() const {return m_maxRadialIndex;}
	int fourierHarmonic() const {return m_fourierHarmonic;}

	Eigen::ArrayXXcd   densityGrid(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const; 
	Eigen::ArrayXXcd potentialGrid(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const;

	Eigen::ArrayXXd   densityGridReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const;
	Eigen::ArrayXXd potentialGridReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const;

	Eigen::VectorXcd potentialResolving(const Eigen::ArrayXXcd &potentialArray, const double rMax) const; // Not great at resolving l =0
	Eigen::VectorXcd densityResolving(const Eigen::ArrayXXcd &densityArray, const double rMax) const;

	Eigen::MatrixXcd getScriptE() {return m_scriptE;}

private:	
	
	std::vector<T> m_potentialDensityContainer;
	int m_maxRadialIndex,  m_fourierHarmonic;
	Eigen::MatrixXcd m_scriptE;
	
	
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
Eigen::ArrayXXcd   PotentialDensityPairContainer<T>::densityGrid(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const
{
	Eigen::ArrayXXcd grid = Eigen::ArrayXXcd::Zero(nGrid, nGrid);
	for (int i = 0; i <= m_maxRadialIndex; ++i){
		grid += coefficents(i) * m_potentialDensityContainer[i].densityGrid(nGrid, rMax);
	}
	return grid;
}


template <class T>
Eigen::ArrayXXcd PotentialDensityPairContainer<T>::potentialGrid(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const
{
	Eigen::ArrayXXcd grid = Eigen::ArrayXXcd::Zero(nGrid, nGrid);
	for (int i = 0; i <= m_maxRadialIndex; ++i){
		grid += coefficents(i) * m_potentialDensityContainer[i].potentialGrid(nGrid, rMax);
	}
	return grid;
}

template <class T>
Eigen::ArrayXXd   PotentialDensityPairContainer<T>::densityGridReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const{
	return (potentialGrid(coefficents, nGrid, rMax) + potentialGrid(coefficents, nGrid, rMax).conj()).real();
}


template <class T>
Eigen::ArrayXXd PotentialDensityPairContainer<T>::potentialGridReal(const Eigen::VectorXcd &coefficents, const int nGrid, const double rMax) const{
	return (densityGrid(coefficents, nGrid, rMax) + densityGrid(coefficents, nGrid, rMax).conj()).real(); 
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
	return -(m_scriptE.inverse())*coefficents; // WE NEED TO DO MULTIPLICATION BY SCRIPTE
	// Do mulitplication by E 
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
	return -(m_scriptE.inverse())*coefficents; // WE NEED TO DO MULTIPLICATION BY SCRIPTE
	// Do mulitplication by E 
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
			out << real(m_scriptE(i,j)) << " ";
		}
		out << '\n';
	}
	out.close();
}

template <class T>
Eigen::MatrixXcd PotentialDensityPairContainer<T>::scriptE() const
{
	Eigen::MatrixXcd scriptE = Eigen::MatrixXcd::Zero(m_maxRadialIndex+1, m_maxRadialIndex+1);
	for (int i = 0; i<= m_maxRadialIndex; ++i)
	{
		for(int j = 0; j<= m_maxRadialIndex; ++j)
		{
			scriptE(i,j) = PotentialDensityPairContainer<T>::scriptEelement(i,j);
		}
	}
	return scriptE; 
}


#endif 
