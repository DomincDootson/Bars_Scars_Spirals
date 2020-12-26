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
	PotentialDensityPairContainer(int maxN, int l) : // Can we make it more efficent so it can figure out which basis it is??
	m_potentialDensityContainer{}, m_maxRadialIndex{maxN}, m_fourierHarmonic{l}
	{
		for (int i = 0; i <= m_maxRadialIndex ; ++i)
		{
			m_potentialDensityContainer.emplace_back(i, m_fourierHarmonic);
		}	
	}

	~PotentialDensityPairContainer() {}


	T operator()(int index) const{return m_potentialDensityContainer[index];}

	double potential(const double radius, const int n) const;
	double   density(const double radius, const int n) const;

	void scriptE(const std::string &filename) const;

	int maxRadialIndex() const {return m_maxRadialIndex;}
	int fourierHarmonic() const {return m_fourierHarmonic;}

private:	
	
	std::vector<T> m_potentialDensityContainer; // Please overload the [] operator so this can be private
	int m_maxRadialIndex,  m_fourierHarmonic;
	
	
	void saveParamters(std::ofstream & out) const;
	double scriptEelement(int k, int j) const;

};


// So down here we have our function definitions

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


// Script generating functions
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
	saveParamters(out);

	for (int i = 0; i<= m_maxRadialIndex; ++i)
	{
		for(int j = 0; j<= m_maxRadialIndex; ++j)
		{
			out << PotentialDensityPairContainer<T>::scriptEelement(i,j) << " ";
		}
		out << '\n';
	}
	out.close();
}

#endif 
