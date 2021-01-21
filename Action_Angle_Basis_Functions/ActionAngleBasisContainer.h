#ifndef ACTIONANGLEBASISCONTAINER
#define ACTIONANGLEBASISCONTAINER

// Add way to read in sciptE 
#include "ActionAngleBasisFunction.h"
#include "../DF_Class/Mestel.h" 
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <typeinfo>
#include <string>


class ActionAngleBasisContainer
{
public:
	ActionAngleBasisContainer(int maxRadialIndex, int fourierHarmonic, int maxFourierHarmonic, int sizeArray, double maxRadius)
	: m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_maxFourierHarmonic{maxFourierHarmonic},
	m_sizeArray{sizeArray}, m_maxRadius{maxRadius}, 
	m_spacingSize{maxRadius/((double) m_sizeArray - 1)}, // The size has minus 1 as we want the end point to be inclusive
	m_basisContainer{}
	{
		for (int np = 0; np <= m_maxRadialIndex; ++np)
		{
			for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++ m1)
			{
				m_basisContainer.emplace_back(np, m1, m_fourierHarmonic, m_sizeArray);
			}
		}
	}


	ActionAngleBasisContainer(std::string dir, int maxRadialIndex, int fourierHarmonic, int maxFourierHarmonic, int sizeArray, double maxRadius)
	: m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_maxFourierHarmonic{maxFourierHarmonic},
	m_sizeArray{sizeArray}, m_maxRadius{maxRadius}, 
	m_spacingSize{maxRadius/((double) m_sizeArray - 1)}, // The size has minus 1 as we want the end point to be inclusive
	m_basisContainer{}
	{
		std::cout << "Reading in basis functions.\n";
		for (int np = 0; np <= m_maxRadialIndex; ++np)
		{
			for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++ m1)
			{
				m_basisContainer.emplace_back(dir, np, m1, m_fourierHarmonic, m_sizeArray);
			}
		}
	}

	~ActionAngleBasisContainer() {}

	template <class Tbf, class Tdf>
	void scriptW(Tbf & basisFunctions, Tdf & df, std::string directory);
	double spacing() const {return m_spacingSize;}
	
	int size(int axis) const {return m_sizeArray;}
	int maxFourierHarmonic() const {return m_maxFourierHarmonic;}

	double operator()(int np, int m1, int i, int j) const {return m_basisContainer[(m_maxFourierHarmonic + m1) + (2*m_maxFourierHarmonic+1) * np](i,j);}

	template <class T>
	Eigen::MatrixXd omega1Grid(const T & distFunction) const;
	template <class T>
	Eigen::MatrixXd omega2Grid(const T & distFunction) const;


private:
	
	const int m_maxRadialIndex, m_fourierHarmonic, m_maxFourierHarmonic, m_sizeArray;
	const double m_maxRadius, m_spacingSize;
	std::vector<ActionAngleBasisFunction> m_basisContainer;

	int index(int np, int m1) const {return m1 + (2*m_maxFourierHarmonic+1) * np;}

	

	

	template <class T>
	std::vector<double>  theta1Vector(T & df, std::vector<double> & radii, double om1);
	template <class T>
	std::vector<double>  theta2Vector(T & df, std::vector<double> & radii, double om2);
	template <class T>
	std::vector<double>  theta1DerivVector(T & df, std::vector<double> & radii, double om1);

	
	
};
std::vector<double> nonLinearRadii(int steps, double rApo, double rPer)
{
	std::vector<double> radii(steps);

	double stepSize{0.5*M_PI/ ((double) steps - 1)}; // CHANGE THIS BACK

	for (int i = 0; i < steps; ++i){
		radii[i] = rPer + (rApo-rPer) * pow(sin(stepSize*i),2);
	}
	return radii;
}

template <class T>
Eigen::MatrixXd ActionAngleBasisContainer::omega1Grid(const T & distFunction) const // Calculates a grid of omega1 to read off
{
	Eigen::MatrixXd grid = Eigen::MatrixXd::Zero(size(0), size(1));
	for (int i = 0; i < size(0); ++i)
		{
			for (int j = 1; j < i; ++j)
			{
				double rApo{i*m_spacingSize}, rPer{j*m_spacingSize};
				grid(i,j) = distFunction.omega1(rApo, rPer);
			}
		}
	return grid;
}

template <class T>
Eigen::MatrixXd ActionAngleBasisContainer::omega2Grid(const T & distFunction) const // Calculates a grid of omega2 to read off
{
	Eigen::MatrixXd grid = Eigen::MatrixXd::Zero(size(0), size(1));
	for (int i = 0; i < size(0); ++i)
		{
			for (int j = 1; j < i; ++j)
			{
				double rApo{i*m_spacingSize}, rPer{j*m_spacingSize};
				grid(i,j) = distFunction.omega2(rApo, rPer);
			}
		}
	return grid;
}


template <class T>
std::vector<double> ActionAngleBasisContainer::theta1Vector(T & df, std::vector<double> & radii, double om1) // Can figure out own spacing form rApo rPer and theta1.size()
{
	std::vector<double> theta1(radii.size());
	double rPer{radii.front()}, rApo{radii.back()};
	theta1.front() = 0;
	for (int i = 1; i < ((radii.size())-1); ++i) {
		theta1[i] = df.theta1(radii[i], rApo, rPer, om1);
	}
	theta1.back() = 0;
	return theta1;
}

template <class T>
std::vector<double> ActionAngleBasisContainer::theta2Vector(T & df, std::vector<double> & radii, double om2)
{
	std::vector<double> theta2(radii.size());
	double rPer{radii.front()}, rApo{radii.back()};
	theta2.front() = 0;
	for (int i = 1; i < ((radii.size())-1); ++i) {
		theta2[i] = df.theta2(radii[i], rApo, rPer, om2);
	}
	theta2.back() = 0;
	return theta2;
}

template <class T>
std::vector<double> ActionAngleBasisContainer::theta1DerivVector(T & df, std::vector<double> & radii, double om1)
{
	std::vector<double> theta1Deriv(radii.size());
	double rPer{radii.front()}, rApo{radii.back()};
	theta1Deriv.front() = 0;
	for (int i = 1; i < ((radii.size())-1); ++i) {
		theta1Deriv[i] = df.theta1Deriv(radii[i], rApo, rPer, om1);
	}
	theta1Deriv.back() = 0;
	return theta1Deriv;
}



template <class Tbf, class Tdf>
void ActionAngleBasisContainer::scriptW(Tbf & basisFunctions, Tdf & df, std::string directory)
{
	assert(m_maxRadialIndex == basisFunctions.maxRadialIndex() && "The radial indices of the BF and Action-angle rep must be the same"); 
	assert(m_fourierHarmonic == basisFunctions.fourierHarmonic() && "The fourier harmonics of the BF and Action-angle rep must be the same");
	int nIntegrationSteps{500};
	Eigen::MatrixXd om1Grid{omega1Grid(df)}, om2Grid{omega2Grid(df)}; 
	std::vector<double> theta1(nIntegrationSteps), theta2(nIntegrationSteps), theta1Deriv(nIntegrationSteps), radii(nIntegrationSteps);

	for (int i = 1; i < om1Grid.rows(); ++i)
	{
		std::cout << "Fraction of rows completed: " << round(100*i/((double)  om1Grid.rows())) << '\n';
		for (int j = 1; j<i; ++j)
		{
			double rApo{i*m_spacingSize}, rPer{j*m_spacingSize};
			radii = nonLinearRadii(nIntegrationSteps, rApo, rPer);

			theta1 = theta1Vector(df, radii, om1Grid(i,j));
			theta2 = theta2Vector(df, radii, om2Grid(i,j));
			theta1Deriv = theta1DerivVector(df, radii, om1Grid(i,j));
		
			int count{0};
			for (int np = 0; np <= m_maxRadialIndex; ++np)
			 {
			 	for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++m1){		 		
			 		(m_basisContainer[count])(i,j) = (basisFunctions(np)).scriptWElement(m1, radii, theta1, theta2, theta1Deriv);	
				 	count += 1;
				 }
			}

		}
	} 
		for (auto bfGrid = m_basisContainer.begin(); bfGrid != m_basisContainer.end(); ++bfGrid){
		bfGrid -> save(directory, m_spacingSize); // We need to use -> as bfGrid is an interator
	}
	basisFunctions.scriptE(directory + "/scriptE_"+"Mestel_" + std::to_string(m_fourierHarmonic) + ".out");
}

#endif 