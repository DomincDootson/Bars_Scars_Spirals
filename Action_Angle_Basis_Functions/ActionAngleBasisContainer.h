#ifndef ACTIONANGLEBASISCONTAINER
#define ACTIONANGLEBASISCONTAINER

#include "ActionAngleBasisFunction.h"
#include "../DF_Class/Mestel.h" 
#include <vector>
#include <Eigen/Dense>
#include <string>

// Could we make this more efficent to have the grids that keep getting passed around, saved as memeber variables? 
class ActionAngleBasisContainer
{
public:
	ActionAngleBasisContainer(std::string prefix, int maxRadialIndex, int fourierHarmonic, int maxFourierHarmonic, int sizeArray, double maxRadius)
	: m_prefix{prefix},
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_maxFourierHarmonic{maxFourierHarmonic},
	m_sizeArray{sizeArray},
	m_spacingSize{maxRadius/((double) m_sizeArray - 1)}, // The size has minus 1 as we want the end point to be inclusive
	m_basisContainer{},
	v_radii(500), v_theta1(500), v_theta2(500), v_theta1Deriv(500), 
	m_scriptE(m_maxRadialIndex+1, m_maxRadialIndex+1)
	{		
		for (int np = 0; np <= m_maxRadialIndex; ++np)
		{
			for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++ m1)
			{
				m_basisContainer.emplace_back(np, m1, m_fourierHarmonic, m_sizeArray);
			}
		}
	}


	ActionAngleBasisContainer(std::string dir, std::string prefix, int maxRadialIndex, int fourierHarmonic, int maxFourierHarmonic, int sizeArray, double maxRadius)
	: m_prefix{prefix},
	m_maxRadialIndex{maxRadialIndex}, m_fourierHarmonic{fourierHarmonic}, m_maxFourierHarmonic{maxFourierHarmonic},
	m_sizeArray{sizeArray}, 
	m_spacingSize{maxRadius/((double) m_sizeArray - 1)}, // The size has minus 1 as we want the end point to be inclusive
	m_basisContainer{},
	m_scriptE{ readInScriptE(dir)}
	{
		std::cout << "Reading in basis functions.\n";
		for (int np = 0; np <= m_maxRadialIndex; ++np)
		{
			for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++ m1)
			{
				m_basisContainer.emplace_back(dir, prefix, np, m1, m_fourierHarmonic, m_sizeArray);
			}
		}
	}

	~ActionAngleBasisContainer() {}

	template <class Tbf, class Tdf>
	void scriptW(Tbf & basisFunctions, Tdf & df, std::string directory);
	double spacing() const {return m_spacingSize;}
	
	int size(int axis) const {return m_sizeArray;}
	int fourierHarmonic() const {return m_fourierHarmonic;}
	int maxFourierHarmonic() const {return m_maxFourierHarmonic;}
	int maxRadialIndex() const {return m_maxRadialIndex;}
	double step() const {return m_spacingSize;}

	double operator()(int np, int m1, int i, int j) const {return m_basisContainer[(m_maxFourierHarmonic + m1) + (2*m_maxFourierHarmonic+1) * np](i,j);}
	Eigen::MatrixXd operator()(int np, int m1) const {return m_basisContainer[(m_maxFourierHarmonic + m1) + (2*m_maxFourierHarmonic+1) * np].get();}

	template <class T>
	Eigen::MatrixXd omega1Grid(const T & distFunction) const;
	template <class T>
	Eigen::MatrixXd omega2Grid(const T & distFunction) const;

	Eigen::MatrixXd inverseScriptE() const {return m_scriptE.inverse();};


private:
	
	const std::string m_prefix; 
	const int m_maxRadialIndex, m_fourierHarmonic, m_maxFourierHarmonic, m_sizeArray;
	const double m_spacingSize;
	std::vector<ActionAngleBasisFunction> m_basisContainer;
	std::vector<double> v_radii, v_theta1, v_theta2, v_theta1Deriv; 
	Eigen::MatrixXd m_scriptE, m_om1Grid, m_om2Grid;

	int index(int np, int m1) const {return m1 + (2*m_maxFourierHarmonic+1) * np;}
	
	void nonLinearRadii(double rApo, double rPer); 
	template <class T>
	void theta1Vector(T & df, std::vector<double> & radii, double om1);
	template <class T>
	void theta2Vector(T & df, std::vector<double> & radii, double om2);
	template <class T>
	void theta1DerivVector(T & df, std::vector<double> & radii, double om1);

	
	void checkParams(int maxRadialIndex, int fourierHarmonic);
	template <class Tdf>
	void omegaGridSetUp(const Tdf & df);
	template <class Tdf>
	void vectorsAtRadii(const Tdf & df, int rApoIndex, int rPerIndex);
	template <class Tbf>
	void scriptWLoop(const Tbf & basisFunctions, int i, int j);
	template <class Tbf>
	void write2file(const std::string & dir, const Tbf & basisFunctions); 


	std::string scriptEfilename(const std::string dir) const 
	{return dir + "/scriptE_" + std::to_string(m_maxRadialIndex) + '_' + std::to_string(m_fourierHarmonic) + ".out";}
	Eigen::MatrixXd readInScriptE(const std::string dir);
	
};
void ActionAngleBasisContainer::nonLinearRadii(double rApo, double rPer)
{
	double stepSize{0.5*M_PI/ ((double) v_radii.size() - 1)};
	for (int i = 0; i < v_radii.size(); ++i){
		v_radii[i] = rPer + (rApo-rPer) * pow(sin(stepSize*i),2);
	}
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
void ActionAngleBasisContainer::theta1Vector(T & df, std::vector<double> & radii, double om1) // Can figure out own spacing form rApo rPer and theta1.size()
{
	double rPer{radii.front()}, rApo{radii.back()};
	v_theta1.front() = 0;
	for (int i = 1; i < ((radii.size())-1); ++i) {
		v_theta1[i] = df.theta1(radii[i], rApo, rPer, om1);
	}
	v_theta1.back() = 0;
}

template <class T>
void ActionAngleBasisContainer::theta2Vector(T & df, std::vector<double> & radii, double om2)
{
	double rPer{radii.front()}, rApo{radii.back()};
	v_theta2.front() = 0;
	for (int i = 1; i < ((radii.size())-1); ++i) {
		v_theta2[i] = df.theta2(radii[i], rApo, rPer, om2);
	}
	v_theta2.back() = 0;
}

template <class T>
void ActionAngleBasisContainer::theta1DerivVector(T & df, std::vector<double> & radii, double om1)
{
	double rPer{radii.front()}, rApo{radii.back()};
	v_theta1Deriv.front() = 0;
	for (int i = 1; i < ((radii.size())-1); ++i) {
		v_theta1Deriv[i] = df.theta1Deriv(radii[i], rApo, rPer, om1);
	}
	v_theta1Deriv.back() = 0;
}

void ActionAngleBasisContainer::checkParams(int maxRadialIndex, int fourierHarmonic) {
	assert(m_maxRadialIndex == maxRadialIndex && "The radial indices of the BF and Action-angle rep must be the same"); 
	assert(m_fourierHarmonic == fourierHarmonic && "The fourier harmonics of the BF and Action-angle rep must be the same");
}

template <class Tdf>
void ActionAngleBasisContainer::omegaGridSetUp(const Tdf & df) {
	m_om1Grid = omega1Grid(df); 
	m_om2Grid = omega2Grid(df);	
}

template <class Tdf>
void ActionAngleBasisContainer::vectorsAtRadii(const Tdf & df, int rApoIndex, int rPerIndex) {
	double rApo{rApoIndex * m_spacingSize}, rPer{rPerIndex * m_spacingSize};
	nonLinearRadii(rApo, rPer);
	theta1Vector(df, v_radii, m_om1Grid(rApoIndex, rPerIndex));
	theta2Vector(df, v_radii, m_om2Grid(rApoIndex, rPerIndex));
	theta1DerivVector(df, v_radii, m_om1Grid(rApoIndex, rPerIndex));
}

template <class Tbf>
void ActionAngleBasisContainer::scriptWLoop(const Tbf & basisFunctions, int i, int j){
	int count{0}; 
	for (int np = 0; np <= m_maxRadialIndex; ++np)
	 {
		for (int m1 = -m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++m1){		 		
			(m_basisContainer[count])(i,j) = (basisFunctions(np)).scriptWElement(m1, v_radii, v_theta1, v_theta2, v_theta1Deriv);	
		 	count += 1;
		}
	}
}

template <class Tbf>
void ActionAngleBasisContainer::write2file(const std::string & dir, const Tbf & basisFunctions) {
	for (auto bfGrid = m_basisContainer.begin(); bfGrid != m_basisContainer.end(); ++bfGrid){ // This to end becomes a function
		bfGrid -> save(dir, m_prefix, m_spacingSize); // We need to use -> as bfGrid is an interator
	}
	basisFunctions.scriptE(scriptEfilename(dir));
}


template <class Tbf, class Tdf> 
void ActionAngleBasisContainer::scriptW(Tbf & basisFunctions, Tdf & df, std::string directory) {
	checkParams(basisFunctions.maxRadialIndex(), basisFunctions.fourierHarmonic());
	omegaGridSetUp(df);

	for (int i = 1; i < m_om1Grid.rows(); ++i)
	{
		std::cout << "Fraction of rows completed: " << .1*round(1000*i/((double)  m_om1Grid.rows())) << "%" << '\n';
		for (int j = 1; j<i; ++j)
		{
			vectorsAtRadii(df, i, j);
			scriptWLoop(basisFunctions, i, j); 
		}
	} 
	write2file(directory, basisFunctions); 	
}

Eigen::MatrixXd ActionAngleBasisContainer::readInScriptE(const std::string dir)
{
	std::ifstream inScriptE(scriptEfilename(dir));
	Eigen::MatrixXd scriptE(m_maxRadialIndex+1, m_maxRadialIndex+1);

	for (int i =0; i < scriptE.rows(); ++i)
	{
		for (int j = 0; j <scriptE.cols(); ++j){
			 inScriptE >> scriptE(i,j);
		}
	}
	inScriptE.close();
	return scriptE;

}




#endif 