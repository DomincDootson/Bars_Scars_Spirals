#ifndef ACTIONANGLEBASISFUNCTION
#define ACTIONANGLEBASISFUNCTION 

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>

#include <complex>


class ActionAngleBasisFunction
{
public:
	ActionAngleBasisFunction(int np, int m1, int m2, int sizeArray) 
	: m_basisGrid{Eigen::MatrixXd::Zero(sizeArray, sizeArray)}, m_radialIndex{np}, m_m1{m1}, m_m2{m2} {}

	~ActionAngleBasisFunction() {}
	
	
	void save(const std::string directory, const std::string DF, double step) const;
	double operator()(int i, int j) const {return m_basisGrid(i,j);}
	double& operator()(int i, int j) {return m_basisGrid(i,j);} // This allows us to assign values.
private:
	Eigen::MatrixXd m_basisGrid;// I'm not sure I want this public... I think the trick is to overload the () operator. 
	const int m_radialIndex, m_m1, m_m2; 
	std::string filename(std::string directory) const;
	
};

std::string ActionAngleBasisFunction::filename(std::string directory) const
{
	return directory + "/" + directory + "_" +std::to_string(m_radialIndex) + "_" + std::to_string(m_m1) + "_" + std::to_string(m_m2) + ".out";
}



void ActionAngleBasisFunction::save(const std::string directory, const std::string DF, double step) const
{
	std::string file = filename(directory);
	std::cout << file << '\n';
	std::ofstream out(file);

	out << DF << " " << m_basisGrid.rows() << " " << step << " " << m_radialIndex << " " << m_m1 << " " << m_m2 << '\n';
	
	for (int i = 1; i < m_basisGrid.rows(); ++i) // Row, First non-zero value on row 2
		{

			for (int j = 1; j < i; ++j) // Col
			{
				if ((j == i-1)&&(i < m_basisGrid.rows()-1))
				{
					out << m_basisGrid(i,j) << '\n';
				}
				else if ((j == i-1)&&(i == m_basisGrid.rows()-1))
				{
					out << m_basisGrid(i,j); 
				} 
				else
				{
					out << m_basisGrid(i,j) << " ";
				}
			}
		}
	out.close();
}

#endif