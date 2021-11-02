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

	ActionAngleBasisFunction(const std::string & dir, const std::string & prefix, int np, int m1, int m2, int sizeArray);

	

	~ActionAngleBasisFunction() {}
	
	
	void save(const std::string & directory, const std::string & prefix, double step) const;
	double operator()(int i, int j) const {return m_basisGrid(i,j);}
	double& operator()(int i, int j) {return m_basisGrid(i,j);} // This allows us to assign values.
	Eigen::MatrixXd get() const {return m_basisGrid;}

	int m1() const {return m_m1;}
private:
	Eigen::MatrixXd m_basisGrid;
	const int m_radialIndex, m_m1, m_m2; 
	std::string filename(std::string directory, std::string prefix) const;
	void checkParam(int basisRows, double step, int radialIndex, int m1, int m2) const;
	
};

ActionAngleBasisFunction::ActionAngleBasisFunction(const std::string & dir, const std::string & prefix, int np, int m1, int m2, int sizeArray)
	: m_basisGrid{Eigen::MatrixXd::Zero(sizeArray, sizeArray)}, m_radialIndex{np}, m_m1{m1}, m_m2{m2} 
{
	std::ifstream inFile;
	inFile.open(filename(dir, prefix));
	int rbasisRows{}, rradialIndex{}, rm1{}, rm2{}; 
	double rstep{};
	inFile >> rbasisRows >> rstep >> rradialIndex >> rm1 >> rm2;
	std::cout << "Reading in action-angle basis function from: " << filename(dir, prefix) << '\n';
	checkParam(rbasisRows, rstep, rradialIndex, rm1, rm2);
	// some test function
	
	for (int i = 1; i<m_basisGrid.rows(); ++i){
		for (int j = 1; j < i; ++j){
			double holding; inFile >> holding;
			if (holding != holding) {std::cout << "We have read in Nan\n";}
			m_basisGrid(i,j) = holding;	
		}
	}
	inFile.close();
}

void ActionAngleBasisFunction::checkParam(int basisRows, double step, int radialIndex, int m1, int m2) const
{
	assert(basisRows == m_basisGrid.rows() && "Wrong number of rows in basis function.\n");
	assert(radialIndex == m_radialIndex && "Wrong radialIndex read in.\n");
	assert(m1 == m_m1 && "Read in wrong m1.\n");
	assert(m2 == m_m2 && "Read in wrong m2.\n");
}



std::string ActionAngleBasisFunction::filename(std::string directory, std::string prefix) const
{
	return directory + "/" + prefix + "_" +std::to_string(m_radialIndex) + "_" + std::to_string(m_m1) + "_" + std::to_string(m_m2) + ".out";
	//GaussianLog
	//Kalnajs
}



void ActionAngleBasisFunction::save(const std::string & directory, const std::string & prefix, double step) const
{
	std::string file = filename(directory, prefix);
	std::cout << file << '\n';
	std::ofstream out(file);

	out << m_basisGrid.rows() << " " << step << " " << m_radialIndex << " " << m_m1 << " " << m_m2 << '\n';
	
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