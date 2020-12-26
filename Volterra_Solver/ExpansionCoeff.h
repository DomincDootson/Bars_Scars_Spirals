#ifndef EXPANSIONCOEFF
#define EXPANSIONCOEFF 

#include <vector>
#include <cmath>

class ExpansionCoeff
{
public:
	ExpansionCoeff(int numbTimeStep) : m_coeff(numbTimeStep) {}
	
	ExpansionCoeff(const std::string & filename, int numbTimeStep) : m_coeff(numbTimeStep) 
	{coefficentReadIn(filename);}
	

	~ExpansionCoeff() {}

	
	
	Eigen::VectorXcd operator()(int timeIndex) const {return m_coeff[timeIndex];}
	
	void coefficentReadIn(const std::string &filename); 
	void write2File(const std::string & filename) const;

private:
	
	std::vector<Eigen::VectorXcd> m_coeff;
	

};
char sign(double number){
	if (number < 0){ return '-';}
	else {return '+';}
}

void ExpansionCoeff::write2File(const std::string & filename) const 
{
	std::ofstream out(filename);
	for (auto it = m_coeff.begin(); it != m_coeff.end(); ++it){
		for (int n = 0; n < (it -> size()); ++n){
			if (n == it ->size() -1){
				out << real((*it)(n)) << sign(imag((*it)(n)))  << abs(imag((*it)(n)))  << "j\n";
			}
			else{
			out << real((*it)(n)) << sign(imag((*it)(n)))  << abs(imag((*it)(n))) << "j,"; // Could we change the syntax t
			}
		}
	}
	out.close();
}

#endif