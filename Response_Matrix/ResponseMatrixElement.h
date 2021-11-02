#ifndef RESPONSEMATRIXELEMENT
#define RESPONSEMATRIXELEMENT

#include <complex>

class ResponseMatrixElement
{
public:
	ResponseMatrixElement() : 
	m_omega1Grid{df.omega1Grid()}, m_omega2Grid{df.omega2Grid()}, 
	m_dfdEGrid{df.dFdEgrid(, m_omega1Grid)}, m_dfdJGrid{df.dFdJgrid( , m_omega2Grid)},
	m_elJacobian{energyAngMomJacobain} 
	{}
	~ResponseMatrixElement() {} 
	
	std::complex<double> responseElement(int i, int j, int m1); // Effectively this calulates aleph  


private:
	
	const double m_spacing; 
	const Eigen::MatrixXd m_omega1Grid, m_omega2Grid, m_dfdEGrid, m_dfdJGrid, m_elJacobian; 
	Eigen::ArrayXXcd m_gGrid, m_hGrid; 
	const ActionAngleBasisContainer m_basisContainer; 
	std::vector<double> v_coefficents; 

	void rescale();

	std::complex<double> GFunction(double x, double y) const;
	std::<double> responseElement(int i, int j, int m1); 

	/*
	I think a vector to contain all the elements that are needed 
	some functions that calculate g and h
	some function that gets the values of a,b,c, |A function that does both of these
	some function that does the rescling 		 | 	

	some function that will caluclate G       
	*/ 
};


void ResponseMatrixElement::rescale() {
	double ag{v_coefficents[0]}, ah{v_coefficents[3]}, deltaR{v_coefficents.back()};
	v_coefficents[1] *= deltaR/ag; v_coefficents[2] *= deltaR/ag; 
	v_coefficents[4] *= deltaR/ah; v_coefficents[5] *= deltaR/ah; v_coefficents[6] *= (1/ah);
}


std::complex<double> ResponseMatrixElement::GFunction(double x, double y) const { // PLEASE, FOR THE LOVE OF GOD, CHECK THIS FUNCTION 
	double b{v_coefficents[1]}, c{v_coefficents[2]}, e{v_coefficents[4]}, f{v_coefficents[5]}, eta{v_coefficents[6]};
	std::complex<double> i{0,1}; 
	return (1/(4*e*e*f*f)) * log(e*e*x*x +2*e*(f*x*y + x) +f*f*y*y + 2*f*y + eta*eta +1) *
	(b*f*(e*e*x*x - pow(f*y+i*eta+1,2)) + 2*e*f*(e*x+i*eta + 1) - c*e*pow(e*x+i*eta+1, 2)) + 
	(i/(2*e*e*f*f)) * (0.5*M_PI-arctan((e*x+f*y+1) / eta)) *
	(b*f*(e*e*x*x - pow(f*y+i*eta+1,2)) + 2*e*f*(e*x+i*eta + 1) - c*e*pow(e*x+i*eta+1, 2)) + 
	(y/(4*e*e*f)) * (f*(-4*e + b*(2*e*x + f*y +2*i*eta +2)) + 
		c*e*(2*e*x-f*y +2*i*eta+2) + 2*e*f*(c*y+2)*log(e*x+f*y+i*eta+1)); 
}

std::complex<double> ResponseMatrixElement::responseElement(int i, int j, int m1) {
	// call some function that generates all the values of the coefficents

	return (v_coefficents[0]/v_coefficents[3]) * v_coefficents.back()* v_coefficents.back() *
	(GFunction(0.5, 0.5) + GFunction(-0.5, -0.5) - GFunction(-0.5, 0.5) - GFunction(0.5, -0.5));
}

#endif