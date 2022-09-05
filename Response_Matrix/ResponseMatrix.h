#ifndef RESPONSEMATRIX
#define RESPONSEMATRIX

#include <complex>
#include <vector> 

#include <Eigen/Eigenvalues> 

#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../Volterra_Solver/EvolutionKernels.h"

class ResponseMatrix
{
public:
	template <class Tdf>
	ResponseMatrix(const ActionAngleBasisContainer & BF, const Tdf & df) : m_BF{BF},
	m_intStep{m_BF.size(0)}, m_fourierHarmonic{m_BF.fourierHarmonic()}, m_maxFourierHarmonic{m_BF.maxFourierHarmonic()}, 
	m_spacing{m_BF.step()}, 
	v_harmonicDotOmega(2 * m_maxFourierHarmonic + 1), v_gPrefactor(2 * m_maxFourierHarmonic + 1), v_taylorCoeff()
	{harmonicDotOmega(df);
	
	v_taylorCoeff.reserve(7);
	for (int i = 0; i < 6+1; ++i) {v_taylorCoeff.emplace_back(m_intStep, m_intStep); v_taylorCoeff[i] *= 0;}
	}

	ResponseMatrix(const ActionAngleBasisContainer & BF) : m_BF{BF}, m_intStep{m_BF.size(0)}, m_fourierHarmonic{m_BF.fourierHarmonic()}, 
	m_maxFourierHarmonic{m_BF.maxFourierHarmonic()}, m_spacing{m_BF.step()} 
	{}

	ResponseMatrix() : m_BF{}, m_intStep{}, m_fourierHarmonic{}, m_maxFourierHarmonic{}, m_spacing{} {}
	
	~ResponseMatrix() {}

	Eigen::MatrixXcd responseMatrix(const std::complex<double> & omega);
	
	template <class Tdf>
	Eigen::MatrixXcd responseMatrix(const std::complex<double> & omegas, EvolutionKernels & kernel, const Tdf & df) {m_responseMatrix = kernel.responseMatrix(omegas, df, m_BF); return m_responseMatrix;}

	Eigen::MatrixXcd responseMatrix(const std::complex<double> & omegas, const EvolutionKernels & kernel) {m_responseMatrix = kernel.responseMatrix(omegas); return m_responseMatrix;}

	void responseMatrixElement(const std::complex<double> & omega, int p, int q);
	
	std::complex<double> det() const;
	std::complex<double> det(const std::complex<double> & omega) {responseMatrix(omega); return det();} 

	void saveResponseMatrix(const std::string & filename) const;
private:
	// Some internal parameters
	const ActionAngleBasisContainer m_BF; 
	const int m_intStep, m_fourierHarmonic, m_maxFourierHarmonic; 
	const double m_spacing;

	std::vector<Eigen::ArrayXXd> v_harmonicDotOmega, v_gPrefactor; // equ 26. & second half of 27. 
	std::vector<Eigen::ArrayXXd> v_taylorCoeff; 
	Eigen::MatrixXcd m_responseMatrix; 

	template <class Tdf>
	void harmonicDotOmega(const Tdf & df); 
	
	
	void calculateNormalisedCoeff(const std::complex<double> & omega, int m1, int p, int q);
	void getDerivCoeff();
	void normaliseCoeff();
	void normaliseTwoArrays(int n, int d); // This devides one array, ingnoring the points we don't need (and no devide by zero!)
	std::complex<double> alephSummedOver() const ;
	std::complex<double> Gfunctions(int i, int j) const; // Checked against JB 

	void dumpCoeff(int i, int j) const {std::cout << "Coefficents: "; for (auto each : v_harmonicDotOmega) {std::cout << each(i,j) << " ";} std::cout <<'\n';}
	
};

template <class Tdf> 
void ResponseMatrix::harmonicDotOmega(const Tdf & df) {
	Eigen::ArrayXXd omega1Grid{df.omega1Grid(m_intStep, m_spacing).array()}, omega2Grid{df.omega2Grid(m_intStep, m_spacing).array()}, 
					elJacobian{df.energyAngMomJacobain(m_intStep, m_spacing).array()}, inverseOmega1{Eigen::ArrayXXd::Zero(m_intStep, m_intStep)},
					dFdE{df.dFdEgrid(m_spacing, omega1Grid).array()}, dFdJgrid{df.dFdJgrid(m_spacing, omega2Grid).array()};
					 
	for (int i = 0; i < m_intStep; ++i){
		for (int j = 1; j < i; ++j) {inverseOmega1(i,j) = 1/omega1Grid(i,j);}
	}

	for (int m1 = - m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++m1) {
		v_harmonicDotOmega[m1 + m_maxFourierHarmonic] = (m1 * omega1Grid + m_fourierHarmonic * omega2Grid);
		v_gPrefactor[m1 + m_maxFourierHarmonic] = pow(2 * M_PI, 2) * elJacobian * (inverseOmega1) *  (m1 * dFdE + m_fourierHarmonic * dFdJgrid);

	} 
}


/* calculating the Response matrix */ 

Eigen::MatrixXcd ResponseMatrix::responseMatrix(const std::complex<double> & omega) {
	//std::cout << "Calculating response matrix for: " << omega << '\n';
	m_responseMatrix = Eigen::MatrixXcd(m_BF.maxRadialIndex()+1, m_BF.maxRadialIndex()+1); 
	for (int p = 0; p < m_responseMatrix.rows(); ++p) {
		for (int q = 0; q <= p; ++q) {
			responseMatrixElement(omega, p, q);
		}
	}
	return m_responseMatrix; 
}


void ResponseMatrix::responseMatrixElement(const std::complex<double> & omega, int p, int q) {
	std::complex<double> element(0,0);
	for (int m1 = - m_maxFourierHarmonic; m1 <= m_maxFourierHarmonic; ++m1) 
	{
		calculateNormalisedCoeff(omega, m1, p, q); 
		element += alephSummedOver(); 
	}  
	m_responseMatrix(p,q) = element;
	m_responseMatrix(q,p) = element;
	//std::cout << element <<'\n'; 
}
void ResponseMatrix::calculateNormalisedCoeff(const std::complex<double> & omega, int m1, int p, int q) { 
	v_taylorCoeff[0] = v_gPrefactor[m1 + m_maxFourierHarmonic] * m_BF(p, m1).array()* m_BF(q, m1).array(); 
	v_taylorCoeff[1] = omega.real() - v_harmonicDotOmega[m1 + m_maxFourierHarmonic]; 
	
	v_taylorCoeff[6] = Eigen::ArrayXXd::Constant(m_intStep, m_intStep, omega.imag());

	getDerivCoeff(); 
	normaliseCoeff(); 
}

void ResponseMatrix::getDerivCoeff() {
	
	for (int i = 4; i < m_intStep-1; i += 2)
	{
		for (int j = 2; j < i; j += 2)  
		{ // Note that \Delta r/2 = m_spacing
 			v_taylorCoeff[2](i,j) = (0.25/m_spacing) * (v_taylorCoeff[0](i , j+1) - v_taylorCoeff[0](i , j-1)); // bg
			v_taylorCoeff[3](i,j) = (0.25/m_spacing) * (v_taylorCoeff[1](i , j+1) - v_taylorCoeff[1](i , j-1)); // bh

			v_taylorCoeff[4](i,j) = (0.25/m_spacing) * (v_taylorCoeff[0](i-1 , j) - v_taylorCoeff[0](i+1 , j)); // cg
			v_taylorCoeff[5](i,j) = (0.25/m_spacing) * (v_taylorCoeff[1](i-1 , j) - v_taylorCoeff[1](i+1 , j)); // ch
		}
		
	} 
}

void ResponseMatrix::normaliseCoeff() {
	normaliseTwoArrays(2,0);
	normaliseTwoArrays(3,1);

	normaliseTwoArrays(4,0);
	normaliseTwoArrays(5,1);

	normaliseTwoArrays(6,1);
}

void ResponseMatrix::normaliseTwoArrays(int n, int d) {
	for (int i = 4; i < m_intStep-1; i += 2)
	{
		for (int j = 2; j < i; j += 2)  { v_taylorCoeff[n](i,j) /= v_taylorCoeff[d](i,j);}
	}
}

std::complex<double> ResponseMatrix::alephSummedOver() const {
	std::complex<double> summedAleph(0,0); 
	for (int i = 4; i < m_intStep-1; i += 2)
	{
		for (int j = 2; j < i; j += 2) {
			summedAleph += (v_taylorCoeff[0](i,j) / v_taylorCoeff[1](i,j)) * pow(2 * m_spacing, 2) * Gfunctions(i, j);
			//dumpCoeff(i,j); 
		}
	} 
	
	return summedAleph; 
}

std::complex<double> ResponseMatrix::Gfunctions(int i, int j) const { // We can check this function 
	std::complex<double> bf{v_taylorCoeff[2](i,j)}, cf{v_taylorCoeff[4](i,j)}, 
	ef{v_taylorCoeff[3](i,j)}, ff{v_taylorCoeff[5](i,j)}, 
	etaf{v_taylorCoeff[6](i,j)}, I(0,1);


	return 1.0/(16.0*ef*ef*ff*ff)*
	(
	-I * (ff*(bf*(ef*ef - (-2.0 + ff - 2.0*I * etaf)*(-2.0 + ff - 2.0*I * etaf)) + 4.0*ef*(2.0 + ef + 2.0*I * etaf)) - cf*ef*(2.0 + ef + 2.0*I * etaf)*(2.0 + ef + 2.0*I * etaf)) * 
    (M_PI - 2.0*atan((2.0 + ef - ff)/(2.0*etaf))) - 
   
   I * (-ff*(bf*(-ef*ef + (2.0 + ff + 2.0*I * etaf)*(2.0 + ff + 2.0*I * etaf)) + 4.0*ef*(-2.0 + ef - 2.0*I * etaf)) - cf*ef*(-2.0 + ef - 2.0*I * etaf)*(-2.0 + ef - 2.0*I * etaf)) *
	(M_PI - 2.0*atan((2.0 - ef + ff)/(2.0*etaf))) + 
   
   I * (ff*(bf*(ef*ef - (-2.0 + ff - 2.0*I * etaf)*(-2.0 + ff - 2.0*I * etaf)) - 4.0*ef*(-2.0 + ef - 2.0*I * etaf)) - cf*ef*(-2.0 + ef - 2.0*I * etaf)*(-2.0 + ef - 2.0*I * etaf)) *
   (M_PI + 2.0*atan((-2.0 + ef + ff)/(2.0*etaf))) + 
   
   I * (ff*(bf*(ef*ef - (2.0 + ff + 2.0*I * etaf)*(2.0 + ff + 2.0*I * etaf)) +  4.0*ef*(2.0 + ef + 2.0*I * etaf)) - cf*ef*(2.0 + ef + 2.0*I * etaf)*(2.0 + ef + 2.0*I * etaf)) *
    (M_PI - 2.0*atan((2.0 + ef + ff)/(2.0*etaf))) + 
   
   ff*(ff*(8.0*ef + bf*(-4.0+2.0*ef + ff - 4.0*I * etaf)) + cf*ef*(-4.0+2.0*ef - ff - 4.0*I * etaf) + 
    2.0*(-4.0+cf)*ef*ff*log(0.5*(2.0 - ef - ff + 2.0*I * etaf))) - 
   
   ff*(ff*(8.0*ef + bf*(-4.0 - 2.0*ef + ff - 4.0*I * etaf)) - cf*ef*(4.0+2.0*ef + ff + 4.0*I * etaf) + 
    2.0*(-4.0+cf)*ef*ff*log(0.5*(2.0 + ef - ff + 2.0*I * etaf))) - 
   
   ff*(ff*(-8.0*ef + bf*(4.0 - 2.0*ef + ff + 4.0*I * etaf)) - cf*ef*(-4.0+2.0*ef + ff - 4.0*I * etaf) + 
   2.0*(4.0+cf)*ef*ff*log(0.5*(2.0 - ef + ff + 2.0*I * etaf))) + 
   
   ff*(ff*(-8.0*ef + bf*(4.0+2.0*ef + ff + 4.0*I * etaf)) + cf*ef*(4.0+2.0*ef - ff + 4.0*I * etaf) + 
   2.0*(4.0+cf)*ef*ff*log(0.5*(2.0 + ef + ff + 2.0*I * etaf))) + 


   (ff*(bf*(-ef*ef + (-2.0 + ff - 2.0*I * etaf)*(-2.0 + ff - 2.0*I * etaf)) - 4.0*ef*(2.0 + ef + 2.0*I * etaf)) + cf*ef*(2.0 + ef + 2.0*I * etaf)*(2.0 + ef + 2.0*I * etaf)) *
   log(0.25*(4.0+ef*ef - 2.0*ef*(-2.0 + ff) - 4.0*ff + ff*ff + 4.0*etaf*etaf)) - 

   (ff*(bf*(-ef*ef + (-2.0 + ff - 2.0*I * etaf)*(-2.0 + ff - 2.0*I * etaf)) + 4.0*ef*(-2.0 + ef - 2.0*I * etaf)) + cf*ef*(-2.0 + ef - 2.0*I * etaf)*(-2.0 + ef - 2.0*I * etaf)) * 
    log(0.25* (4.0+ef*ef + 2.0*ef*(-2.0 + ff) - 4.0*ff + ff*ff + 4.0*etaf*etaf)) + 

	(ff*(bf*(-ef*ef + (2.0 + ff + 2.0*I * etaf)*(2.0 + ff + 2.0*I * etaf)) + 4.0*ef*(-2.0 + ef - 2.0*I * etaf)) + cf*ef*(-2.0 + ef - 2.0*I * etaf)*(-2.0 + ef - 2.0*I * etaf)) *
    log(0.25* (4.0+ef*ef + 4.0*ff + ff*ff - 2.0*ef*(2.0 + ff) + 4.0*etaf*etaf)) - 

    (ff*(bf*(-ef*ef + (2.0 + ff + 2.0*I * etaf)*(2.0 + ff + 2.0*I * etaf)) - 4.0*ef*(2.0 + ef + 2.0*I * etaf)) + cf*ef*(2.0 + ef + 2.0*I * etaf)*(2.0 + ef + 2.0*I * etaf)) *
    log(0.25*(4.0+ef*ef + 4.0*ff + ff*ff + 2.0*ef*(2.0 + ff) + 4.0*etaf*etaf))

    );
}

/* Calculation Functions */

std::complex<double> ResponseMatrix::det() const {
	Eigen::MatrixXcd id = Eigen::MatrixXcd::Identity(m_responseMatrix.rows(), m_responseMatrix.rows());
	return (m_responseMatrix - id).determinant(); 
}

std::string signC(double number) {if (number < 0) {return std::to_string(number) + "j";} else {return"+" + std::to_string(number) + "j";}}
void ResponseMatrix::saveResponseMatrix(const std::string & filename) const {
	std::ofstream out(filename);

	//auto sign = [] (double number) {if (number < 0) {return std::to_string(number);} else {return"+" + std::to_string(number);}};

	for (int i = 0; i < m_responseMatrix.rows(); ++i) {
		for (int j = 0; j < m_responseMatrix.cols(); ++j) {
			std::complex<double> n{m_responseMatrix(i,j)}; 
			if (i == (m_responseMatrix.rows()-1) && j == (m_responseMatrix.cols()-1)) {out << std::to_string(n.real()) + signC(n.imag());}

			else if (j == (m_responseMatrix.cols()-1)) {out << std::to_string(n.real()) + signC(n.imag()) <<'\n';}

			else {out << std::to_string(n.real()) + signC(n.imag()) + ',';}
		}
	}


	out.close(); 

}




#endif
