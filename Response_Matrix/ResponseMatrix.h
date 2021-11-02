#ifndef RESPONSEMATRIX
#define RESPONSEMATRIX

#include <ResponseMatrixElement.h>

class ResponseMatrix
{
public:
	ResponseMatrix() {}
	~ResponseMatrix() {}

	Eigen::MatrixXcd responseMatrix(const double omega);

private:
	// Some internal parameters
	ResponseMatrixElement m_elementCalculator;

	Eigen::MatrixXcd responseMatrixM1(const double omega, const int m1);
	
};

#endif