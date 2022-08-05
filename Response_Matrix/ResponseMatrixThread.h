#ifndef RESPONSEMATRIXTHREAD
#define RESPONSEMATRIXTHREAD 

#include <vector>
#include <thread>
#include <complex>

#include "ResponseMatrix.h"
#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"

using ComplexVector = std::vector<std::complex<double>>; // We really want to define this in local scope 

class ResponseMatrixThread
{
public:

	template <class Tdf>
	ResponseMatrixThread(const ActionAngleBasisContainer & BF, const Tdf & df, int nCores) : m_cores{nCores}
	{
		/*std::vector<std::thread> threads; // PLEASE FIGURE OUT MORE ABOUT THREAD SAFTEY
		for (int i = 0; i < m_cores; ++i) {threads.emplace_back( [this, BF, df] () {constructRM(BF, df);});}
		for (int i = 0; i < m_cores; ++i) {std::cout << threads[i].get_id() <<'\n';}
		

		for (auto& thread : threads) {thread.join();}*/

		for (int i = 0; i < m_cores; ++i) {v_responseMatrices.emplace_back(BF, df);}
	}
	
	~ResponseMatrixThread() {}

	ComplexVector det(const ComplexVector & omegas);

	void modeGridSearch(const std::string & filename, double eMin, double eMax, int eN, double oMin, double oMax, int oN); 


private:
	int m_cores; 
	std::vector<ResponseMatrix> v_responseMatrices; 

	template <class Tdf>
	void constructRM(const ActionAngleBasisContainer & BF, const Tdf & df) {std::cout << "Constructing RM\n"; v_responseMatrices.emplace_back(BF, df);}

	ComplexVector uniformOmegaVector(double eMin, double eMax, int eN, double oMin, double oMax, int oN);
	void saveVector(std::ofstream & out, const ComplexVector & vec, int eN, int oN ) const;
};

ComplexVector ResponseMatrixThread::det(const ComplexVector& omegas) {
	ComplexVector dets; 
	for (int i = 0; i < omegas.size(); i += m_cores) {
		std::vector<std::thread> threads; int nUsed{m_cores}; // Incase omeags.size()//m_cores !=0
		std::cout << "Percentage of modes calculated: " << int(100*i / ((double) omegas.size())) << '%' << '\n';
		for (int c = 0; c < m_cores; ++c) {
			if (i+c >= omegas.size()) {nUsed = c+1; break;}
			std::complex<double> omega{omegas[i+c]};
			threads.emplace_back( [this, omega, c] () {v_responseMatrices[c].responseMatrix(omega);}); 
		}

		for (auto & th : threads) {th.join();}
		for (int n = 0; n < nUsed; ++n) {dets.emplace_back(v_responseMatrices[n].det());}  
	}
	return dets; 
}

/* modeGridSearch Functions */ 
/* ------------------------ */
std::string complex2String(const std::complex<double> & cn) {
	if (cn.imag() >= 0) {return std::to_string(cn.real()) + '+' + std::to_string(cn.imag()) +'j';}
	else {return std::to_string(cn.real()) + std::to_string(cn.imag()) +'j';}
}

void ResponseMatrixThread::modeGridSearch(const std::string & filename, double eMin, double eMax, int eN, double oMin, double oMax, int oN) {
	ComplexVector omegaVec{uniformOmegaVector(eMin, eMax, eN, oMin, oMax, oN)}, detVec(eN*oN);
	detVec = det(omegaVec); 

	std::ofstream out(filename);
	saveVector(out, omegaVec, eN, oN);
	out << '\n' <<'\n'; // Put space in so we know where the grid stops and det starts  
	saveVector(out, detVec, eN, oN);
	out.close(); 
	std::cout << "Written determinants to: " << filename << '\n';
}

ComplexVector ResponseMatrixThread::uniformOmegaVector(double eMin, double eMax, int eN, double oMin, double oMax, int oN) {
	ComplexVector holding; std::complex<double> unitComplex(0,1);
	double oSpacing{(oMax-oMin)/ ((double) oN-1)}, eSpacing{(eMax-eMin)/ ((double) eN-1)}; 
	for (int i = 0; i < eN; ++i) {
		for (int j = 0; j < oN; ++j) {
			holding.emplace_back((j * oSpacing + oMin) + unitComplex * (i * eSpacing +eMin)); // Access with [j + i * oN]
		}
	}

	return holding;
}


void ResponseMatrixThread::saveVector(std::ofstream & out, const ComplexVector & vec, int eN, int oN ) const {
	for (int i = 0; i < eN; ++i) {
		for (int j = 0; j < oN-1; ++j) {
			out << complex2String(vec[j + i * oN]) <<',';
		}
		if (i < (eN-1)) {out << complex2String(vec[oN-1 + i * oN]) <<'\n';}
		else {out << complex2String(vec[oN-1 + i * oN]);}
	}	

}


#endif