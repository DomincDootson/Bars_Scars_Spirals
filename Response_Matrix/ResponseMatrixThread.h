#ifndef RESPONSEMATRIXTHREAD
#define RESPONSEMATRIXTHREAD 

#include <vector>
#include <thread>
#include <algorithm>
#include <complex>

#include "ResponseMatrix.h"
#include "../Action_Angle_Basis_Functions/ActionAngleBasisContainer.h"
#include "../Volterra_Solver/EvolutionKernels.h"


using ComplexVector = std::vector<std::complex<double>>; // We really want to define this in local scope 

class ResponseMatrixThread
{
public:

	template <class Tdf>
	ResponseMatrixThread(const ActionAngleBasisContainer & BF, const Tdf & df, const int nCores) : m_cores{nCores}
	{
		/*std::vector<std::thread> threads; // PLEASE FIGURE OUT MORE ABOUT THREAD SAFTEY
		for (int i = 0; i < m_cores; ++i) {threads.emplace_back( [this, BF, df] () {constructRM(BF, df);});}
		for (int i = 0; i < m_cores; ++i) {std::cout << threads[i].get_id() <<'\n';}
		

		for (auto& thread : threads) {thread.join();}*/

		for (int i = 0; i < m_cores; ++i) {v_responseMatrices.emplace_back(BF, df);}
	}

	ResponseMatrixThread(const ActionAngleBasisContainer & BF, const int nCores = 1) : m_cores{1} // This constructor will be used to generate modes via kernel method
	{for (auto i = 0; i < m_cores; ++i) {v_responseMatrices.emplace_back(BF);} }

	ResponseMatrixThread(const int nCores = 1) : m_cores{nCores} // This constructor will be used to generate modes via kernel method
	{for (auto i = 0; i < m_cores; ++i) {v_responseMatrices.emplace_back();} }
		
	~ResponseMatrixThread() {}

	ComplexVector det(const ComplexVector & omegas);
	ComplexVector det(const ComplexVector & omegas, const EvolutionKernels & kernel);
	
	template <class Tdf>
	ComplexVector det(const ComplexVector & omegas, EvolutionKernels & kernel, const Tdf & df);

	void modeGridSearch(const std::string & filename, double eMin, double eMax, int eN, double oMin, double oMax, int oN); // These points should really be in a struct
	template <class Tdf>
	void modeGridSearch(const std::string & filename, EvolutionKernels & kernel, const Tdf & df, double eMin, double eMax, int eN, double oMin, double oMax, int oN); 

	void modeGridSearch(const std::string & filename, const EvolutionKernels & kernel, double eMin, double eMax, int eN, double oMin, double oMax, int oN); 
	


private:
	int m_cores; 
	std::vector<ResponseMatrix> v_responseMatrices; 

	template <class Tdf>
	void constructRM(const ActionAngleBasisContainer & BF, const Tdf & df) {std::cout << "Constructing RM\n"; v_responseMatrices.emplace_back(BF, df);}

	ComplexVector uniformOmegaVector(double eMin, double eMax, int eN, double oMin, double oMax, int oN);
	void saveVector(std::ofstream & out, const ComplexVector & vec, int eN, int oN ) const;
	void searchSave(const std::string & filename, const ComplexVector & omegaVec, const ComplexVector & detVec, const int eN, const int oN) const;
};

ComplexVector ResponseMatrixThread::det(const ComplexVector& omegas) { // CAN WE TIDY UP THIS AND THE FUNCTION BELOW IT? 
	ComplexVector dets(omegas.size()); 
	for (int i = 0; i < omegas.size(); i += m_cores) {
		 std::vector<std::thread> threads; int nUsed{m_cores}; // Incase omeags.size()//m_cores !=0
		std::cout << "Percentage of modes calculated: " << int(100*i / ((double) omegas.size())) << '%' << '\n';
		for (int c = 0; c < m_cores; ++c) {
			if (i+c >= omegas.size()) {nUsed = c+1; break;}
			const std::complex<double> omega{omegas[i+c]};
			threads.emplace_back( [this] (const int core, const std::complex<double> & om) {v_responseMatrices[core].responseMatrix(om);}, c, omega); 
		}

		for (auto & th : threads) {th.join();}
		for (int n = 0; n < nUsed; ++n) {dets[i + n] = v_responseMatrices[n].det();}  
	
	}
	return dets; 
}

ComplexVector ResponseMatrixThread::det(const ComplexVector& omegas, const EvolutionKernels & kernel) { // CAN WE TIDY UP THIS AND THE FUNCTION BELOW IT? 
	ComplexVector dets(omegas.size()); 
	for (int i = 0; i < omegas.size(); i += m_cores) {
		 std::vector<std::thread> threads; int nUsed{m_cores}; // Incase omeags.size()//m_cores !=0
		
		if (i%1000 ==0) {std::cout << "Percentage of modes calculated: " << int(100*i / ((double) omegas.size())) << '%' << '\n';}
	
		for (int c = 0; c < m_cores; ++c) {
			if (i+c >= omegas.size()) {nUsed = c+1; break;}
			const std::complex<double> omega{omegas[i+c]};
			threads.emplace_back( [this, kernel] (const int core, const std::complex<double> & om) {v_responseMatrices[core].responseMatrix(om, kernel);}, c, omega); 
		}

		for (auto & th : threads) {th.join();}
		for (int n = 0; n < nUsed; ++n) {dets[i + n] = v_responseMatrices[n].det();}  
	
	}
	return dets; 
}

template <class Tdf>
ComplexVector ResponseMatrixThread::det(const ComplexVector & omegas, EvolutionKernels & kernel, const Tdf & df) {
	ComplexVector dets; std::cout << "Calculating Modes\n";
	for (int i = 0; i < omegas.size(); i += m_cores) {
		std::cout << "Percentage of modes calculated: " << int(100*i / ((double) omegas.size())) << '%' << '\n';
		
		std::complex<double> omega{omegas[i]}; 
		v_responseMatrices[i].responseMatrix(omega, kernel, df); 
		dets.emplace_back(v_responseMatrices[i].det());
	}
	return dets; 
}


/* modeGridSearch Functions */ 
/* ------------------------ */
std::string complex2String(const std::complex<double> & cn) {
	if (cn.imag() >= 0) {return std::to_string(cn.real()) + '+' + std::to_string(cn.imag()) +'j';}
	else {return std::to_string(cn.real()) + std::to_string(cn.imag()) +'j';}
}

void ResponseMatrixThread::modeGridSearch(const std::string & filename, double eMin, double eMax, int eN, double oMin, double oMax, int oN) { //change the way we have done this (i.e. o first )
	ComplexVector omegaVec{uniformOmegaVector(eMin, eMax, eN, oMin, oMax, oN)}, detVec(eN*oN);
	detVec = det(omegaVec); 
	searchSave(filename, omegaVec, detVec, eN, oN); 
}

template <class Tdf>
void ResponseMatrixThread::modeGridSearch(const std::string & filename, EvolutionKernels & kernel, const Tdf & df,double eMin, double eMax, int eN, double oMin, double oMax, int oN) {
	ComplexVector omegaVec{uniformOmegaVector(eMin, eMax, eN, oMin, oMax, oN)}, detVec(eN*oN);
 	detVec = det(omegaVec, kernel, df); 
	searchSave(filename, omegaVec, detVec, eN, oN);
}


void ResponseMatrixThread::modeGridSearch(const std::string & filename, const EvolutionKernels & kernel, double eMin, double eMax, int eN, double oMin, double oMax, int oN) {
	ComplexVector omegaVec{uniformOmegaVector(eMin, eMax, eN, oMin, oMax, oN)}, detVec(eN*oN);
 	detVec = det(omegaVec, kernel); 
	searchSave(filename, omegaVec, detVec, eN, oN);
}

void ResponseMatrixThread::searchSave(const std::string & filename, const ComplexVector & omegaVec, const ComplexVector & detVec, const int eN, const int oN) const 
{
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
	if (eN ==1) {eSpacing = 0;}
	
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