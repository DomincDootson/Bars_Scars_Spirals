#ifndef THEORETICALPDCONTAINER
#define THEORETICALPDCONTAINER

#include <vector>


template <class T> 
class TheoreticalPDContainer : public PotentialDensityPairContainer<T> {
public: 
	TheoreticalPDContainer(const std::vector<double> & params, const int maxN, const int l) : PotentialDensityPairContainer<T>{params, maxN, l}
	{}
	
	~TheoreticalPDContainer() {}

	void analyticKernel(const std::string & filename, const int nTimeStep, const double timeStep, const double sigmaR) const; 

private:  
	Eigen::MatrixXcd analyticKernel(const int nCols, const double timeDelta, const double sigmaR) const; 
	void writeKernels2File(const std::string & filename, const std::vector<Eigen::MatrixXcd> & kernels, const double timeStep) const;
};


template <class T>
Eigen::MatrixXcd TheoreticalPDContainer<T>::analyticKernel(const int nCols, const double timeDelta, const double sigmaR ) const {
	Eigen::MatrixXcd holding = Eigen::MatrixXcd::Zero(nCols, nCols);
	for (int i = 0; i < holding.cols(); ++i) {holding(i,i) = this->m_potentialDensityContainer[i].analyticKernel(timeDelta, sigmaR);}
	return holding; 
}

template <class T>
void TheoreticalPDContainer<T>::writeKernels2File(const std::string & filename, const std::vector<Eigen::MatrixXcd> & kernels, const double timeStep) const {
	std::ofstream out(filename); int nCols = kernels[0].cols();
	out << nCols-1 << " " << this->m_fourierHarmonic << " " << kernels.size() << " " << 0 << " " << timeStep << " " << 0 << '\n';
	for (int time = 0; time < kernels.size(); ++ time){
		for (int i = 0; i< nCols; ++i){
			if (i == (nCols-1)) {out << real(kernels[time](i,i)) << " " << imag(kernels[time](i,i)) <<'\n';}
			else {out << real(kernels[time](i,i)) << " " << imag(kernels[time](i,i)) << " ";}
		}
	}
	out.close();
	std::cout << "Kernel saved to: " << filename << '\n';
}

template <class T>
void TheoreticalPDContainer<T>::analyticKernel(const std::string & filename, const int nTimeStep, const double timeStep, const double sigmaR) const {
	std::vector<Eigen::MatrixXcd> kernels(nTimeStep);
	for (int timeIndex = 0; timeIndex < nTimeStep; timeIndex++) {kernels[timeIndex] = analyticKernel(this->m_maxRadialIndex+1, timeIndex * timeStep, sigmaR);}
	writeKernels2File(filename, kernels, timeStep);
}



#endif
