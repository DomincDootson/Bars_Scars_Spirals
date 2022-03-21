#ifndef DFCLASS
#define DFCLASS 

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <random>
#include <fstream>

#include <Eigen/Dense>


class DFClass
{
public:
	DFClass() {}
	~DFClass() {}

	virtual double potential(double radius) const = 0;
	virtual double distFunc(double E, double J) const = 0;

	double rad2Energy(double rApo, double rPer) const;
	double rad2AngMom(double rApo, double rPer) const;

	double omega1(double rApo, double rPer) const;
	double omega2(double rApo, double rPer) const;

	Eigen::MatrixXd omega1Grid(int dim, double spacing) const;
	Eigen::MatrixXd omega2Grid(int dim, double spacing) const;


	double theta1(double radius, double rApo, double rPer, double omega1) const;
	double theta2(double radius, double rApo, double rPer, double omega2) const; // Note this is really theta2 - phi 
	double theta1Deriv(double radius, double rApo, double rPer, double omega1) const;

	Eigen::MatrixXd dFdEgrid(const double spacing, const Eigen::MatrixXd & om1) const; // Note that this includes om1 
	Eigen::MatrixXd dFdJgrid(const double spacing, const Eigen::MatrixXd & om2) const; // Note that this includes om2

	Eigen::MatrixXd energyAngMomJacobain(const int size, const double spacing) const;

	double density(const double radius) const;
	double maxDensityRadius() const; 
	void cumulativeDensity(const std::string fileName) const;
	double diskMass() const; 

	void savingDensity(const std::string fileName, const double spacing = 0.03, const int nSteps = 1000);
	void densityVector(const int nSteps, const double spacing);
	
	double radiusSampling(const double uniformDistNumb, const std::vector<double> & cumulativeDensity) const; 
	std::vector<double> readInCumulativeDensity(const std::string & fileName) const;

	double vPhiSampling(const double radius, const double vr) const;//  I want to put this in the DF_Class
	virtual double vRSampling() const = 0; // I WANT THIS TO BE VIRTUAL, BUT IT IS NOT WORKING
	
	double xAccel(const double xPos, const double yPos) const;
	double yAccel(const double xPos, const double yPos) const;


	double densityFromVector(const double radius) const; 
	double maxDensityFromVector() const { 
    auto isGreater = [](auto x, auto y) {return x < y ? true : isnan(x);}; return *std::max_element(v_densityVector.begin(),v_densityVector.end(), isGreater);}
	void readInDensityVector(const std::string & filename);  
	Eigen::ArrayXXcd densityArray(const double rMax, const int nGrid) const;
protected: 	
	std::vector<double> v_densityVector;
	double m_densitySpacing; 
	
	virtual double vRScale() const = 0; // Give the upper limit for marganlising over
	virtual double vPhiScale() const = 0;
	
	virtual double jMax(const double radius) const = 0;
	virtual double dfMax(const double radius, const double vR) const = 0;



};

double DFClass::rad2Energy(double rApo, double rPer) const 
{
	assert(rPer <= rApo);
	if (rPer == rApo)
	{
		return 0.5 + potential(rPer);
	}
	return (pow(rApo,2)*potential(rApo) -pow(rPer,2)*potential(rPer)) / (rApo*rApo - rPer*rPer);
}

double DFClass::rad2AngMom(double rApo, double rPer) const
{
	assert(rPer <= rApo);
	if (rPer == rApo)
	{
		return rPer;
	}
	return pow(2*(potential(rApo) - potential(rPer))/ (pow(rPer ,-2) - pow(rApo ,-2)), 0.5);	
}

// Function for linear calculation //
// ------------------------------- //


double DFClass::omega1(double rApo, double rPer) const // what to do for circular orbits?
{
	double E{rad2Energy(rApo, rPer)}, J{rad2AngMom(rApo, rPer)};
	double integral{0}, upperU{0.5 * M_PI}, nstep{1000}, stepSize{upperU/(nstep -1)}, rstep{};
	for (int i = 1; i < nstep-1; ++i)
	{
		rstep = rPer +(rApo-rPer)*pow(sin(stepSize * i),2);
		integral += sin(2 * i * stepSize) * pow(2*(E- potential(rstep)) -J*J * pow(rstep, -2) , -0.5) * stepSize; 
	}

	return M_PI / ((rApo - rPer) * integral);
}


double DFClass::omega2(double rApo, double rPer) const // what to do for circular orbits?
{
	double E{rad2Energy(rApo, rPer)}, J{rad2AngMom(rApo, rPer)};
	double integral{0}, upperU{0.5 * M_PI}, nstep{1000}, stepSize{upperU/(nstep-1)}, rstep{};
	for (int i = 1; i < nstep-1; ++i)
	{
		rstep = rPer +(rApo-rPer)*pow(sin(stepSize * i),2);
		integral += sin(2 * i * stepSize) * pow(2*(E- potential(rstep)) -J*J * pow(rstep, -2) , -0.5) * (stepSize / (rstep*rstep)); 	
	}
	
	return DFClass::omega1(rApo, rPer) * J * ((rApo - rPer) * integral)/M_PI;
}
  
Eigen::MatrixXd DFClass::omega1Grid(int dim, double spacing) const {
	Eigen::MatrixXd grid = Eigen::MatrixXd::Zero(dim, dim);
	for (int i = 0; i < dim; ++i){
			for (int j = 1; j < i; ++j){double rApo{i*spacing}, rPer{j*spacing}; grid(i,j) = omega1(rApo, rPer);}}
	return grid;
}

Eigen::MatrixXd DFClass::omega2Grid(int dim, double spacing) const {
	Eigen::MatrixXd grid = Eigen::MatrixXd::Zero(dim, dim);
	for (int i = 0; i < dim; ++i){
			for (int j = 1; j < i; ++j){double rApo{i*spacing}, rPer{j*spacing}; grid(i,j) = omega2(rApo, rPer);}}
	return grid;
}


double DFClass::theta1(double radius, double rApo, double rPer, double omega1) const
{
	int nstep{1000};
	double E{rad2Energy(rApo, rPer)}, J{rad2AngMom(rApo, rPer)};
	double integral{0}, upperU{asin(pow((radius- rPer)/(rApo-rPer), 0.5))}, stepSize{upperU/(nstep-1)}, rstep{};
	
	for (int i = 1; i < nstep; ++i)
	{
		rstep = rPer +(rApo-rPer)*pow(sin(stepSize *i),2);
		integral += sin(2 * i * stepSize)*pow(2*(E- potential(rstep)) -J*J * pow(rstep, -2) , -0.5)*stepSize;// *  
	}
	double th1 = (rApo - rPer) * integral * omega1;

	if (th1 != th1) {th1 = 0;} // For very small integrals you get numberical error that cause infinites. 
	return th1;
}

double DFClass::theta2(double radius, double rApo, double rPer, double omega2) const 
{
	int nstep{1000};
	double E{rad2Energy(rApo, rPer)}, J{rad2AngMom(rApo, rPer)};
	double integral{0}, upperU{asin(pow((radius- rPer)/(rApo-rPer), 0.5))}, stepSize{upperU/(nstep-1)}, rstep{};

	for (int i = 1; i < nstep; ++i)
	{
		rstep = rPer +(rApo-rPer)*pow(sin(stepSize *i),2);
		integral += (omega2 - J*pow(rstep,-2))* sin(2 * i * stepSize) * pow(2*(E- potential(rstep)) -J*J * pow(rstep, -2) , -0.5) * stepSize; 
	}

	if (integral != integral) {return 0;}

	return  (rApo - rPer) * integral;
}

double DFClass::theta1Deriv(double radius, double rApo, double rPer, double omega1) const
{
	double E{rad2Energy(rApo, rPer)}, J{rad2AngMom(rApo, rPer)};
	return omega1*pow(2 * (E - potential(radius)) - pow(J/radius, 2), -.5);  
}

Eigen::MatrixXd DFClass::dFdEgrid(const double spacing, const Eigen::MatrixXd & om1) const
{
	Eigen::MatrixXd dFdE{Eigen::MatrixXd::Zero(om1.rows(), om1.cols())};
	double E{}, J{};
	for (int i = 1; i < dFdE.rows(); ++i)
	{
		for (int j = 1; j < i; ++j)
		{
			E = rad2Energy(i * spacing,  j * spacing); 
			J = rad2AngMom(i * spacing,  j * spacing); 
			dFdE(i,j) = 0.5  * ((distFunc(E+0.001, J) - distFunc(E-0.001,J))/0.001)* om1(i,j);		
		}
	}
	return dFdE;
}


Eigen::MatrixXd DFClass::dFdJgrid(const double spacing, const Eigen::MatrixXd & om2) const
{
	Eigen::MatrixXd dFdJ{Eigen::MatrixXd::Zero(om2.rows(), om2.cols())};
	double E{}, J{};

	for (int i = 1; i < dFdJ.rows(); ++i)
	{
		for (int j = 1; j < i; ++j)
		{
			E = rad2Energy(i * spacing,  j * spacing); 
			J = rad2AngMom(i * spacing,  j * spacing); 
			dFdJ(i,j) =
			om2(i,j) * (distFunc(E+0.00001, J) - distFunc(E-0.00001,J)) / (2*0.00001) 
						+ (distFunc(E, J+0.00001) - distFunc(E,J-0.00001)) / (2*0.00001);

		}
	}
	return dFdJ;
}

Eigen::MatrixXd DFClass::energyAngMomJacobain(const int size, const double spacing) const
{
	Eigen::MatrixXd jacobian{Eigen::MatrixXd::Zero(size, size)};
	double step{0.001}, rApo{}, rPer{}, dErApo{}, dErPer{}, dLrApo{}, dLrPer{};
	for (int i = 1; i < size; ++i)
	{
		for (int j = 1; j < i; ++j)
		{
			rApo = i * spacing;
			rPer = j * spacing; 
			
			dErApo = (rad2Energy(rApo + step, rPer) - rad2Energy(rApo - step, rPer))/(2*step);
			dErPer = (rad2Energy(rApo, rPer + step) - rad2Energy(rApo, rPer - step))/(2*step);

			dLrApo = (rad2AngMom(rApo + step, rPer) - rad2AngMom(rApo - step, rPer))/(2*step);
			dLrPer = (rad2AngMom(rApo, rPer + step) - rad2AngMom(rApo, rPer - step))/(2*step);

			jacobian(i,j) = abs((dErApo)*(dLrPer) - (dErPer)*(dLrApo));
		}
	}
	return jacobian;
}

// Functions for n Body //
// -------------------- //

double DFClass::density(const double radius) const // Maybe include some integration constants
{
	int nSteps{300}, nStepsRadial{300}; 
	double density{0}, stepVPHI{vPhiScale()/nSteps}, stepVR{vRScale()/nSteps}, energy{0}, vPhi{0}, vrIntegral{}; // FIGURE OUT THE STEP SIZE
	//0.26112 // 5.49443e+14
	for (int i = 1; i<nSteps; ++i) // Outer integral is over J (don't start at J = 0, as f(J) = 0)
	{
		vPhi = (i * stepVPHI);
		vrIntegral = 0;
		for (int k = 0; k < nStepsRadial; ++k) // Remember df is even in vr. 
		{
			double df = distFunc(energy, vPhi*radius);
			energy = 0.5* (vPhi*vPhi) + 0.5*(k*stepVR)*(k*stepVR) + potential(radius);
			vrIntegral += 2*stepVR*df;

		}

		
		density += stepVPHI*vrIntegral;
	}
	return density;
}



void DFClass::cumulativeDensity(const std::string fileName) const 
{
	std::vector<double> densityArray, cumulative;
	double spacing{.1};

	densityArray.push_back(0); cumulative.push_back(0); int n{1};
	do
	{
		densityArray.push_back(2*M_PI * (spacing* n) * density(n*spacing));	
		cumulative.push_back(spacing * densityArray.back() + cumulative.back());
		n += 1;
		//std::cout << n * spacing << " " << densityArray.back() << '\n';
		if ( (n * spacing) ==static_cast<int>(n * spacing)) {std::cout << "Finding cumulative density for point: " << n*spacing << " Density: " <<  cumulative.back()<< '\n';}
	} while (n*spacing < 10 || n*spacing <20); //|| (densityArray.back()/cumulative.back() > .001)


	std::ofstream out(fileName);
	out << n << " " << spacing << " " << cumulative.back() <<'\n';
	for (int i = 0; i < n; ++i)	{out << i * spacing << " " << cumulative[i]/cumulative.back() << '\n';}
	out.close();
}

double DFClass::diskMass() const {
	double mass{0}, spacing{0.1}; std::cout << "Calculating Disk Mass\n";
	for (double r = spacing; r < 25; r += spacing) {mass += 2*M_PI * (r) * density(r) * spacing;}
	return mass; 
}

std::vector<double> DFClass::readInCumulativeDensity(const std::string & fileName) const{ // I am not sure why this must be a memeber function.... 
	std::ifstream inFile; inFile.open(fileName);
	int n; double spacing, normalisation;
	inFile >> n >> spacing >> normalisation;

	std::vector<double> cumulative = {spacing};
	for (int i = 0; i < n; ++i){
		double radius, value;
		inFile >> radius >>value;
		cumulative.push_back(value);
	}
	inFile.close();
	return cumulative;
}

double DFClass::maxDensityRadius() const // Assume that there is only one turning point 
{
	std::vector<double> densityArray = {0};
	double stepSize{0.1}; int size{(int) densityArray.size()};
	do
	{
		densityArray.push_back(density(size*stepSize));
		size = densityArray.size();
	} while (densityArray[size-1] >  densityArray[size-2]);
	
	double gradBeforeTP{abs(densityArray[size-2] - densityArray[size-3])}, gradAfterTP{abs(densityArray[size-1] - densityArray[size-2])};
	return stepSize * (size-1.5) + stepSize * (gradBeforeTP)/(gradBeforeTP + gradAfterTP); // Linearly interpolate
}


double DFClass::radiusSampling(const double uniformDistNumb, const std::vector<double> & cumulativeDensity) const 
{
	int nLower, nUpper;
	double spacing = cumulativeDensity.front();
	for (int i = 1; i < cumulativeDensity.size(); ++i)
	{
		if ((cumulativeDensity[i] < uniformDistNumb) && (cumulativeDensity[i+1] > uniformDistNumb))
		{nLower = i; nUpper = i+1;
		return (nLower + (uniformDistNumb - cumulativeDensity[nLower])/(cumulativeDensity[nUpper] - cumulativeDensity[nLower]))*spacing; // Put into the loop to kill it once we have the required values
		}
	}
}

// Some virtual function that does the vr sampling 

double DFClass::vPhiSampling(const double radius, const double vr) const 
{
	
	double J, f, E, maxJ{jMax(radius)}, maxF{dfMax(radius, vr)};
	
	std::random_device generator;
	std::uniform_real_distribution<double> uniform(0,1);

	do // Rejection sampling
	{
		J = uniform(generator) * maxJ;
		f = uniform(generator) * maxF;

		E = 0.5*vr*vr + 0.5*(J/radius)*(J/radius) + potential(radius);
	} while (f > distFunc(E,J)); 
	
	return (J/radius);
	
}


// Some accel functions for the nnBody

double DFClass::xAccel(const double xPos, const double yPos) const{
	double spacing{0.001}, rInner{sqrt(pow(xPos-spacing ,2) + yPos*yPos)}, rOuter{sqrt(pow(xPos+spacing ,2) + yPos*yPos)};
	return -(1/(2*spacing)) * (potential(rOuter) - potential(rInner));
}

double DFClass::yAccel(const double xPos, const double yPos) const{
	double spacing{0.001}, rInner{sqrt(pow(yPos-spacing ,2) + xPos*xPos)}, rOuter{sqrt(pow(yPos+spacing ,2) + xPos*xPos)};
	return -(1/(2*spacing)) * (potential(rOuter) - potential(rInner));
}

/* Density Grid Function */
/* --------------------- */

void DFClass::readInDensityVector(const std::string & filename) {
	std::ifstream inFile(filename); if (!inFile.good()) {std::cout<< "Couldn't find file: " << filename <<'\n'; exit(0);}
	int n{};
	inFile >>n >> m_densitySpacing; v_densityVector.resize(n);
	for (int i =0; i<n; ++i) {inFile >> v_densityVector[i];}

	inFile.close();
}

void DFClass::densityVector(const int nSteps, const double spacing) {
	v_densityVector.resize(nSteps); m_densitySpacing = spacing;
	for (int i =0; i < nSteps; ++i) {v_densityVector[i] = density(i *spacing);}
}

void DFClass::savingDensity(const std::string fileName, const double spacing, const int nSteps) {
 	std::ofstream out(fileName); std::cout << "Saving density to: " << fileName << '\n';
	densityVector(nSteps, spacing);

	out << v_densityVector.size() << " " << m_densitySpacing << '\n';
	out << 0 << '\n';

	for (int i = 1; i < nSteps; i++) {out << v_densityVector[i] << '\n'; if(i%1==0) {std::cout << "Density for: " << i*spacing << '\n';}}
	out.close();
}

double DFClass::densityFromVector(const double radius) const {
	double index{radius/m_densitySpacing}; int lowerIndex{(int) floor(index)}; 
	if (lowerIndex+1 >= v_densityVector.size()) {return 0;}
	return v_densityVector[lowerIndex] + (index-lowerIndex) *(v_densityVector[lowerIndex+1] - v_densityVector[lowerIndex]); 
}

Eigen::ArrayXXcd DFClass::densityArray(const double rMax, const int nGrid) const {
	Eigen::ArrayXXcd densityArray = Eigen::ArrayXXcd::Zero(nGrid, nGrid);
	double x{}, y{}, r{}, spacing{2*rMax/((double) nGrid-1)}, centre{0.5*(nGrid-1)};
	for (int i = 0; i < densityArray.rows(); ++i)
	{
		for (int j = 0; j < densityArray.cols(); ++j)
		{
			x = spacing * (i - centre); y = spacing * (j - centre);
			r = sqrt(x*x+y*y); 
			if (r!= 0){
				densityArray(i,j) = densityFromVector(r);
			}
		}
	}
	return densityArray;
}

#endif