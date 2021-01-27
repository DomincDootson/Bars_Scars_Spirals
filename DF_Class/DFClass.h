#ifndef DFCLASS
#define DFCLASS 

#include <cmath>
#include <Eigen/Dense>


class DFClass
{
public:
	DFClass() {}
	~DFClass() {}

	double omega1(double rApo, double rPer) const;
	double omega2(double rApo, double rPer) const;
	double theta1(double radius, double rApo, double rPer, double omega1) const;
	double theta2(double radius, double rApo, double rPer, double omega2) const;
	double theta1Deriv(double radius, double rApo, double rPer, double omega1) const;

	Eigen::MatrixXd dFdEgrid(const double spacing, const Eigen::MatrixXd & om1) const;
	Eigen::MatrixXd dFdJgrid(const double spacing, const Eigen::MatrixXd & om2) const; 

	Eigen::MatrixXd energyAngMomJacobain(const int size, const double spacing) const;

protected: 	
	virtual double potential(double radius) const = 0;
	virtual double distFunc(double E, double J) const = 0;

	double rad2Energy(double rApo, double rPer) const 
	{
		assert(rPer <= rApo);
		if (rPer == rApo)
		{
			return 0.5 + potential(rPer); // PROPER FUNCTION THAT USES THIS
		}
		return (pow(rApo,2)*potential(rApo) -pow(rPer,2)*potential(rPer)) / (rApo*rApo - rPer*rPer);
	}

	double rad2AngMom(double rApo, double rPer) const
	{
		assert(rPer <= rApo);
		if (rPer == rApo)
		{
			return rPer; // WE NEED TO INCLUDE A PROPER FUNCITON THAT DOES THIS
		}
		return pow(2*(potential(rApo) - potential(rPer))/ (pow(rPer ,-2) - pow(rApo ,-2)), 0.5);	
	}

};


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
			//0.5  * ((distFunc(E+0.001, J) - distFunc(E-0.001,J))/0.001)* om1(i,j);	
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
			dFdJ(i,j) = ///om2(i,j) * (distFunc(E+0.00001, J) - distFunc(E-0.00001,J)) / (2*0.00001) 
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


#endif