#ifndef ORBITSECTIONS
#define ORBITSECTIONS

#include <vector>
#include <iostream>
#include <valarray>

#include "../Bodies/Bodies.h"
#include "../Perturbation_Grids/BarGrid.h"


class OrbitSections
{
public:
	OrbitSections(const int nPtle = 1, const double upperX = 3) : m_nParticles{nPtle}, m_nPoints{250}, m_oldBarAngle{0},
	v_xValues(m_nParticles * m_nPoints), v_pxValues(m_nParticles * m_nPoints), v_indicies(m_nParticles),
	m_ptle(nPtle, upperX),
	m_oldXY(2*m_nParticles), m_oldVXVY(2*m_nParticles), m_oldAXAY(2*m_nParticles)
	{
		for (int i = 0; i < m_nParticles; ++i) {v_indicies[i] = 0;}
	}

	OrbitSections(const BarGrid & barGrid, const int nPtle = 20) : m_nParticles{nPtle}, m_nPoints{250}, m_oldBarAngle{0},
	v_xValues(m_nParticles * m_nPoints), v_pxValues(m_nParticles * m_nPoints), v_indicies(m_nParticles),
	m_ptle("Bodies/particleSamplesSections.out", nPtle, 1),
	m_oldXY(2*m_nParticles), m_oldVXVY(2*m_nParticles), m_oldAXAY(2*m_nParticles)
	{
		std::cout << m_ptle.xy[0] << " " << m_ptle.xy[1] << '\n';
	}



	~OrbitSections() {}
	
	Bodies m_ptle;  
	
	void driftStep(const double timeStep);
	void kickStep(const std::valarray<double> & accels, const double timeStep);
	void kickStep(const BarGrid & barGrid, const std::valarray<double> & accels, const double timeStep);
	void checkSections();
	void checkSections(const BarGrid & barGrid); 
	void angularMomentumSections(double barAngle);

	bool continueSections(const BarGrid & barGrid);
	bool continueSections();

	int minIndex() {
		int minValue{m_nPoints};
		for (auto it : v_indicies) {
			if (it < minValue) {minValue = it;} 
		}
	return minValue;
	}

	void printIndex() const {for (auto it : v_indicies) {std::cout << it << " ";} std::cout << '\n';}

	void outputSections(const std::string & filename); 


private:
	const int m_nParticles, m_nPoints;
	double m_oldBarAngle; 
	std::vector<double> v_xValues, v_pxValues; std::vector<int> v_indicies; 
	std::valarray<double> m_oldXY, m_oldVXVY, m_oldAXAY; 

	double radius(const double x, const double y) {return sqrt(x*x +y*y);}
	//bool crossAxis(double previousAng, double newAngle, double newBarAngle) const;
	bool crossAxis(int i, double newAngle, double patternSpeed) const;

};


void OrbitSections::driftStep(const double timeStep) {
	m_ptle.xy += m_ptle.vxvy * timeStep * 0.5; 
}

void OrbitSections::kickStep(const std::valarray<double> & accels, const double timeStep) {
	m_oldVXVY = m_ptle.vxvy;
	m_ptle.vxvy += 0.5 * (m_oldAXAY + accels) * timeStep; 
	m_oldAXAY = accels; 
	checkSections();
}

void OrbitSections::kickStep(const BarGrid & barGrid, const std::valarray<double> & accels, const double timeStep) {
	m_oldVXVY = m_ptle.vxvy;
	m_ptle.vxvy += (m_oldAXAY + accels) * timeStep; 
	m_oldAXAY = accels; 
	checkSections(barGrid);
}


void OrbitSections::checkSections() {
	for (int i =0; i <m_nParticles; ++i) {
		if ((m_ptle.xy[2*i +1] > 0 && m_oldXY[2*i +1] < 0)   && v_indicies[i] < m_nPoints) {
			v_xValues[i + v_indicies[i] * m_nParticles] = m_oldXY[2*i] + (m_ptle.xy[2*i] - m_oldXY[2*i]) / (1 - m_ptle.xy[2*i +1]/m_oldXY[2*i +1]);
			v_pxValues[i + v_indicies[i] * m_nParticles] = m_oldVXVY[2*i] + (m_ptle.vxvy[2*i] - m_oldVXVY[2*i]) / (1 - m_ptle.xy[2*i +1]/m_oldXY[2*i +1]);
			v_indicies[i] += 1; 
		}
	}
}


bool OrbitSections::crossAxis(int i, double newAngle, double patternSpeed) const {
	double oldYprime{m_oldXY[2*i+1] * cos(m_oldBarAngle) - m_oldXY[2*i] * sin(m_oldBarAngle)}, newYprime{m_ptle.xy[2*i+1] * cos(newAngle) - m_ptle.xy[2*i] * sin(newAngle)};
	double radius = sqrt(m_oldXY[2*i]*m_oldXY[2*i] + m_oldXY[2*i+1]*m_oldXY[2*i+1]);
	double vYprime = ((m_oldXY[2*i] * m_oldVXVY[2*i+1] - m_oldXY[2*i+1] * m_oldVXVY[2*i]))/radius - patternSpeed*radius;
	if ((oldYprime*newYprime <= 0) && vYprime < 0) {return true;} // Please change back! 
	return false;
}

void OrbitSections::checkSections(const BarGrid & barGrid) {
	double barAngle{barGrid.angle() -  2*M_PI * floor(barGrid.angle()/(2*M_PI))};

	for (int i = 0; i < m_nParticles; ++i) {
 		if (crossAxis(i, barAngle, barGrid.patternSpeed()) && (v_indicies[i] < m_nPoints)) 
		{
			auto xPrime = [] (double x, double y, double ang) {return x*cos(ang) + y * sin(ang);};
									
			double x0{xPrime(m_oldXY[2*i], m_oldXY[2*i+1], m_oldBarAngle)}, x1{xPrime(m_ptle.xy[2*i], m_ptle.xy[2*i+1], barAngle)};
			double px0{barGrid.corotatingPx(m_oldVXVY[2*i], m_oldXY[2*i+1])}, px1{barGrid.corotatingPx(m_ptle.vxvy[2*i], m_ptle.xy[2*i+1])};
		
			px0 = cos(m_oldBarAngle) * m_oldVXVY[2*i] + sin(m_oldBarAngle) * m_oldVXVY[2*i+1];
			px1 = cos(barAngle) * m_ptle.vxvy[2*i] + sin(barAngle) * m_ptle.vxvy[2*i+1];

			v_xValues[i + v_indicies[i] * m_nParticles] = 0.5*(x0 + x1); //+ (x1 -x0) * frac;
			v_pxValues[i + v_indicies[i] * m_nParticles] = 0.5*(px0 + px1);  //+ (px1 - px0) * frac;
			v_indicies[i] += 1; 
		}
	}
	m_oldXY = m_ptle.xy; m_oldVXVY = m_ptle.vxvy; m_oldBarAngle = barAngle;
}

void OrbitSections::angularMomentumSections(double barAngle) {
	// Can we save the previous A.M and angles as m_oldXY and m_oldVXVY??
	std::valarray<double> angles = m_ptle.angle();
	std::valarray<double> angMom = m_ptle.angularMomentum(); 
	barAngle =barAngle -  2*M_PI * floor(barAngle/(2*M_PI));


	auto pr = [] (double vx, double vy, double angle) {return vx*cos(angle) + vy*sin(angle);};

	for (int i = 0; i < m_nParticles; ++i) {
		double ang{atan2(m_ptle.xy[2*i+1], m_ptle.xy[2*i])};
		if ((pr(m_ptle.vxvy[2*i], m_ptle.vxvy[2*i+1], ang) < 0) && ((pr(m_oldVXVY[2*i], m_oldVXVY[2*i+1], atan2(m_oldXY[2*i+1], m_oldXY[2*i])) > 0))) {
			
			if (barAngle > M_PI) {barAngle -= 2*M_PI;}
			double holding{ang - (barAngle)};

			if (holding >  M_PI) {holding -= 2*M_PI;}
			if (holding < -M_PI) {holding += 2*M_PI;}

			//if(holding > 0) {holding -= M_PI;}
			//if(holding < 0) {holding += M_PI;}
			v_xValues[i + v_indicies[i] * m_nParticles] = holding;// - (barAngle- M_PI);
			if (i ==0) {std::cout <<holding <<  '\n';}
			v_pxValues[i + v_indicies[i] * m_nParticles] = angMom[i];
			v_indicies[i] += 1; 
		}
	}
	m_oldXY = m_ptle.xy; m_oldVXVY = m_ptle.vxvy;
}

bool OrbitSections::continueSections() {
		if (minIndex() == m_nPoints) {return false;}	
		return true;
	}

bool OrbitSections::continueSections(const BarGrid & barGrid) 
	{
		checkSections(barGrid);
		for (auto it = v_indicies.begin(); it != v_indicies.end(); ++it) {
			if (*it != m_nPoints) {return true;}
		}
		return false;
	}

void OrbitSections::outputSections(const std::string & filename) {
	std::ofstream out(filename);
	for (int i = 0; i < v_xValues.size(); i++){
		if ((i+1) % m_nParticles == 0) {out << v_xValues[i] << ',' << v_pxValues[i] << '\n';}
		else {out << v_xValues[i] << ',' << v_pxValues[i] << ',' ;}
	}
	std::cout << "Do we get here\n";
	out.close();
}
#endif