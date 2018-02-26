#ifndef _Model_h_
#define _Model_h_

#include <string>
#include <cmath>

class Model
{
	public:

	std::string name;
	std::string label;
	double ranking;

	virtual double si_p_p(double s) const = 0;
	virtual double si_p_ap(double s) const = 0;

	virtual double rho_p_p(double s) const = 0;
	virtual double rho_p_ap(double s) const = 0;
};

#endif
