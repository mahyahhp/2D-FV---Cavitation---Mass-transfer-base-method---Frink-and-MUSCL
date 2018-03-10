#ifndef _numerical_flux_
#define _numerical_flux_
extern void numerical_flux(int const num_flag, double const nx, double const ny, double const beta,
	double const uL, double const vL, double const rhoL, double const pL, double const alpha_lL,
	double const uR, double const vR, double const rhoR, double const pR, double const alpha_lR,
	double &num_flux_1, double &num_flux_2, double &num_flux_3, double &num_flux_4);
double maxval(double const a, double const b);
double minval(double const a, double const b);
#endif

