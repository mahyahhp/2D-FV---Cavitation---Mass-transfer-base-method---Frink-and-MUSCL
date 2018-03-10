
#ifndef _InitialCondition_
#define _InitialCondition_
extern void InitialCondition(int const K, //inputs
	double &C_prod, double &C_dest, double &tinf, double &p_v, double &rho_v, double &beta, double &sigma, double &AOA,
	double *rho, double *u, double *v, double *p, double *alpha_l, double *m_dot_minus, double *m_dot_plus,
	double *Qo_1, double *Qo_2, double *Qo_3, double *Qo_4,
	double *F1, double *F2, double *F3, double *F4,
	double *G1, double *G2, double *G3, double *G4,
	double *Source1, double *Source2, double *Source3, double *Source4,
	double & uinf, double &  vinf, double &  rho_m_inf, double &  alpha_l_inf, double &  pinf);
#endif