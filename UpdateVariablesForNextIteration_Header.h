#ifndef _UpdateVariablesForNextIteration_
#define _UpdateVariablesForNextIteration_
extern void UpdateVariablesForNextIteration(int const K,
	double *Qn_1, double *Qn_2, double *Qn_3, double *Qn_4,
	double const C_prod, double const C_dest, double const tinf, double const p_v, double const rho_v, //inputs
	double *rho, double *u, double *v, double *p, double *alpha_l, double *m_dot_minus, double *m_dot_plus,
	double *Qo_1, double *Qo_2, double *Qo_3, double *Qo_4,
	double *F1, double *F2, double *F3, double *F4,
	double *G1, double *G2, double *G3, double *G4,
	double *Source1, double *Source2, double *Source3, double *Source4
	);
#endif