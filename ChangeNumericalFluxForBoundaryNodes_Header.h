#ifndef _ChangeNumericalFluxForBoundaryNodes_
#define _ChangeNumericalFluxForBoundaryNodes_
extern void ChangeNumericalFluxForBoundaryNodes(int const K, int const num_flag, double const beta,
	double * nx1, double * ny1, double * nx2, double * ny2, double * nx3, double * ny3,
	double **u, double **v, double **rho, double **p, double **alpha_l,
	int const NBCs, int *BC_NumElem, int **BC_Elem, int **BC_Edge,
	double const uinf, double const vinf, double const rho_m_inf, double const pinf, double const alpha_l_inf, //inputs
	double **NumFlux_1, double **NumFlux_2, double **NumFlux_3, double **NumFlux_4);

extern void BoundaryValues(double const nx, double const ny, int const BC_Counter,
	double const uinf, double const vinf, double const pinf, double const rho_m_inf, double const alpha_l_inf,
	double const u_minus, double const v_minus, double const p_minus, double const rho_minus, double const alpha_l_minus,
	double &uBC, double &vBC, double &pBC, double &rhoBC, double &alpha_lBC);
#endif