#ifndef _NumericalFluxInALLNodes_
#define _NumericalFluxInALLNodes_
extern void NumericalFluxInALLNodes(int const K, int const num_flag, double const beta,
	double * nx1, double * ny1, double * nx2, double * ny2, double * nx3, double * ny3,
	double **u, double **v, double **rho, double **p, double **alpha_l, int **EToE, int **EToF,//inputs
	double **NumFlux_1, double **NumFlux_2, double **NumFlux_3, double **NumFlux_4);
#endif