#include "Define_All_Includes_Header.h"
#include "NumericalFlux_Header.h"


void NumericalFluxInALLNodes(int const K, int const num_flag, double const beta,
	double * nx1, double * ny1, double * nx2, double * ny2, double * nx3, double * ny3,
	double **u, double **v, double **rho, double **p, double **alpha_l, int **EToE, int **EToF,//inputs
	double **NumFlux_1, double **NumFlux_2, double **NumFlux_3, double **NumFlux_4)

{

	parallel_for(1, K + 1, [&](int j)
	{
		int jM, jP, iM, iP;
		double num_flux1, num_flux2, num_flux3, num_flux4;
		double n_x1 = nx1[j];
		double n_x2 = nx2[j];
		double n_x3 = nx3[j];
		double n_y1 = ny1[j];
		double n_y2 = ny2[j];
		double n_y3 = ny3[j];

		//face#1
		jM = j;
		iM = 1;
		jP = EToE[j][iM];
		iP = EToF[j][iM];

		numerical_flux(num_flag, n_x1, n_y1, beta,
			u[jM][iM], v[jM][iM], rho[jM][iM], p[jM][iM], alpha_l[jM][iM],
			u[jP][iP], v[jP][iP], rho[jP][iP], p[jP][iP], alpha_l[jP][iP],
			num_flux1, num_flux2, num_flux3, num_flux4);

		NumFlux_1[j][iM] = num_flux1;
		NumFlux_2[j][iM] = num_flux2;
		NumFlux_3[j][iM] = num_flux3;
		NumFlux_4[j][iM] = num_flux4;

		//face#2
		jM = j;
		iM = 2;
		jP = EToE[j][iM];
		iP = EToF[j][iM];

		numerical_flux(num_flag, n_x2, n_y2, beta,
			u[jM][iM], v[jM][iM], rho[jM][iM], p[jM][iM], alpha_l[jM][iM],
			u[jP][iP], v[jP][iP], rho[jP][iP], p[jP][iP], alpha_l[jP][iP],
			num_flux1, num_flux2, num_flux3, num_flux4);

		NumFlux_1[j][iM] = num_flux1;
		NumFlux_2[j][iM] = num_flux2;
		NumFlux_3[j][iM] = num_flux3;
		NumFlux_4[j][iM] = num_flux4;

		//face#3
		jM = j;
		iM = 3;
		jP = EToE[j][iM];
		iP = EToF[j][iM];

		numerical_flux(num_flag, n_x3, n_y3, beta,
			u[jM][iM], v[jM][iM], rho[jM][iM], p[jM][iM], alpha_l[jM][iM],
			u[jP][iP], v[jP][iP], rho[jP][iP], p[jP][iP], alpha_l[jP][iP],
			num_flux1, num_flux2, num_flux3, num_flux4);

		NumFlux_1[j][iM] = num_flux1;
		NumFlux_2[j][iM] = num_flux2;
		NumFlux_3[j][iM] = num_flux3;
		NumFlux_4[j][iM] = num_flux4;
	});
}

