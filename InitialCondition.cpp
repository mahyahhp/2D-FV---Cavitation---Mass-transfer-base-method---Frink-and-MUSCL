#include "Define_All_Includes_Header.h"
#include "PrimitiveVariables_Header.h"
#include "FluxVectors_Header.h"
#include "SourceVector_Header.h"
#include "SolutionVector_Header.h"

void InitialCondition(int const K, //inputs
	double &C_prod, double &C_dest, double &tinf, double &p_v, double &rho_v, double &beta, double &sigma, double &AOA,
	double *rho, double *u, double *v, double *p, double *alpha_l, double *m_dot_minus, double *m_dot_plus,
	double *Qo_1, double *Qo_2, double *Qo_3, double *Qo_4,
	double *F1, double *F2, double *F3, double *F4,
	double *G1, double *G2, double *G3, double *G4,
	double *Source1, double *Source2, double *Source3, double *Source4,
	double & uinf, double &  vinf, double &  rho_m_inf, double &  alpha_l_inf, double &  pinf)
{
	int i, j, ini_set;;
	double ini_coef, rho_l_sat, p_sat, rho_v_sat, drho_l_sat_dT, dp_sat_dT, drho_l_sat_dp, dp_sat_dp,
		rho_g, rho_l, ev, el, yv, yl, C3l, C3v, temp_RHS, pi=acos(-1);

	cout << "Cavitation number= ";
	cin >> sigma;
	//sigma = 0.9;
	cout << endl;

	cout << "AOA= ";
	cin >> AOA;
	//AOA = 5;
	cout << endl;

	rho_l = 0.01;//non-dimensional
	double rho_v_over_rho_l = 1;//non-dimensional
	rho_v = rho_v_over_rho_l*rho_l;//non-dimensional
	double p_inf = 1;//non-dimensional
	p_v = p_inf - sigma / 2; //non-dimensional
	double 	U_inf = 1; //m/sec
	double L_scale = 1; //airfoil chord in meter
	tinf = L_scale / U_inf;

	////Merkle
	cout << "C_dest(1)= ";
	cin >> C_dest;
	//C_dest = 1;
	cout << endl;
	cout << "C_prod(80)= ";
	cin >> C_prod;
	//C_prod = 80;
	cout << endl;
	//C_prod = 0.9;

	//////Kunz et al.
	//C_dest = 10e5;
	//C_prod = 200;
	//////C_dest = 0;
	//////C_prod = 0;

	beta = 1.4;

	for (j = 1; j <= K; j++)
	{

				p[j] = p_inf;
				u[j] = 1 * cos(AOA*pi / 180);
				v[j] = 1 * sin(AOA*pi / 180);
				alpha_l[j] = 1;
				rho[j] = alpha_l[j] + (1 - alpha_l[j]) * rho_v;

				SolutionVector(p[j], u[j], v[j], alpha_l[j],  //inputs
					Qo_1[j], Qo_2[j], Qo_3[j], Qo_4[j]);


				FluxVectors(p[j], u[j], v[j], rho[j], alpha_l[j],  //inputs
					F1[j], F2[j], F3[j], F4[j], G1[j], G2[j], G3[j], G4[j]);

				PrimitiveVariables(Qo_1[j], Qo_2[j], Qo_3[j], Qo_4[j],
					C_prod, C_dest, tinf, p_v, rho_v, //inputs
					p[j], u[j], v[j], rho[j], alpha_l[j], m_dot_minus[j], m_dot_plus[j]);

				SourceVector(m_dot_minus[j], m_dot_plus[j], rho_v, //inputs
					Source1[j], Source2[j], Source3[j], Source4[j]);


				uinf = u[j];
				vinf = v[j];
				rho_m_inf = rho[j];
				alpha_l_inf = alpha_l[j];
				pinf = p[j];
	}


}
