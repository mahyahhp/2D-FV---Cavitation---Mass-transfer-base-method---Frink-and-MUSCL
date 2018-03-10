#include "Define_All_Includes_Header.h"
#include "NumericalFlux_Header.h"

void PrimitiveVariables(double const Q1, double Q2, double const Q3, double const Q4, 
	double const C_prod, double const C_dest, double const tinf, double const p_v, double const rho_v,
	double &p, double &u, double &v, double &rho, double &alpha_l, double &m_dot_minus, double &m_dot_plus)
{
	int MaxSubiteration = 10000, subiteration;
	double T_old, T_new, p_new, omega, subiteration_error, 
		rho_l_sat, p_sat, rho_v_sat, drho_l_sat_dT, dp_sat_dT,
		drho_l_sat_dp, dp_sat_dp, C3l, C3v, temp_RHS;

	p = Q1;
	u = Q2;
	//if (u < 0)
	//{
	//	u = 0.001;
	//	Q2 = u;
	//}
	v = Q3;
	alpha_l = Q4;

	rho = alpha_l + (1 - alpha_l)*rho_v;

	//Merkle
	m_dot_minus = (C_dest / tinf)*alpha_l*minval(0, p - p_v) / rho_v;
	m_dot_plus = (C_prod / tinf)*(1 - alpha_l)*maxval(0, p - p_v);

	////Kunz et al.
	//m_dot_minus = C_dest*rho_v*alpha_l*minval(0,p - p_v);
	//m_dot_plus = C_prod*rho_v*(1 - alpha_l)*alpha_l*alpha_l;
}
