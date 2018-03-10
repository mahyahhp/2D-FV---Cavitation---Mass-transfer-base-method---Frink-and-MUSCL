#include "Define_All_Includes_Header.h"
#include "PrimitiveVariables_Header.h"
#include "FluxVectors_Header.h"
#include "SourceVector_Header.h"
#include "SolutionVector_Header.h"

void UpdateVariablesForNextIteration(int const K,
	double *Qn_1, double *Qn_2, double *Qn_3, double *Qn_4,
	double const C_prod, double const C_dest, double const tinf, double const p_v, double const rho_v, //inputs
	double *rho, double *u, double *v, double *p, double *alpha_l, double *m_dot_minus, double *m_dot_plus,
	double *Qo_1, double *Qo_2, double *Qo_3, double *Qo_4,
	double *F1, double *F2, double *F3, double *F4,
	double *G1, double *G2, double *G3, double *G4,
	double *Source1, double *Source2, double *Source3, double *Source4
	)
{
	parallel_for(1, K + 1, [&](int j)
	{
		PrimitiveVariables(Qn_1[j], Qn_2[j], Qn_3[j], Qn_4[j],
			C_prod, C_dest, tinf, p_v, rho_v, //inputs
			p[j], u[j], v[j], rho[j], alpha_l[j], m_dot_minus[j], m_dot_plus[j]);

		FluxVectors(p[j], u[j], v[j], rho[j], alpha_l[j],  //inputs
			F1[j], F2[j], F3[j], F4[j], G1[j], G2[j], G3[j], G4[j]);

		SourceVector(m_dot_minus[j], m_dot_plus[j], rho_v, //inputs
			Source1[j], Source2[j], Source3[j], Source4[j]);

		SolutionVector(p[j], u[j], v[j], alpha_l[j],  //inputs
			Qo_1[j], Qo_2[j], Qo_3[j], Qo_4[j]);
	});


}

