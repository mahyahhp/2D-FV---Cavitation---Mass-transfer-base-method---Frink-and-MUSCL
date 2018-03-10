#include "Define_All_Includes_Header.h"
#include "NumericalFlux_Header.h"


void BoundaryValues(double const nx, double const ny, int const BC_Counter,
	double const uinf, double const vinf, double const pinf, double const rho_m_inf, double const alpha_l_inf,
	double const u_minus, double const v_minus, double const p_minus, double const rho_minus, double const alpha_l_minus,
	double &uBC, double &vBC, double &pBC, double &rhoBC, double &alpha_lBC)
{


	if (BC_Counter == 1) //Subsonic Inflow Boundary Condition
	{

		uBC = uinf;
		vBC = vinf;
		rhoBC = rho_m_inf;
		alpha_lBC = alpha_l_inf;
		pBC = p_minus;
	}
	else if (BC_Counter == 2) //Subsonic Outflow Boundary Condition
	{
		uBC = u_minus;
		vBC = v_minus;
		rhoBC = rho_minus;
		alpha_lBC = alpha_l_minus;
		pBC = pinf;
	}
	else if (BC_Counter == 3) //Slip Boundary Condition
	{
		pBC = p_minus;
		rhoBC = rho_minus;
		alpha_lBC = alpha_l_minus;

		uBC = u_minus * (ny * ny - nx * nx) - 2 * v_minus * nx * ny;
		vBC = v_minus * (nx * nx - ny * ny) - 2 * u_minus * nx * ny; //derived by myself
	}

	else if (BC_Counter == 4) //Slip Boundary Condition
	{
		pBC = p_minus;
		rhoBC = rho_minus;
		alpha_lBC = alpha_l_minus;

		uBC = u_minus * (ny * ny - nx * nx) - 2 * v_minus * nx * ny;
		vBC = v_minus * (nx * nx - ny * ny) - 2 * u_minus * nx * ny; //derived by myself
	}

}



void ChangeNumericalFluxForBoundaryNodes(int const K, int const num_flag, double const beta, 
	double * nx1, double * ny1, double * nx2, double * ny2, double * nx3, double * ny3,
	double **u, double **v, double **rho, double **p, double **alpha_l,
	int const NBCs, int *BC_NumElem, int **BC_Elem, int **BC_Edge,
	double const uinf, double const vinf, double const rho_m_inf, double const pinf, double const alpha_l_inf, //inputs
	double **NumFlux_1, double **NumFlux_2, double **NumFlux_3, double **NumFlux_4)
{
	//double u_minus, v_minus, p_minus, rho_minus, E_minus, H_minus, SoundSpeed_minus;
	//double uBC, vBC, pBC, rhoBC, EBC, HBC, SoundSpeedBC;
	//int i, j, ii, BC_Counter, k;
	//double num_flux1, num_flux2, num_flux3, num_flux4;

	for (int BC_Counter = 1; BC_Counter <= NBCs; BC_Counter++)
	{
		parallel_for(1, BC_NumElem[BC_Counter] + 1, [&](int BC_Elem_Counter)
		{
				double u_minus, v_minus, p_minus, rho_minus, alpha_l_minus;
				double uBC, vBC, pBC, rhoBC, alpha_lBC;
				double num_flux1, num_flux2, num_flux3, num_flux4;
				double nx, ny;

				int j = BC_Elem[BC_Elem_Counter][BC_Counter];
				int i = BC_Edge[BC_Elem_Counter][BC_Counter];

				u_minus = u[j][i];
				v_minus = v[j][i];
				p_minus = p[j][i];
				rho_minus = rho[j][i];
				alpha_l_minus = alpha_l[j][i];


				if (i == 1)
				{
					nx = nx1[j];
					ny = ny1[j];
					BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
						, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
						uBC, vBC, pBC, rhoBC, alpha_lBC);
				}
				else if (i == 2)
				{
					nx = nx2[j];
					ny = ny2[j];
					BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
						, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
						uBC, vBC, pBC, rhoBC, alpha_lBC);
				}
				else
				{
					nx = nx3[j];
					ny = ny3[j];
					BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
						, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
						uBC, vBC, pBC, rhoBC, alpha_lBC);
				}


				numerical_flux(num_flag, nx, ny, beta,
					u_minus, v_minus, rho_minus, p_minus, alpha_l_minus,
					uBC, vBC, rhoBC, pBC, alpha_lBC,
					num_flux1, num_flux2, num_flux3, num_flux4);

				NumFlux_1[j][i] = num_flux1;	//[Nfp*Nfaces) x K
				NumFlux_2[j][i] = num_flux2;
				NumFlux_3[j][i] = num_flux3;
				NumFlux_4[j][i] = num_flux4;

			//system("pause");
		});
	}
}

