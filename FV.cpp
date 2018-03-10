#include "Define_All_Includes_Header.h"
#include "Define_Arrays_Header.h"
#include "Read_Gambit_Data_Header.h"
#include "ElementToElement_and_ElementToFace_Connectivity_Header.h"
#include "Jacobian_Matrices_NormalVectors_Header.h"
#include "NumericalFlux_Header.h"
#include "InitialCondition_Header.h"
#include "PrimitiveVariables_Header.h"
#include "PrimitiveVariables_Header.h"
#include "FluxVectors_Header.h"
#include "SourceVector_Header.h" 
#include "SolutionVector_Header.h"
#include "NumericalFluxInALLNodes_Header.h"
#include "ChangeNumericalFluxForBoundaryNodes_Header.h"
#include"UpdateVariablesForNextIteration_Header.h"
#include "VariablesInVertices_Header.h"


void InversionOfPreconditioningMatrix(int const K, double *u, double *v, double *rho, double *alpha_l, double const beta, double const drho_l,
	double *P_Matrix_inv11, double *P_Matrix_inv12, double *P_Matrix_inv13, double *P_Matrix_inv14,
	double *P_Matrix_inv21, double *P_Matrix_inv22, double *P_Matrix_inv23, double *P_Matrix_inv24,
	double *P_Matrix_inv31, double *P_Matrix_inv32, double *P_Matrix_inv33, double *P_Matrix_inv34,
	double *P_Matrix_inv41, double *P_Matrix_inv42, double *P_Matrix_inv43, double *P_Matrix_inv44)
{


	parallel_for(1, K + 1, [&](int j)
	{
		P_Matrix_inv11[j] = rho[j] * beta*beta;
		P_Matrix_inv12[j] = 0;
		P_Matrix_inv13[j] = 0;
		P_Matrix_inv14[j] = 0;

		P_Matrix_inv21[j] = u[j] * drho_l * alpha_l[j] / rho[j];
		P_Matrix_inv22[j] = 1 / rho[j];
		P_Matrix_inv23[j] = 0;
		P_Matrix_inv24[j] = -u[j] * drho_l / rho[j];

		P_Matrix_inv31[j] = v[j] * drho_l*alpha_l[j] / rho[j];
		P_Matrix_inv32[j] = 0;
		P_Matrix_inv33[j] = 1 / rho[j];
		P_Matrix_inv34[j] = -v[j] * drho_l / rho[j];

		P_Matrix_inv41[j] = -alpha_l[j];
		P_Matrix_inv42[j] = 0;
		P_Matrix_inv43[j] = 0;
		P_Matrix_inv44[j] = 1;
	});






}


double del2Q(int const Nfaces, double const Q, double **Q_neighbor, int const j)
{
	//double d2Q = 0;
	//for (int i = 1; i <= Nfaces; i++)
	//{
	//	d2Q = d2Q + (Q_neighbor[j][i] - Q[j]);
	//}

	double d2Q = (Q_neighbor[j][1] - Q) + (Q_neighbor[j][2] - Q) + (Q_neighbor[j][3] - Q);
	return d2Q;
}


void DistanceElementCenterToFaceCenter(int const K, int ** EToV, double *Vx, double *Vy, double **a)
{
	for (int j = 1; j <= K; j++)
	{
		double x_ElementCenter = (Vx[EToV[j][1]] + Vx[EToV[j][2]] + Vx[EToV[j][3]]) / 3;
		double y_ElementCenter = (Vx[EToV[j][1]] + Vx[EToV[j][2]] + Vx[EToV[j][3]]) / 3;

		double x_Face_1_center = (Vx[EToV[j][1]] + Vx[EToV[j][2]]) / 2;
		double x_Face_2_center = (Vx[EToV[j][2]] + Vx[EToV[j][3]]) / 2;
		double x_Face_3_center = (Vx[EToV[j][1]] + Vx[EToV[j][3]]) / 2;

		double y_Face_1_center = (Vy[EToV[j][1]] + Vy[EToV[j][2]]) / 2;
		double y_Face_2_center = (Vy[EToV[j][2]] + Vy[EToV[j][3]]) / 2;
		double y_Face_3_center = (Vy[EToV[j][1]] + Vy[EToV[j][3]]) / 2;

		a[j][1] = sqrt((x_ElementCenter - x_Face_1_center)*(x_ElementCenter - x_Face_1_center) + (y_ElementCenter - y_Face_1_center)*(y_ElementCenter - y_Face_1_center));
		a[j][2] = sqrt((x_ElementCenter - x_Face_2_center)*(x_ElementCenter - x_Face_2_center) + (y_ElementCenter - y_Face_2_center)*(y_ElementCenter - y_Face_2_center));
		a[j][3] = sqrt((x_ElementCenter - x_Face_3_center)*(x_ElementCenter - x_Face_3_center) + (y_ElementCenter - y_Face_3_center)*(y_ElementCenter - y_Face_3_center));
	}
}



void ContourPlot(int const Nv, int const K, double *Vx, double *Vy, double *u, double *v, double *p, double *rho, double *alpha_l,
	double *m_dot_minus, double *m_dot_plus, double *Source1, double *Source4, int **EToV)
{
	ofstream Contour("Contour.plt");
	Contour << "variables = \"x\" , \"y\" , \"u\" , \"v\" ,\"p\"  ,\"rho\"  ,\"alpha_l\"  ,\"m_dot_minus\"  ,\"m_dot_plus\"  ,\"Source1\"  ,\"Source4\" " << endl;
	Contour << "zone T=\"Numerical\" N= " << Nv << " , E= " << K << " , ZONETYPE=FETriangle" << endl;
	Contour << "DATAPACKING=BLOCK" << endl;
	Contour << "VARLOCATION=([1-2]=nodal,[3-11]=CELLCENTERED,)" << endl;

	for (int i = 1; i <= Nv; i = i + 1)
	{
		Contour << Vx[i] << " " << endl;
	}

	for (int i = 1; i <= Nv; i = i + 1)
	{
		Contour << Vy[i] << " " << endl;
	}

	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << u[i] << " " << endl;
	}

	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << v[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << p[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << rho[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << alpha_l[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << m_dot_minus[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << m_dot_plus[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << Source1[i] << " " << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << Source4[i] << " " << endl;
	}

	for (int i = 1; i <= K; i = i + 1)
	{
		Contour << EToV[i][1] << "	" << EToV[i][2] << "	" << EToV[i][3] << endl;
	}

}



void MUSCL(double const eta, double const Kappa, double const eps, int const j, int const FaceId, int const NodeID
	, int **EToE, int **EToF, int **EToV, double **a, double *u, double **u_Neighbor, double *u_node, double &u_f)
{
	int NeighborElement = EToE[j][FaceId];
	int NeighborElementFace = EToF[j][FaceId];
	double aa = a[j][FaceId];
	double bb = a[NeighborElement][NeighborElementFace];
	double IncrementCoef = (2 * aa) / (aa + bb);

	double Increment_u = IncrementCoef*(u_Neighbor[j][FaceId] - u[j]);

	double Nabla_u = u[j] - u_node[EToV[j][NodeID]];

	double s_u = (2 * Increment_u*Nabla_u + eps) / (Increment_u*Increment_u + Nabla_u*Nabla_u + eps);

	u_f = u[j] + 0.25*eta*s_u*((1 - Kappa*s_u)*Nabla_u + (1 + Kappa*s_u)*Increment_u);
}





void PrimitivVariablesForFluxComputation_MUSCL(double const Kappa, double const eta, int const K, double const rho_v, int **EToE, int **EToF, int **EToV, double **a,
	double *p, double *u, double *v, double *alpha_l, double **p_Neighbor, double **u_Neighbor, double **v_Neighbor, double **alpha_l_Neighbor,
	double *u_node, double *v_node, double *p_node, double *rho_node, double *alpha_l_node,
	double **u_flux, double **v_flux, double **p_flux, double **rho_flux, double **alpha_l_flux)
{
	double const eps = 1e-10;

	parallel_for(1, K + 1, [&](int j)
		//for (int j = 1; j < K + 1; j++)
	{
		double p_f, u_f, v_f, alpha_l_f, rho_f;

		//for 1st face
		int FaceId = 1;
		int NodeID = 3;

		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, p, p_Neighbor, p_node, p_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, u, u_Neighbor, u_node, u_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, v, v_Neighbor, v_node, v_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, alpha_l, alpha_l_Neighbor, alpha_l_node, alpha_l_f);

		p_flux[j][FaceId] = p_f;
		u_flux[j][FaceId] = u_f;
		v_flux[j][FaceId] = v_f;
		alpha_l_flux[j][FaceId] = alpha_l_f;
		rho_flux[j][FaceId] = alpha_l_f + (1 - alpha_l_f)*rho_v;

		//for 2nd face
		FaceId = 2;
		NodeID = 1;

		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, p, p_Neighbor, p_node, p_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, u, u_Neighbor, u_node, u_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, v, v_Neighbor, v_node, v_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, alpha_l, alpha_l_Neighbor, alpha_l_node, alpha_l_f);

		p_flux[j][FaceId] = p_f;
		u_flux[j][FaceId] = u_f;
		v_flux[j][FaceId] = v_f;
		alpha_l_flux[j][FaceId] = alpha_l_f;
		rho_flux[j][FaceId] = alpha_l_f + (1 - alpha_l_f)*rho_v;

		//for 3th face
		FaceId = 3;
		NodeID = 2;

		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, p, p_Neighbor, p_node, p_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, u, u_Neighbor, u_node, u_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, v, v_Neighbor, v_node, v_f);
		MUSCL(eta, Kappa, eps, j, FaceId, NodeID, EToE, EToF, EToV, a, alpha_l, alpha_l_Neighbor, alpha_l_node, alpha_l_f);

		p_flux[j][FaceId] = p_f;
		u_flux[j][FaceId] = u_f;
		v_flux[j][FaceId] = v_f;
		alpha_l_flux[j][FaceId] = alpha_l_f;
		rho_flux[j][FaceId] = alpha_l_f + (1 - alpha_l_f)*rho_v;
	});
	//}
}





void PrimitivVariablesForFluxComputation_Frink(double const eta, int const K, int **EToV, double *u, double *v, double *p, double *rho, double *alpha_l,
	double *u_node, double *v_node, double *p_node, double *rho_node, double *alpha_l_node,
	double **u_flux, double **v_flux, double **p_flux, double **rho_flux, double **alpha_l_flux)
{
	parallel_for(1, K + 1, [&](int j)
	{
		//for 1st face
		u_flux[j][1] = u[j] + eta*(0.5*(u_node[EToV[j][1]] + u_node[EToV[j][2]]) - u_node[EToV[j][3]]) / 3;
		v_flux[j][1] = v[j] + eta*(0.5*(v_node[EToV[j][1]] + v_node[EToV[j][2]]) - v_node[EToV[j][3]]) / 3;
		p_flux[j][1] = p[j] + eta*(0.5*(p_node[EToV[j][1]] + p_node[EToV[j][2]]) - p_node[EToV[j][3]]) / 3;
		rho_flux[j][1] = rho[j] + eta*(0.5*(rho_node[EToV[j][1]] + rho_node[EToV[j][2]]) - rho_node[EToV[j][3]]) / 3;
		alpha_l_flux[j][1] = alpha_l[j] + eta*(0.5*(alpha_l_node[EToV[j][1]] + alpha_l_node[EToV[j][2]]) - alpha_l_node[EToV[j][3]]) / 3;

		//for 2nd face
		u_flux[j][2] = u[j] + eta*(0.5*(u_node[EToV[j][3]] + u_node[EToV[j][2]]) - u_node[EToV[j][1]]) / 3;
		v_flux[j][2] = v[j] + eta*(0.5*(v_node[EToV[j][3]] + v_node[EToV[j][2]]) - v_node[EToV[j][1]]) / 3;
		p_flux[j][2] = p[j] + eta*(0.5*(p_node[EToV[j][3]] + p_node[EToV[j][2]]) - p_node[EToV[j][1]]) / 3;
		rho_flux[j][2] = rho[j] + eta*(0.5*(rho_node[EToV[j][3]] + rho_node[EToV[j][2]]) - rho_node[EToV[j][1]]) / 3;
		alpha_l_flux[j][2] = alpha_l[j] + eta*(0.5*(alpha_l_node[EToV[j][3]] + alpha_l_node[EToV[j][2]]) - alpha_l_node[EToV[j][1]]) / 3;

		//for 3th face
		u_flux[j][3] = u[j] + eta*(0.5*(u_node[EToV[j][1]] + u_node[EToV[j][3]]) - u_node[EToV[j][2]]) / 3;
		v_flux[j][3] = v[j] + eta*(0.5*(v_node[EToV[j][1]] + v_node[EToV[j][3]]) - v_node[EToV[j][2]]) / 3;
		p_flux[j][3] = p[j] + eta*(0.5*(p_node[EToV[j][1]] + p_node[EToV[j][3]]) - p_node[EToV[j][2]]) / 3;
		rho_flux[j][3] = rho[j] + eta*(0.5*(rho_node[EToV[j][1]] + rho_node[EToV[j][3]]) - rho_node[EToV[j][2]]) / 3;
		alpha_l_flux[j][3] = alpha_l[j] + eta*(0.5*(alpha_l_node[EToV[j][1]] + alpha_l_node[EToV[j][3]]) - alpha_l_node[EToV[j][2]]) / 3;
	});
}


void ReadDataFunction(int const Nv, int const K, double *Vx, double *Vy, double *u, double *v, double *p, double *rho, double *alpha_l, double *m_dot_plus, double *m_dot_minus, double *Source1, double *Source4)
{
	ifstream ReadData;
	ReadData.open("input.plt");
	char A[116];
	ReadData >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A >> A;

	for (int i = 1; i <= Nv; i = i + 1)
	{
		ReadData >> Vx[i];
	}

	for (int i = 1; i <= Nv; i = i + 1)
	{
		ReadData >> Vy[i];
	}

	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> u[i];
	}

	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> v[i];
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> p[i];
		//cout << p[i] << endl;
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> rho[i];
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> alpha_l[i];
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> m_dot_minus[i];
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> m_dot_plus[i];
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> Source1[i];
	}
	for (int i = 1; i <= K; i = i + 1)
	{
		ReadData >> Source4[i];
	}
}



void NeighborValues(int const K, int **EToE, int const NBCs, int *BC_NumElem, int **BC_Type_Elem, double const rho_v,
	double *nx1, double *nx2, double *nx3, double *ny1, double *ny2, double *ny3,
	double const uinf, double const vinf, double const pinf, double const rho_m_inf, double const alpha_l_inf,
	int const GhostElementCounter, int *BC_Ghost_Edge, int *BC_Ghost_Elem, int *BC_Type_Ghost, int **EToE_Ghost,
	int **BC_Type_Edge, double *Qo_1, double *Qo_2, double *Qo_3, double *Qo_4,
	double **Q_neighbour1, double **Q_neighbour2, double **Q_neighbour3, double **Q_neighbour4)
{


	parallel_for(1, K + 1, [&](int j)
		//for (int j = 1; j < K + 1; j++)
	{
		//for (int j = 1; j <= K; j++)
		//{
		Q_neighbour1[j][1] = Qo_1[EToE[j][1]];
		Q_neighbour2[j][1] = Qo_2[EToE[j][1]];
		Q_neighbour3[j][1] = Qo_3[EToE[j][1]];
		Q_neighbour4[j][1] = Qo_4[EToE[j][1]];

		Q_neighbour1[j][2] = Qo_1[EToE[j][2]];
		Q_neighbour2[j][2] = Qo_2[EToE[j][2]];
		Q_neighbour3[j][2] = Qo_3[EToE[j][2]];
		Q_neighbour4[j][2] = Qo_4[EToE[j][2]];

		Q_neighbour1[j][3] = Qo_1[EToE[j][3]];
		Q_neighbour2[j][3] = Qo_2[EToE[j][3]];
		Q_neighbour3[j][3] = Qo_3[EToE[j][3]];
		Q_neighbour4[j][3] = Qo_4[EToE[j][3]];

		//cout <<"j = "<<j<< ", j(1) = " << EToE[j][1] << ", j(2) = " << EToE[j][2] << ", j(3) = " << EToE[j][3] << endl;

		//Qo_1_W_Ghost[j] = Qo_1[j];
		//Qo_2_W_Ghost[j] = Qo_2[j];
		//Qo_3_W_Ghost[j] = Qo_3[j];
		//Qo_4_W_Ghost[j] = Qo_4[j];
		//}
	});
	//}
	//revising for real boundary elements

	for (int BC_Counter = 1; BC_Counter <= NBCs; BC_Counter++)
	{
		parallel_for(1, BC_NumElem[BC_Counter] + 1, [&](int BC_Elem_Counter)
			//for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter] + 1; BC_Elem_Counter++)
		{
			double u_minus, v_minus, p_minus, rho_minus, alpha_l_minus;
			double uBC, vBC, pBC, rhoBC, alpha_lBC;

			int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter];
			int i = BC_Type_Edge[BC_Elem_Counter][BC_Counter];

			p_minus = Qo_1[j];
			u_minus = Qo_2[j];
			v_minus = Qo_3[j];
			alpha_l_minus = Qo_4[j];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;


			if (i == 1)
			{
				double nx = nx1[j];
				double ny = ny1[j];
				BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
			}
			else if (i == 2)
			{
				double nx = nx2[j];
				double ny = ny2[j];
				BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
			}
			else
			{
				double nx = nx3[j];
				double ny = ny3[j];
				BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
			}

			Q_neighbour1[j][i] = pBC;
			Q_neighbour2[j][i] = uBC;
			Q_neighbour3[j][i] = vBC;
			Q_neighbour4[j][i] = alpha_lBC;
		});
		//}
	}


	//Values for ghost elements
	for (int GEC = K + 1; GEC <= GhostElementCounter; GEC++)
	{
		double u_minus, v_minus, p_minus, rho_minus, alpha_l_minus;
		double uBC, vBC, pBC, rhoBC, alpha_lBC;

		int i = BC_Ghost_Edge[GEC];
		int j = BC_Ghost_Elem[GEC]; //elemane vagheii ke ba aan dar ertebat ast!
		int k = BC_Type_Ghost[GEC];

		if (i == 1)
		{
			double nx = nx1[j];
			double ny = ny1[j];

			Q_neighbour1[GEC][1] = Qo_1[j];
			Q_neighbour2[GEC][1] = Qo_2[j];
			Q_neighbour3[GEC][1] = Qo_3[j];
			Q_neighbour4[GEC][1] = Qo_4[j];

			p_minus = Qo_1[EToE_Ghost[GEC][2]];
			u_minus = Qo_2[EToE_Ghost[GEC][2]];
			v_minus = Qo_3[EToE_Ghost[GEC][2]];
			alpha_l_minus = Qo_4[EToE_Ghost[GEC][2]];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
			BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
				, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
				uBC, vBC, pBC, rhoBC, alpha_lBC);
			Q_neighbour1[GEC][2] = pBC;
			Q_neighbour2[GEC][2] = uBC;
			Q_neighbour3[GEC][2] = vBC;
			Q_neighbour4[GEC][2] = alpha_lBC;

			p_minus = Qo_1[EToE_Ghost[GEC][3]];
			u_minus = Qo_2[EToE_Ghost[GEC][3]];
			v_minus = Qo_3[EToE_Ghost[GEC][3]];
			alpha_l_minus = Qo_4[EToE_Ghost[GEC][3]];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
			BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
				, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
				uBC, vBC, pBC, rhoBC, alpha_lBC);
			Q_neighbour1[GEC][3] = pBC;
			Q_neighbour2[GEC][3] = uBC;
			Q_neighbour3[GEC][3] = vBC;
			Q_neighbour4[GEC][3] = alpha_lBC;

		}
		else if (i == 2)
		{
			double nx = nx2[j];
			double ny = ny2[j];

			Q_neighbour1[GEC][2] = Qo_1[j];
			Q_neighbour2[GEC][2] = Qo_2[j];
			Q_neighbour3[GEC][2] = Qo_3[j];
			Q_neighbour4[GEC][2] = Qo_4[j];

			p_minus = Qo_1[EToE_Ghost[GEC][1]];
			u_minus = Qo_2[EToE_Ghost[GEC][1]];
			v_minus = Qo_3[EToE_Ghost[GEC][1]];
			alpha_l_minus = Qo_4[EToE_Ghost[GEC][1]];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
			BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
				, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
				uBC, vBC, pBC, rhoBC, alpha_lBC);
			Q_neighbour1[GEC][1] = pBC;
			Q_neighbour2[GEC][1] = uBC;
			Q_neighbour3[GEC][1] = vBC;
			Q_neighbour4[GEC][1] = alpha_lBC;

			p_minus = Qo_1[EToE_Ghost[GEC][3]];
			u_minus = Qo_2[EToE_Ghost[GEC][3]];
			v_minus = Qo_3[EToE_Ghost[GEC][3]];
			alpha_l_minus = Qo_4[EToE_Ghost[GEC][3]];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
			BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
				, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
				uBC, vBC, pBC, rhoBC, alpha_lBC);
			Q_neighbour1[GEC][3] = pBC;
			Q_neighbour2[GEC][3] = uBC;
			Q_neighbour3[GEC][3] = vBC;
			Q_neighbour4[GEC][3] = alpha_lBC;
		}
		else
		{
			double nx = nx3[j];
			double ny = ny3[j];

			Q_neighbour1[GEC][3] = Qo_1[j];
			Q_neighbour2[GEC][3] = Qo_2[j];
			Q_neighbour3[GEC][3] = Qo_3[j];
			Q_neighbour4[GEC][3] = Qo_4[j];

			p_minus = Qo_1[EToE_Ghost[GEC][2]];
			u_minus = Qo_2[EToE_Ghost[GEC][2]];
			v_minus = Qo_3[EToE_Ghost[GEC][2]];
			alpha_l_minus = Qo_4[EToE_Ghost[GEC][2]];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
			BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
				, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
				uBC, vBC, pBC, rhoBC, alpha_lBC);
			Q_neighbour1[GEC][2] = pBC;
			Q_neighbour2[GEC][2] = uBC;
			Q_neighbour3[GEC][2] = vBC;
			Q_neighbour4[GEC][2] = alpha_lBC;

			p_minus = Qo_1[EToE_Ghost[GEC][1]];
			u_minus = Qo_2[EToE_Ghost[GEC][1]];
			v_minus = Qo_3[EToE_Ghost[GEC][1]];
			alpha_l_minus = Qo_4[EToE_Ghost[GEC][1]];
			rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
			BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
				, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
				uBC, vBC, pBC, rhoBC, alpha_lBC);
			Q_neighbour1[GEC][1] = pBC;
			Q_neighbour2[GEC][1] = uBC;
			Q_neighbour3[GEC][1] = vBC;
			Q_neighbour4[GEC][1] = alpha_lBC;
		}

	}


}




void main()
{
	//declare variables
	int num_flag, Nv, K, NBCs;
	double EndTime, dt, *Vx, *Vy;
	//=== input data ===//
	int const Nfaces = 3;
	double const NODETOL = 1e-10;
	double const pi = acos(-1.0);

	cout << "EndTime (100)= "; //final time
	cin >> EndTime;
	//EndTime = 100;
	cout << endl;

	cout << "dt (1e-4)= ";
	cin >> dt;
	//dt = 1e-4;
	cout << endl;

	cout << "numerical flax: 1= Roe, 2=Lax, 3=AUSM+, 4=HLL, 5=HLLC, 6=LDFSS, 7=Average " << endl;
	cin >> num_flag;
	//num_flag = 7;
	cout << endl;

	double k2p = 0.0;

	cout << "k2rho (1.0/10.0) = " << endl;
	double k2rho;
	cin >> k2rho;
	//k2rho = 0.1;

	cout << "k4 (1.0/32.0) = " << endl;
	double k4;
	cin >> k4;
	//k4 = 0.03125; 

	//double eta = 1; //eta=0 first order and eta=1 for high order methods
	//double const Kappa = 1.l/3.l; //-1: fully upwind scheme, 0: Fromm's scheme, 1: centeral differencing, and 1/3: third order upwind-biased scheme
	//string HighOrderMethod = "MUSCL"; //"MUSCL" or "Frink"

	double eta = 0; //eta=0 first order and eta=1 for high order methods 
	double const Kappa = 1.l / 3.l; //-1: fully upwind scheme, 0: Fromm's scheme, 1: centeral differencing, and 1/3: third order upwind-biased scheme
	string HighOrderMethod = "Frink"; //"MUSCL" or "Frink"
	 
	///****************///
	///*** Start up ***///
	///****************///

	//=== Read data ===//
	//Nv=number of nodes in the solution domain
	//K=number of elements in the solution domain
	//NBCs= number of the prescribed boundary conditions
	read_data_control_info(Nv, K, NBCs);
	//=== Determining of array dimensiones ===//
	//Vx=x position of nodes
	//Vx=y position of nodes
	Vx = d1array(Nv + 1), Vy = d1array(Nv + 1);

	int **EToV, **EToE, **EToF, **BC_Type_Edge, *BC_NumElem, **BC_Type_Elem, NIter, *Num_NodeCell, **NodeCell, **NodeNode;
	double time, *Area, **LengFace, *nx1, *nx2, *nx3, *ny1, *ny2, *ny3;
	double *u, *v, *p, *rho, *alpha_l, *m_dot_minus, *m_dot_plus, *Source1, *Source2, *Source3, *Source4,
		*Qn_1, *Qn_2, *Qn_3, *Qn_4, *Qo_1, *Qo_2, *Qo_3, *Qo_4, *F1, *F2, *F3, *F4, *G1, *G2, *G3, *G4;
	double uinf, vinf, rho_m_inf, pinf, alpha_l_inf;
	double **NumFlux_1, **NumFlux_2, **NumFlux_3, **NumFlux_4;
	double *R_1, *R_2, *R_3, *R_4, error;
	double *P_Matrix_inv11, *P_Matrix_inv12, *P_Matrix_inv13, *P_Matrix_inv14,
		*P_Matrix_inv21, *P_Matrix_inv22, *P_Matrix_inv23, *P_Matrix_inv24,
		*P_Matrix_inv31, *P_Matrix_inv32, *P_Matrix_inv33, *P_Matrix_inv34,
		*P_Matrix_inv41, *P_Matrix_inv42, *P_Matrix_inv43, *P_Matrix_inv44;
	double **Q_neighbour1, **Q_neighbour2, **Q_neighbour3, **Q_neighbour4;
	double *d1, *d2, *d3, *d4;
	double *D1, *D2, *D3, *D4;

	int *BC_Elem_ID, *BC_Edge_Counter, **BC_Edge, **BC_Edge_type;

	double *u_node, *v_node, *p_node, *rho_node, *alpha_l_node;
	double **u_flux, **v_flux, **p_flux, **rho_flux, **alpha_l_flux;

	double CavNum, AOA;

	double **a;

	//=== Allocating arrays ===//

	Area = d1array(K + 1);
	LengFace = d2array(K + 1, Nfaces + 1);
	//EToE[j][i]= Moshakhas mikonad ke elemane "j"-om, az tarighe face "i"-om khod ba elemane shomare chand dar ertebat ast
	//EToE[j][i]= Moshakhas mikonad ke elemane "j"-om, az tarighe face "i"-om khod ba face shomare chand az elemane hamsaye dar ertebat ast
	//EToV[j][i]= Moshakhas mikonad ke elemane "j"-om, gerehe "i"-om aan gereh shomare chand mibashand
	EToV = i2array(K + 1, Nfaces + 1); EToE = i2array(K + 1, Nfaces + 1); EToF = i2array(K + 1, Nfaces + 1);
	nx1 = d1array(K + 1); ny1 = d1array(K + 1); nx2 = d1array(K + 1); ny2 = d1array(K + 1); nx3 = d1array(K + 1); ny3 = d1array(K + 1);

	u = d1array(K + 1); v = d1array(K + 1); p = d1array(K + 1); rho = d1array(K + 1); alpha_l = d1array(K + 1);
	m_dot_minus = d1array(K + 1); m_dot_plus = d1array(K + 1); Source1 = d1array(K + 1); Source2 = d1array(K + 1); Source3 = d1array(K + 1); Source4 = d1array(K + 1);
	Qn_1 = d1array(K + 1); Qn_2 = d1array(K + 1); Qn_3 = d1array(K + 1); Qn_4 = d1array(K + 1);
	Qo_1 = d1array(K + 1); Qo_2 = d1array(K + 1); Qo_3 = d1array(K + 1); Qo_4 = d1array(K + 1);
	F1 = d1array(K + 1); F2 = d1array(K + 1); F3 = d1array(K + 1); F4 = d1array(K + 1);
	G1 = d1array(K + 1); G2 = d1array(K + 1); G3 = d1array(K + 1); G4 = d1array(K + 1);

	//Num_NodeCell[i]: moshakhas mikonad ke har gerehe "i"-ome tolid shodah dar Gambit ba chand eleman dargir ast
	//NodeCell[i][j]: moshakhas mikonad ke baraye gerehe "i"-ome tolid shodeh dar Gambit, "j"-omin elemani ke ba an dar artebat ast elemane chandom ast
	//NodeNode[i][j]: moshakhas mikonad ke baraye gerehe "i"-ome tolid shodeh dar Gambit, noghte chandom az elemane "j"-om ba an dar artebat ast 
	Num_NodeCell = i1array(Nv + 1);
	NodeCell = i2array(Nv + 1, 21); //A node will be shareed with maximum 21 elements.
	NodeNode = i2array(Nv + 1, 21);

	//BC_NumElem[i]=tedade kole edge haye marzi az dasteye i om
	//BC_Type_Edge[j][i]=moshakhas mikonad ke face/edge shomare chand az "j"-omin eleman az "i"-omin dasteye sharte marzi rooye marz gharar darad
	//BC_Type_Elem[j][i]=moshakhas mikonad ke elemane shomare chand az "j"-omin eleman az "i"-omin dasteye sharte marzi rooye marz gharar darad
	BC_NumElem = i1array(NBCs + 1); BC_Type_Elem = i2array(Nv + 1, NBCs + 1); BC_Type_Edge = i2array(Nv + 1, NBCs + 1);
	NumFlux_1 = d2array(K + 1, Nfaces + 1); NumFlux_2 = d2array(K + 1, Nfaces + 1); NumFlux_3 = d2array(K + 1, Nfaces + 1); NumFlux_4 = d2array(K + 1, Nfaces + 1);
	R_1 = d1array(K + 1); R_2 = d1array(K + 1); R_3 = d1array(K + 1); R_4 = d1array(K + 1);

	P_Matrix_inv11 = d1array(K + 1); P_Matrix_inv12 = d1array(K + 1); P_Matrix_inv13 = d1array(K + 1); P_Matrix_inv14 = d1array(K + 1);
	P_Matrix_inv21 = d1array(K + 1); P_Matrix_inv22 = d1array(K + 1); P_Matrix_inv23 = d1array(K + 1); P_Matrix_inv24 = d1array(K + 1);
	P_Matrix_inv31 = d1array(K + 1); P_Matrix_inv32 = d1array(K + 1); P_Matrix_inv33 = d1array(K + 1); P_Matrix_inv34 = d1array(K + 1);
	P_Matrix_inv41 = d1array(K + 1); P_Matrix_inv42 = d1array(K + 1); P_Matrix_inv43 = d1array(K + 1); P_Matrix_inv44 = d1array(K + 1);

	d1 = d1array(K + 1); d2 = d1array(K + 1); d3 = d1array(K + 1); d4 = d1array(K + 1);
	D1 = d1array(K + 1); D2 = d1array(K + 1); D3 = d1array(K + 1); D4 = d1array(K + 1);

	BC_Elem_ID = i1array(K + 1); BC_Edge_Counter = i1array(K + 1); BC_Edge = i2array(K + 1, 3 + 1); BC_Edge_type = i2array(K + 1, 10 + 1);

	u_node = d1array(Nv + 1); v_node = d1array(Nv + 1); p_node = d1array(Nv + 1); rho_node = d1array(Nv + 1); alpha_l_node = d1array(Nv + 1);
	u_flux = d2array(K + 1, Nfaces + 1); v_flux = d2array(K + 1, Nfaces + 1); p_flux = d2array(K + 1, Nfaces + 1); rho_flux = d2array(K + 1, Nfaces + 1); alpha_l_flux = d2array(K + 1, Nfaces + 1);

	a = d2array(K + 1, Nfaces + 1); //a[j][i]=moshakhas mikonad ke fasele az markaze elemane "j"-om ta markaze face "i"-om cheghadr ast

	//==========================================


	read_main_data(Vx, Vy, EToV, Area, LengFace, BC_NumElem, BC_Type_Elem, BC_Type_Edge, Num_NodeCell, NodeCell, NodeNode,
		BC_Elem_ID, BC_Edge_Counter, BC_Edge, BC_Edge_type);
	ElementToElement_and_ElementToFace_Connectivity(K, Nv, Nfaces, EToV, EToE, EToF);
	Jacobian_Matrices_NormalVectors(K, EToV, Vx, Vy, LengFace, nx1, ny1, nx2, ny2, nx3, ny3);


	DistanceElementCenterToFaceCenter(K, EToV, Vx, Vy, a);


	int Total_BC_Elem;
	if (NBCs == 0) Total_BC_Elem = 0;
	if (NBCs == 1) Total_BC_Elem = BC_NumElem[1];
	if (NBCs == 2) Total_BC_Elem = BC_NumElem[1] + BC_NumElem[2];
	if (NBCs == 3) Total_BC_Elem = BC_NumElem[1] + BC_NumElem[2] + BC_NumElem[3];
	if (NBCs == 4) Total_BC_Elem = BC_NumElem[1] + BC_NumElem[2] + BC_NumElem[3] + BC_NumElem[4];


	//Finding ghost elements neighbors
	int **EToE_Ghost, *BC_Ghost_Elem, *BC_Type_Ghost, *BC_Ghost_Edge;
	EToE_Ghost = i2array(K + Total_BC_Elem + 1, 3 + 1);
	BC_Ghost_Elem = i1array(K + Total_BC_Elem + 1); BC_Type_Ghost = i1array(K + Total_BC_Elem + 1); BC_Ghost_Edge = i1array(K + Total_BC_Elem + 1);
	Q_neighbour1 = d2array(K + Total_BC_Elem + 1, Nfaces + 1); Q_neighbour2 = d2array(K + Total_BC_Elem + 1, Nfaces + 1); Q_neighbour3 = d2array(K + Total_BC_Elem + 1, Nfaces + 1); Q_neighbour4 = d2array(K + Total_BC_Elem + 1, Nfaces + 1);

	for (int j = 1; j <= K; j++)
	{
		for (int i = 1; i <= 3; i++)
		{
			EToE_Ghost[j][i] = EToE[j][i];
		}
	}
	//revising element to element conectivity for both the boundary elements and ghost elements
	int GhostElementCounter = K;
	for (int BC_Counter = 1; BC_Counter <= NBCs; BC_Counter++)
	{
		for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter] + 1; BC_Elem_Counter++)
		{
			GhostElementCounter = GhostElementCounter + 1;
			int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter];

			BC_Ghost_Elem[GhostElementCounter] = j;
			BC_Type_Ghost[GhostElementCounter] = BC_Counter;

			if (BC_Type_Edge[BC_Elem_Counter][BC_Counter] == 1)
			{
				EToE_Ghost[j][1] = GhostElementCounter;
				BC_Ghost_Edge[GhostElementCounter] = 1;
				EToE_Ghost[GhostElementCounter][1] = j;
				EToE_Ghost[GhostElementCounter][2] = EToE[j][2];
				EToE_Ghost[GhostElementCounter][3] = EToE[j][3];
			}
			else if (BC_Type_Edge[BC_Elem_Counter][BC_Counter] == 2)
			{
				EToE_Ghost[j][2] = GhostElementCounter;
				BC_Ghost_Edge[GhostElementCounter] = 2;
				EToE_Ghost[GhostElementCounter][1] = EToE[j][1];
				EToE_Ghost[GhostElementCounter][2] = j;
				EToE_Ghost[GhostElementCounter][3] = EToE[j][3];
			}
			else if (BC_Type_Edge[BC_Elem_Counter][BC_Counter] == 3)
			{
				EToE_Ghost[j][3] = GhostElementCounter;
				BC_Ghost_Edge[GhostElementCounter] = 3;
				EToE_Ghost[GhostElementCounter][1] = EToE[j][1];
				EToE_Ghost[GhostElementCounter][2] = EToE[j][2];
				EToE_Ghost[GhostElementCounter][3] = j;
			}
		}
	}



	//for (int j = 1; j <= GhostElementCounter; j++)
	//{
	//	//cout << "j = " << j << ", i = 1 -> " << EToE_Ghost[j][1] << ", i = 2 -> " << EToE_Ghost[j][2] << ", i = 3 -> " << EToE_Ghost[j][3] << endl;
	//	cout << "j = " << j << ", BC_Ghost_Elem = " << BC_Ghost_Elem[j] << ", BC_Type_Ghost = " << BC_Type_Ghost[j] << ", BC_Ghost_Edge = " << BC_Ghost_Edge[j] << endl;
	//}



	//=== initial condition ===//
	double C_prod, C_dest, tinf, p_v, rho_v, beta;

	parallel_for(1, K + 1, [&](int j)
	{
		Source1[j] = 0;
		Source2[j] = 0;
		Source3[j] = 0;
		Source4[j] = 0;
	});

	InitialCondition(K, //inputs
		C_prod, C_dest, tinf, p_v, rho_v, beta, CavNum, AOA, rho, u, v, p, alpha_l, m_dot_minus, m_dot_plus,
		Qo_1, Qo_2, Qo_3, Qo_4, F1, F2, F3, F4, G1, G2, G3, G4, Source1, Source2, Source3, Source4,
		uinf, vinf, rho_m_inf, alpha_l_inf, pinf);


	//Read data as initial condition
	cout << "Do you want to use the previous results (yes=\"y\" and no=\"n\" )" << endl;
	char A[1];
	cin >> A[1];
	if (A[1] == 'y')
	{
		ReadDataFunction(Nv, K, Vx, Vy, u, v, p, rho, alpha_l, m_dot_plus, m_dot_minus, Source1, Source4);
	}



	parallel_for(1, K + 1, [&](int j)
	{
		D1[j] = 0;
		D2[j] = 0;
		D3[j] = 0;
		D4[j] = 0;
	});

	double drho_l = 1 - rho_v;

	InversionOfPreconditioningMatrix(K, u, v, rho, alpha_l, beta, drho_l,
		P_Matrix_inv11, P_Matrix_inv12, P_Matrix_inv13, P_Matrix_inv14,
		P_Matrix_inv21, P_Matrix_inv22, P_Matrix_inv23, P_Matrix_inv24,
		P_Matrix_inv31, P_Matrix_inv32, P_Matrix_inv33, P_Matrix_inv34,
		P_Matrix_inv41, P_Matrix_inv42, P_Matrix_inv43, P_Matrix_inv44);



	VariablesInVertices(Nv, Num_NodeCell, NodeNode, NodeCell, Area, u, v, p, rho, alpha_l,
		u_node, v_node, p_node, rho_node, alpha_l_node);




	NeighborValues(K, EToE, NBCs, BC_NumElem, BC_Type_Elem, rho_v, nx1, nx2, nx3, ny1, ny2, ny3,
		uinf, vinf, pinf, rho_m_inf, alpha_l_inf, GhostElementCounter,
		BC_Ghost_Edge, BC_Ghost_Elem, BC_Type_Ghost, EToE_Ghost,
		BC_Type_Edge, Qo_1, Qo_2, Qo_3, Qo_4,
		Q_neighbour1, Q_neighbour2, Q_neighbour3, Q_neighbour4);

	if (HighOrderMethod == "MUSCL")
	{
		PrimitivVariablesForFluxComputation_MUSCL(Kappa, eta, K, rho_v, EToE, EToF, EToV, a,
			Qo_1, Qo_2, Qo_3, Qo_4, Q_neighbour1, Q_neighbour2, Q_neighbour3, Q_neighbour4,
			u_node, v_node, p_node, rho_node, alpha_l_node,
			u_flux, v_flux, p_flux, rho_flux, alpha_l_flux);
	}
	else //Frink
	{
		PrimitivVariablesForFluxComputation_Frink(eta, K, EToV, u, v, p, rho, alpha_l,
			u_node, v_node, p_node, rho_node, alpha_l_node,
			u_flux, v_flux, p_flux, rho_flux, alpha_l_flux);
	}



	//set dt and requared iterations
	NIter = floor(EndTime / dt) + 1;
	//plot error
	ofstream Error_Plot("Error.plt");
	Error_Plot << "variables = \"Iteration\" , \"Error\" " << endl;

	//**** solution starting ****//
	//NIter = 1;

	clock_t t;
	t = clock();
	time = 0;
	int iter = 0;

	ofstream Information("Information.txt");
	Information << " *** Grid informations ***" << endl;
	Information << " Total Elements (K) = " << K << endl;
	Information << " Total global nodes (Nv) = " << Nv << endl;
	Information << " Number of specified boundary conditions = " << NBCs << endl;
	Information << " Total boundary elements = " << Total_BC_Elem << endl;
	Information << "=======================" << endl;
	Information << " *** All information about computational cost ***" << endl;
	Information << " End time = " << EndTime << endl;
	Information << " dt = " << dt << endl;
	Information << " Number of iterations = " << iter << endl;
	Information << "=======================" << endl;
	Information << " *** Informations about flow conditions ***" << endl;
	Information << " beta = " << beta << endl;
	Information << " rho_v = " << rho_v << endl;
	Information << " Cavitation number = " << 2 * (pinf - p_v) << endl;
	Information << " C_dest = " << C_dest << endl;
	Information << " C_prod = " << C_prod << endl;
	Information << " Angle of attack = " << AOA << endl;
	Information << "=======================" << endl;
	Information << " *** Informations about solver ***" << endl;
	Information << " Numerical Flux ID (1= Roe, 2=Lax, 3=AUSM+, 4=HLL, 5=HLLC, 6=LDFSS, 7=Average) = " << num_flag << endl;
	Information << " k2p = " << k2p << endl;
	Information << " k2rho = " << k2rho << endl;
	Information << " k4 = " << k4 << endl;
	Information << " High order method = " << HighOrderMethod << endl;
	Information << " eta (0=1st order and 1=2nd order) = " << eta << endl;
	Information << " Kappa (-1: fully upwind scheme, 0: Fromm's scheme, 1: centeral differencing, and 1/3: third order upwind-biased scheme) = " << Kappa << endl;
	Information << "=======================" << endl;





	while (time < EndTime)
	{
		iter = iter + 1;
		time = time + dt;

		//=== numerical fluxes through all faces ===//
		NumericalFluxInALLNodes(K, num_flag, beta, nx1, ny1, nx2, ny2, nx3, ny3,
			u_flux, v_flux, rho_flux, p_flux, alpha_l_flux, EToE, EToF, //inputs
			NumFlux_1, NumFlux_2, NumFlux_3, NumFlux_4);


		//=== Change numerical fluxes for Boundary faces ===//
		ChangeNumericalFluxForBoundaryNodes(K, num_flag, beta, nx1, ny1, nx2, ny2, nx3, ny3,
			u_flux, v_flux, rho_flux, p_flux, alpha_l_flux, //inputs
			NBCs, BC_NumElem, BC_Type_Elem, BC_Type_Edge,
			uinf, vinf, rho_m_inf, pinf, alpha_l_inf,
			NumFlux_1, NumFlux_2, NumFlux_3, NumFlux_4);

		//=== Neighbor values ===//
		//for all real elements
		//if (num_flag == 7)
		//{



		parallel_for(1, K + 1, [&](int j)
			//for (int j = 1; j < K + 1; j++)
		{
			//for (int j = 1; j <= K; j++)
			//{
			Q_neighbour1[j][1] = Qo_1[EToE[j][1]];
			Q_neighbour2[j][1] = Qo_2[EToE[j][1]];
			Q_neighbour3[j][1] = Qo_3[EToE[j][1]];
			Q_neighbour4[j][1] = Qo_4[EToE[j][1]];

			Q_neighbour1[j][2] = Qo_1[EToE[j][2]];
			Q_neighbour2[j][2] = Qo_2[EToE[j][2]];
			Q_neighbour3[j][2] = Qo_3[EToE[j][2]];
			Q_neighbour4[j][2] = Qo_4[EToE[j][2]];

			Q_neighbour1[j][3] = Qo_1[EToE[j][3]];
			Q_neighbour2[j][3] = Qo_2[EToE[j][3]];
			Q_neighbour3[j][3] = Qo_3[EToE[j][3]];
			Q_neighbour4[j][3] = Qo_4[EToE[j][3]];

			//cout <<"j = "<<j<< ", j(1) = " << EToE[j][1] << ", j(2) = " << EToE[j][2] << ", j(3) = " << EToE[j][3] << endl;

			//Qo_1_W_Ghost[j] = Qo_1[j];
			//Qo_2_W_Ghost[j] = Qo_2[j];
			//Qo_3_W_Ghost[j] = Qo_3[j];
			//Qo_4_W_Ghost[j] = Qo_4[j];
			//}
		});
		//}
		//revising for real boundary elements

		for (int BC_Counter = 1; BC_Counter <= NBCs; BC_Counter++)
		{
			parallel_for(1, BC_NumElem[BC_Counter] + 1, [&](int BC_Elem_Counter)
				//for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter] + 1; BC_Elem_Counter++)
			{
				double u_minus, v_minus, p_minus, rho_minus, alpha_l_minus;
				double uBC, vBC, pBC, rhoBC, alpha_lBC;

				int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter];
				int i = BC_Type_Edge[BC_Elem_Counter][BC_Counter];

				p_minus = Qo_1[j];
				u_minus = Qo_2[j];
				v_minus = Qo_3[j];
				alpha_l_minus = Qo_4[j];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;


				if (i == 1)
				{
					double nx = nx1[j];
					double ny = ny1[j];
					BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
						, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
						uBC, vBC, pBC, rhoBC, alpha_lBC);
				}
				else if (i == 2)
				{
					double nx = nx2[j];
					double ny = ny2[j];
					BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
						, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
						uBC, vBC, pBC, rhoBC, alpha_lBC);
				}
				else
				{
					double nx = nx3[j];
					double ny = ny3[j];
					BoundaryValues(nx, ny, BC_Counter, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
						, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
						uBC, vBC, pBC, rhoBC, alpha_lBC);
				}

				Q_neighbour1[j][i] = pBC;
				Q_neighbour2[j][i] = uBC;
				Q_neighbour3[j][i] = vBC;
				Q_neighbour4[j][i] = alpha_lBC;
			});
			//}
		}


		//Values for ghost elements
		for (int GEC = K + 1; GEC <= GhostElementCounter; GEC++)
		{
			double u_minus, v_minus, p_minus, rho_minus, alpha_l_minus;
			double uBC, vBC, pBC, rhoBC, alpha_lBC;

			int i = BC_Ghost_Edge[GEC];
			int j = BC_Ghost_Elem[GEC]; //elemane vagheii ke ba aan dar ertebat ast!
			int k = BC_Type_Ghost[GEC];

			if (i == 1)
			{
				double nx = nx1[j];
				double ny = ny1[j];

				Q_neighbour1[GEC][1] = Qo_1[j];
				Q_neighbour2[GEC][1] = Qo_2[j];
				Q_neighbour3[GEC][1] = Qo_3[j];
				Q_neighbour4[GEC][1] = Qo_4[j];

				p_minus = Qo_1[EToE_Ghost[GEC][2]];
				u_minus = Qo_2[EToE_Ghost[GEC][2]];
				v_minus = Qo_3[EToE_Ghost[GEC][2]];
				alpha_l_minus = Qo_4[EToE_Ghost[GEC][2]];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
				BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
				Q_neighbour1[GEC][2] = pBC;
				Q_neighbour2[GEC][2] = uBC;
				Q_neighbour3[GEC][2] = vBC;
				Q_neighbour4[GEC][2] = alpha_lBC;

				p_minus = Qo_1[EToE_Ghost[GEC][3]];
				u_minus = Qo_2[EToE_Ghost[GEC][3]];
				v_minus = Qo_3[EToE_Ghost[GEC][3]];
				alpha_l_minus = Qo_4[EToE_Ghost[GEC][3]];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
				BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
				Q_neighbour1[GEC][3] = pBC;
				Q_neighbour2[GEC][3] = uBC;
				Q_neighbour3[GEC][3] = vBC;
				Q_neighbour4[GEC][3] = alpha_lBC;

			}
			else if (i == 2)
			{
				double nx = nx2[j];
				double ny = ny2[j];

				Q_neighbour1[GEC][2] = Qo_1[j];
				Q_neighbour2[GEC][2] = Qo_2[j];
				Q_neighbour3[GEC][2] = Qo_3[j];
				Q_neighbour4[GEC][2] = Qo_4[j];

				p_minus = Qo_1[EToE_Ghost[GEC][1]];
				u_minus = Qo_2[EToE_Ghost[GEC][1]];
				v_minus = Qo_3[EToE_Ghost[GEC][1]];
				alpha_l_minus = Qo_4[EToE_Ghost[GEC][1]];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
				BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
				Q_neighbour1[GEC][1] = pBC;
				Q_neighbour2[GEC][1] = uBC;
				Q_neighbour3[GEC][1] = vBC;
				Q_neighbour4[GEC][1] = alpha_lBC;

				p_minus = Qo_1[EToE_Ghost[GEC][3]];
				u_minus = Qo_2[EToE_Ghost[GEC][3]];
				v_minus = Qo_3[EToE_Ghost[GEC][3]];
				alpha_l_minus = Qo_4[EToE_Ghost[GEC][3]];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
				BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
				Q_neighbour1[GEC][3] = pBC;
				Q_neighbour2[GEC][3] = uBC;
				Q_neighbour3[GEC][3] = vBC;
				Q_neighbour4[GEC][3] = alpha_lBC;
			}
			else
			{
				double nx = nx3[j];
				double ny = ny3[j];

				Q_neighbour1[GEC][3] = Qo_1[j];
				Q_neighbour2[GEC][3] = Qo_2[j];
				Q_neighbour3[GEC][3] = Qo_3[j];
				Q_neighbour4[GEC][3] = Qo_4[j];

				p_minus = Qo_1[EToE_Ghost[GEC][2]];
				u_minus = Qo_2[EToE_Ghost[GEC][2]];
				v_minus = Qo_3[EToE_Ghost[GEC][2]];
				alpha_l_minus = Qo_4[EToE_Ghost[GEC][2]];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
				BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
				Q_neighbour1[GEC][2] = pBC;
				Q_neighbour2[GEC][2] = uBC;
				Q_neighbour3[GEC][2] = vBC;
				Q_neighbour4[GEC][2] = alpha_lBC;

				p_minus = Qo_1[EToE_Ghost[GEC][1]];
				u_minus = Qo_2[EToE_Ghost[GEC][1]];
				v_minus = Qo_3[EToE_Ghost[GEC][1]];
				alpha_l_minus = Qo_4[EToE_Ghost[GEC][1]];
				rho_minus = alpha_l_minus + (1 - alpha_l_minus)*rho_v;
				BoundaryValues(nx, ny, k, uinf, vinf, pinf, rho_m_inf, alpha_l_inf
					, u_minus, v_minus, p_minus, rho_minus, alpha_l_minus,
					uBC, vBC, pBC, rhoBC, alpha_lBC);
				Q_neighbour1[GEC][1] = pBC;
				Q_neighbour2[GEC][1] = uBC;
				Q_neighbour3[GEC][1] = vBC;
				Q_neighbour4[GEC][1] = alpha_lBC;
			}

		}



		//for (int j = 1; j <= K + Total_BC_Elem; j++)
		//{
		//	for (int i = 1; i <= Nfaces; i++)
		//	{
		//		cout << "j = " << j << ", i = " << i << ", Q_neighbour1[j][i] = " << Q_neighbour1[j][i] << endl;
		//	}

		//}


		//=== Dissipation term ===//
		parallel_for(1, K + 1, [&](int j)
			//for (int j = 1; j < K + 1; j++)
		{

			//if (j == 4929)
			//{
			//	cout << "j = " << j << ", Q1 = " << Qn_1[j] << ", Q2 = " << Qn_2[j] << ", Q3 = " << Qn_3[j] << ", Q4 = " << Qn_4[j] << endl;
			//}



			int k;
			double d1_1, d2_1, d3_1, d4_1, d1_2, d2_2, d3_2, d4_2, d1_3, d2_3, d3_3, d4_3;
			double nu_k, gamma_k, eps2, eps4;
			k = 1;
			nu_k = abs((p[j] - p[EToE[j][k]]) / (p[j] + p[EToE[j][k]]));
			gamma_k = abs((rho[j] - rho[EToE[j][k]]) / (rho[j] + rho[EToE[j][k]]));
			eps2 = k2p*nu_k + k2rho*gamma_k;
			eps4 = maxval(0, k4 - eps2);
			d1_1 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_1[j], Q_neighbour1, j) - del2Q(Nfaces, Q_neighbour1[j][k], Q_neighbour1, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_1[j] - Q_neighbour1[j][k])) / dt;
			d2_1 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_2[j], Q_neighbour2, j) - del2Q(Nfaces, Q_neighbour2[j][k], Q_neighbour2, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_2[j] - Q_neighbour2[j][k])) / dt;
			d3_1 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_3[j], Q_neighbour3, j) - del2Q(Nfaces, Q_neighbour3[j][k], Q_neighbour3, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_3[j] - Q_neighbour3[j][k])) / dt;
			d4_1 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_4[j], Q_neighbour4, j) - del2Q(Nfaces, Q_neighbour4[j][k], Q_neighbour4, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_4[j] - Q_neighbour4[j][k])) / dt;

			k = 2;
			nu_k = abs((p[j] - p[EToE[j][k]]) / (p[j] + p[EToE[j][k]]));
			gamma_k = abs((rho[j] - rho[EToE[j][k]]) / (rho[j] + rho[EToE[j][k]]));
			eps2 = k2p*nu_k + k2rho*gamma_k;
			eps4 = maxval(0, k4 - eps2);
			d1_2 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_1[j], Q_neighbour1, j) - del2Q(Nfaces, Q_neighbour1[j][k], Q_neighbour1, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_1[j] - Q_neighbour1[j][k])) / dt;
			d2_2 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_2[j], Q_neighbour2, j) - del2Q(Nfaces, Q_neighbour2[j][k], Q_neighbour2, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_2[j] - Q_neighbour2[j][k])) / dt;
			d3_2 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_3[j], Q_neighbour3, j) - del2Q(Nfaces, Q_neighbour3[j][k], Q_neighbour3, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_3[j] - Q_neighbour3[j][k])) / dt;
			d4_2 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_4[j], Q_neighbour4, j) - del2Q(Nfaces, Q_neighbour4[j][k], Q_neighbour4, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_4[j] - Q_neighbour4[j][k])) / dt;

			k = 3;
			nu_k = abs((p[j] - p[EToE[j][k]]) / (p[j] + p[EToE[j][k]]));
			gamma_k = abs((rho[j] - rho[EToE[j][k]]) / (rho[j] + rho[EToE[j][k]]));
			eps2 = k2p*nu_k + k2rho*gamma_k;
			eps4 = maxval(0, k4 - eps2);
			d1_3 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_1[j], Q_neighbour1, j) - del2Q(Nfaces, Q_neighbour1[j][k], Q_neighbour1, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_1[j] - Q_neighbour1[j][k])) / dt;
			d2_3 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_2[j], Q_neighbour2, j) - del2Q(Nfaces, Q_neighbour2[j][k], Q_neighbour2, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_2[j] - Q_neighbour2[j][k])) / dt;
			d3_3 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_3[j], Q_neighbour3, j) - del2Q(Nfaces, Q_neighbour3[j][k], Q_neighbour3, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_3[j] - Q_neighbour3[j][k])) / dt;
			d4_3 = (Area[j] + Area[EToE[j][k]]) *(0.5*eps4*(del2Q(Nfaces, Qo_4[j], Q_neighbour4, j) - del2Q(Nfaces, Q_neighbour4[j][k], Q_neighbour4, EToE_Ghost[j][k])) - 0.5*eps2*(Qo_4[j] - Q_neighbour4[j][k])) / dt;

			D1[j] = d1_1 + d1_2 + d1_3;
			D2[j] = d2_1 + d2_2 + d2_3;
			D3[j] = d3_1 + d3_2 + d3_3;
			D4[j] = d4_1 + d4_2 + d4_3;

		});
		//}
		//}

		//=== Finding new variables ===//
		//From Hejranfar and Ezzatneshan paper
		parallel_for(1, K + 1, [&](int j)
		{
			double BounderyIntegral1 = NumFlux_1[j][1] * LengFace[j][1] + NumFlux_1[j][2] * LengFace[j][2] + NumFlux_1[j][3] * LengFace[j][3];
			double BounderyIntegral2 = NumFlux_2[j][1] * LengFace[j][1] + NumFlux_2[j][2] * LengFace[j][2] + NumFlux_2[j][3] * LengFace[j][3];
			double BounderyIntegral3 = NumFlux_3[j][1] * LengFace[j][1] + NumFlux_3[j][2] * LengFace[j][2] + NumFlux_3[j][3] * LengFace[j][3];
			double BounderyIntegral4 = NumFlux_4[j][1] * LengFace[j][1] + NumFlux_4[j][2] * LengFace[j][2] + NumFlux_4[j][3] * LengFace[j][3];

			R_1[j] = P_Matrix_inv11[j] * (BounderyIntegral1 - Source1[j] * Area[j]) + P_Matrix_inv12[j] * (BounderyIntegral2 - Source2[j] * Area[j])
				+ P_Matrix_inv13[j] * (BounderyIntegral3 - Source3[j] * Area[j]) + P_Matrix_inv14[j] * (BounderyIntegral4 - Source4[j] * Area[j]);
			R_2[j] = P_Matrix_inv21[j] * (BounderyIntegral1 - Source1[j] * Area[j]) + P_Matrix_inv22[j] * (BounderyIntegral2 - Source2[j] * Area[j])
				+ P_Matrix_inv23[j] * (BounderyIntegral3 - Source3[j] * Area[j]) + P_Matrix_inv24[j] * (BounderyIntegral4 - Source4[j] * Area[j]);
			R_3[j] = P_Matrix_inv31[j] * (BounderyIntegral1 - Source1[j] * Area[j]) + P_Matrix_inv32[j] * (BounderyIntegral2 - Source2[j] * Area[j])
				+ P_Matrix_inv33[j] * (BounderyIntegral3 - Source3[j] * Area[j]) + P_Matrix_inv34[j] * (BounderyIntegral4 - Source4[j] * Area[j]);
			R_4[j] = P_Matrix_inv41[j] * (BounderyIntegral1 - Source1[j] * Area[j]) + P_Matrix_inv42[j] * (BounderyIntegral2 - Source2[j] * Area[j])
				+ P_Matrix_inv43[j] * (BounderyIntegral3 - Source3[j] * Area[j]) + P_Matrix_inv44[j] * (BounderyIntegral4 - Source4[j] * Area[j]);

			Qn_1[j] = Qo_1[j] + (-R_1[j] + D1[j]) * dt / Area[j];
			Qn_2[j] = Qo_2[j] + (-R_2[j] + D2[j]) * dt / Area[j];
			Qn_3[j] = Qo_3[j] + (-R_3[j] + D3[j]) * dt / Area[j];
			Qn_4[j] = Qo_4[j] + (-R_4[j] + D4[j]) * dt / Area[j];

			//if (j == 4929)
			//{
			//	cout << "j = " << j << ", Q1 = " << Qn_1[j] << ", Q2 = " << Qn_2[j] << ", Q3 = " << Qn_3[j] << ", Q4 = " << Qn_4[j] << endl;
			//}


		});


		//=== Error calculating ===//
		error = 0;
		for (int j = 1; j < K + 1; j++)
		{
			error = error + (Qn_1[j] - Qo_1[j])*(Qn_1[j] - Qo_1[j]);
		}
		error = sqrt(error / K);


		//=== Update variables for next iteration ===//
		UpdateVariablesForNextIteration(K, Qn_1, Qn_2, Qn_3, Qn_4, C_prod, C_dest, tinf, p_v, rho_v, //inputs
			rho, u, v, p, alpha_l, m_dot_minus, m_dot_plus, Qo_1, Qo_2, Qo_3, Qo_4, F1, F2, F3, F4, G1, G2, G3, G4,
			Source1, Source2, Source3, Source4);

		InversionOfPreconditioningMatrix(K, u, v, rho, alpha_l, beta, drho_l, //nputs
			P_Matrix_inv11, P_Matrix_inv12, P_Matrix_inv13, P_Matrix_inv14,
			P_Matrix_inv21, P_Matrix_inv22, P_Matrix_inv23, P_Matrix_inv24,
			P_Matrix_inv31, P_Matrix_inv32, P_Matrix_inv33, P_Matrix_inv34,
			P_Matrix_inv41, P_Matrix_inv42, P_Matrix_inv43, P_Matrix_inv44);

		VariablesInVertices(Nv, Num_NodeCell, NodeNode, NodeCell, Area, u, v, p, rho, alpha_l,
			u_node, v_node, p_node, rho_node, alpha_l_node);

		NeighborValues(K, EToE, NBCs, BC_NumElem, BC_Type_Elem, rho_v, nx1, nx2, nx3, ny1, ny2, ny3,
			uinf, vinf, pinf, rho_m_inf, alpha_l_inf, GhostElementCounter,
			BC_Ghost_Edge, BC_Ghost_Elem, BC_Type_Ghost, EToE_Ghost,
			BC_Type_Edge, Qo_1, Qo_2, Qo_3, Qo_4,
			Q_neighbour1, Q_neighbour2, Q_neighbour3, Q_neighbour4);

		if (HighOrderMethod == "MUSCL")
		{
			PrimitivVariablesForFluxComputation_MUSCL(Kappa, eta, K, rho_v, EToE, EToF, EToV, a,
				Qo_1, Qo_2, Qo_3, Qo_4, Q_neighbour1, Q_neighbour2, Q_neighbour3, Q_neighbour4,
				u_node, v_node, p_node, rho_node, alpha_l_node,
				u_flux, v_flux, p_flux, rho_flux, alpha_l_flux);
		}
		else //Frink
		{
			PrimitivVariablesForFluxComputation_Frink(eta, K, EToV, u, v, p, rho, alpha_l,
				u_node, v_node, p_node, rho_node, alpha_l_node,
				u_flux, v_flux, p_flux, rho_flux, alpha_l_flux);
		}

		//=== Screen reports ===//
		if (iter % 200 == 0)
		{
			cout << "time = " << time << " , error = " << error << endl;
			Error_Plot << iter << "	" << error << endl;
		}



		if (iter % 10000 == 0)
		{
			cout << "***************" << endl;
			cout << "Ploting Contour" << endl;
			cout << "***************" << endl;
			ContourPlot(Nv, K, Vx, Vy, u, v, p, rho, alpha_l, m_dot_minus, m_dot_plus, Source1, Source4, EToV);




			ofstream Wall_plot("Wall.plt");
			//Wall_plot << "variables = \"x\" , \"y\" , \"u\" , \"v\" ,\"p\"  ,\"rho\"  ,\"alpha_l\"  ,\"m_dot_minus\"  ,\"m_dot_plus\"  ,\"Source1\"  ,\"Source4\" " << endl;
			Wall_plot << "variables = \"x\" , \"y\" , \"u\" , \"v\" ,\"p\"  ,\"rho\"  ,\"alpha_l\"  ,\"Cp\" " << endl;

			int BC_Counter_new = 3;
			Wall_plot << "zone T=\"Lower\"" << endl;
			for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter_new] + 1; BC_Elem_Counter++)
			{
				int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter_new];
				int i = BC_Type_Edge[BC_Elem_Counter][BC_Counter_new];
				/*Wall_plot << (Vx[EToV[j][1]] + Vx[EToV[j][1]] + Vx[EToV[j][1]]) / 3 << "	" << (Vy[EToV[j][1]] + Vy[EToV[j][1]] + Vy[EToV[j][1]]) / 3 << "	" << u[j] << "	" << v[j] << "	" << p[j] << "	" << rho[j] << "	" << alpha_l[j] << "	" << m_dot_minus[j] << "	" << m_dot_plus[j] << "	" << Source1[j] << "	" << Source4[j] << endl;*/
				if (i == 1)
				{
					i = 1;
					double p_wall = p_node[EToV[j][i]];
					double u_wall = u_node[EToV[j][i]];
					double v_wall = v_node[EToV[j][i]];
					double alpha_l_wall = alpha_l_node[EToV[j][i]];
					double rho_wall = rho_node[EToV[j][i]];
					double Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					i = 2;
					p_wall = p_node[EToV[j][i]];
					u_wall = u_node[EToV[j][i]];
					v_wall = v_node[EToV[j][i]];
					alpha_l_wall = alpha_l_node[EToV[j][i]];
					rho_wall = rho_node[EToV[j][i]];
					Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
				}
				else if (i == 2)
				{

					i = 2;
					double p_wall = p_node[EToV[j][i]];
					double u_wall = u_node[EToV[j][i]];
					double v_wall = v_node[EToV[j][i]];
					double alpha_l_wall = alpha_l_node[EToV[j][i]];
					double rho_wall = rho_node[EToV[j][i]];
					double Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					i = 3;
					p_wall = p_node[EToV[j][i]];
					u_wall = u_node[EToV[j][i]];
					v_wall = v_node[EToV[j][i]];
					alpha_l_wall = alpha_l_node[EToV[j][i]];
					rho_wall = rho_node[EToV[j][i]];
					Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
				}
				else if (i == 3)
				{
					i = 3;
					double p_wall = p_node[EToV[j][i]];
					double u_wall = u_node[EToV[j][i]];
					double v_wall = v_node[EToV[j][i]];
					double alpha_l_wall = alpha_l_node[EToV[j][i]];
					double rho_wall = rho_node[EToV[j][i]];
					double Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					i = 1;
					p_wall = p_node[EToV[j][i]];
					u_wall = u_node[EToV[j][i]];
					v_wall = v_node[EToV[j][i]];
					alpha_l_wall = alpha_l_node[EToV[j][i]];
					rho_wall = rho_node[EToV[j][i]];
					Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					//Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_node[EToV[j][i]] << "	" << v_node[EToV[j][i]] << "	" << p_node[EToV[j][i]] << "	" << rho_node[EToV[j][i]] << "	" << alpha_l_node[EToV[j][i]] << endl;
				}

			}
			BC_Counter_new = 4;
			Wall_plot << "zone T=\"Upper\"" << endl;
			for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter_new] + 1; BC_Elem_Counter++)
			{
				int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter_new];
				int i = BC_Type_Edge[BC_Elem_Counter][BC_Counter_new];
				//Wall_plot << (Vx[EToV[j][1]] + Vx[EToV[j][1]] + Vx[EToV[j][1]]) / 3 << "	" << (Vy[EToV[j][1]] + Vy[EToV[j][1]] + Vy[EToV[j][1]]) / 3 << "	" << u[j] << "	" << v[j] << "	" << p[j] << "	" << rho[j] << "	" << alpha_l[j] << "	" << m_dot_minus[j] << "	" << m_dot_plus[j] << "	" << Source1[j] << "	" << Source4[j] << endl;
				if (i == 1)
				{
					i = 1;
					double p_wall = p_node[EToV[j][i]];
					double u_wall = u_node[EToV[j][i]];
					double v_wall = v_node[EToV[j][i]];
					double alpha_l_wall = alpha_l_node[EToV[j][i]];
					double rho_wall = rho_node[EToV[j][i]];
					double Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					i = 2;
					p_wall = p_node[EToV[j][i]];
					u_wall = u_node[EToV[j][i]];
					v_wall = v_node[EToV[j][i]];
					alpha_l_wall = alpha_l_node[EToV[j][i]];
					rho_wall = rho_node[EToV[j][i]];
					Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
				}
				else if (i == 2)
				{
					i = 2;
					double p_wall = p_node[EToV[j][i]];
					double u_wall = u_node[EToV[j][i]];
					double v_wall = v_node[EToV[j][i]];
					double alpha_l_wall = alpha_l_node[EToV[j][i]];
					double rho_wall = rho_node[EToV[j][i]];
					double Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					i = 3;
					p_wall = p_node[EToV[j][i]];
					u_wall = u_node[EToV[j][i]];
					v_wall = v_node[EToV[j][i]];
					alpha_l_wall = alpha_l_node[EToV[j][i]];
					rho_wall = rho_node[EToV[j][i]];
					Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
				}
				else if (i == 3)
				{
					i = 3;
					double p_wall = p_node[EToV[j][i]];
					double u_wall = u_node[EToV[j][i]];
					double v_wall = v_node[EToV[j][i]];
					double alpha_l_wall = alpha_l_node[EToV[j][i]];
					double rho_wall = rho_node[EToV[j][i]];
					double Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
					i = 1;
					p_wall = p_node[EToV[j][i]];
					u_wall = u_node[EToV[j][i]];
					v_wall = v_node[EToV[j][i]];
					alpha_l_wall = alpha_l_node[EToV[j][i]];
					rho_wall = rho_node[EToV[j][i]];
					Cp_wall = (p_wall - 1) * 2;
					Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
				}

			}






		}

		if (error < 1e-8)
		{
			break;
		}

	}

	////===========================//
	////=== Write final results ===//
	////===========================//
	ContourPlot(Nv, K, Vx, Vy, u, v, p, rho, alpha_l, m_dot_minus, m_dot_plus, Source1, Source4, EToV);


	ofstream Wall_plot("Wall.plt");
	//Wall_plot << "variables = \"x\" , \"y\" , \"u\" , \"v\" ,\"p\"  ,\"rho\"  ,\"alpha_l\"  ,\"m_dot_minus\"  ,\"m_dot_plus\"  ,\"Source1\"  ,\"Source4\" " << endl;
	Wall_plot << "variables = \"x\" , \"y\" , \"u\" , \"v\" ,\"p\"  ,\"rho\"  ,\"alpha_l\"  ,\"Cp\" " << endl;

	int BC_Counter_new = 3;
	Wall_plot << "zone T=\"Lower\"" << endl;
	for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter_new] + 1; BC_Elem_Counter++)
	{
		int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter_new];
		int i = BC_Type_Edge[BC_Elem_Counter][BC_Counter_new];
		/*Wall_plot << (Vx[EToV[j][1]] + Vx[EToV[j][1]] + Vx[EToV[j][1]]) / 3 << "	" << (Vy[EToV[j][1]] + Vy[EToV[j][1]] + Vy[EToV[j][1]]) / 3 << "	" << u[j] << "	" << v[j] << "	" << p[j] << "	" << rho[j] << "	" << alpha_l[j] << "	" << m_dot_minus[j] << "	" << m_dot_plus[j] << "	" << Source1[j] << "	" << Source4[j] << endl;*/
		if (i == 1)
		{
			i = 1;
			double p_wall = p_node[EToV[j][i]];
			double u_wall = u_node[EToV[j][i]];
			double v_wall = v_node[EToV[j][i]];
			double alpha_l_wall = alpha_l_node[EToV[j][i]];
			double rho_wall = rho_node[EToV[j][i]];
			double Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			i = 2;
			p_wall = p_node[EToV[j][i]];
			u_wall = u_node[EToV[j][i]];
			v_wall = v_node[EToV[j][i]];
			alpha_l_wall = alpha_l_node[EToV[j][i]];
			rho_wall = rho_node[EToV[j][i]];
			Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
		}
		else if (i == 2)
		{

			i = 2;
			double p_wall = p_node[EToV[j][i]];
			double u_wall = u_node[EToV[j][i]];
			double v_wall = v_node[EToV[j][i]];
			double alpha_l_wall = alpha_l_node[EToV[j][i]];
			double rho_wall = rho_node[EToV[j][i]];
			double Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			i = 3;
			p_wall = p_node[EToV[j][i]];
			u_wall = u_node[EToV[j][i]];
			v_wall = v_node[EToV[j][i]];
			alpha_l_wall = alpha_l_node[EToV[j][i]];
			rho_wall = rho_node[EToV[j][i]];
			Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
		}
		else if (i == 3)
		{
			i = 3;
			double p_wall = p_node[EToV[j][i]];
			double u_wall = u_node[EToV[j][i]];
			double v_wall = v_node[EToV[j][i]];
			double alpha_l_wall = alpha_l_node[EToV[j][i]];
			double rho_wall = rho_node[EToV[j][i]];
			double Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			i = 1;
			p_wall = p_node[EToV[j][i]];
			u_wall = u_node[EToV[j][i]];
			v_wall = v_node[EToV[j][i]];
			alpha_l_wall = alpha_l_node[EToV[j][i]];
			rho_wall = rho_node[EToV[j][i]];
			Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			//Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_node[EToV[j][i]] << "	" << v_node[EToV[j][i]] << "	" << p_node[EToV[j][i]] << "	" << rho_node[EToV[j][i]] << "	" << alpha_l_node[EToV[j][i]] << endl;
		}

	}
	BC_Counter_new = 4;
	Wall_plot << "zone T=\"Upper\"" << endl;
	for (int BC_Elem_Counter = 1; BC_Elem_Counter < BC_NumElem[BC_Counter_new] + 1; BC_Elem_Counter++)
	{
		int j = BC_Type_Elem[BC_Elem_Counter][BC_Counter_new];
		int i = BC_Type_Edge[BC_Elem_Counter][BC_Counter_new];
		//Wall_plot << (Vx[EToV[j][1]] + Vx[EToV[j][1]] + Vx[EToV[j][1]]) / 3 << "	" << (Vy[EToV[j][1]] + Vy[EToV[j][1]] + Vy[EToV[j][1]]) / 3 << "	" << u[j] << "	" << v[j] << "	" << p[j] << "	" << rho[j] << "	" << alpha_l[j] << "	" << m_dot_minus[j] << "	" << m_dot_plus[j] << "	" << Source1[j] << "	" << Source4[j] << endl;
		if (i == 1)
		{
			i = 1;
			double p_wall = p_node[EToV[j][i]];
			double u_wall = u_node[EToV[j][i]];
			double v_wall = v_node[EToV[j][i]];
			double alpha_l_wall = alpha_l_node[EToV[j][i]];
			double rho_wall = rho_node[EToV[j][i]];
			double Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			i = 2;
			p_wall = p_node[EToV[j][i]];
			u_wall = u_node[EToV[j][i]];
			v_wall = v_node[EToV[j][i]];
			alpha_l_wall = alpha_l_node[EToV[j][i]];
			rho_wall = rho_node[EToV[j][i]];
			Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
		}
		else if (i == 2)
		{
			i = 2;
			double p_wall = p_node[EToV[j][i]];
			double u_wall = u_node[EToV[j][i]];
			double v_wall = v_node[EToV[j][i]];
			double alpha_l_wall = alpha_l_node[EToV[j][i]];
			double rho_wall = rho_node[EToV[j][i]];
			double Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			i = 3;
			p_wall = p_node[EToV[j][i]];
			u_wall = u_node[EToV[j][i]];
			v_wall = v_node[EToV[j][i]];
			alpha_l_wall = alpha_l_node[EToV[j][i]];
			rho_wall = rho_node[EToV[j][i]];
			Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
		}
		else if (i == 3)
		{
			i = 3;
			double p_wall = p_node[EToV[j][i]];
			double u_wall = u_node[EToV[j][i]];
			double v_wall = v_node[EToV[j][i]];
			double alpha_l_wall = alpha_l_node[EToV[j][i]];
			double rho_wall = rho_node[EToV[j][i]];
			double Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
			i = 1;
			p_wall = p_node[EToV[j][i]];
			u_wall = u_node[EToV[j][i]];
			v_wall = v_node[EToV[j][i]];
			alpha_l_wall = alpha_l_node[EToV[j][i]];
			rho_wall = rho_node[EToV[j][i]];
			Cp_wall = (p_wall - 1) * 2;
			Wall_plot << Vx[EToV[j][i]] << "	" << Vy[EToV[j][i]] << "	" << u_wall << "	" << v_wall << "	" << p_wall << "	" << rho_wall << "	" << alpha_l_wall << "	" << Cp_wall << endl;
		}

	}


	t = clock() - t;
	double CPU_time = double(t) / CLOCKS_PER_SEC;

	////write all informations


	ofstream Information2("Information.txt");
	Information2 << " *** Grid informations ***" << endl;
	Information2 << " Total Elements (K) = " << K << endl;
	Information2 << " Total global nodes (Nv) = " << Nv << endl;
	Information2 << " Number of specified boundary conditions = " << NBCs << endl;
	Information2 << " Total boundary elements = " << Total_BC_Elem << endl;
	Information2 << "=======================" << endl;
	Information2 << " *** All information about computational cost ***" << endl;
	Information2 << " End time = " << EndTime << endl;
	Information2 << " dt = " << dt << endl;
	Information2 << " Number of iterations = " << iter << endl;
	Information2 << " Convergence criteria = " << error << endl;
	Information2 << " Converged time = " << time << endl;
	Information2 << " Cpu time = " << CPU_time << endl;
	Information2 << "=======================" << endl;
	Information2 << " *** Informations about flow conditions ***" << endl;
	Information2 << " beta = " << beta << endl;
	Information2 << " rho_v = " << rho_v << endl;
	Information2 << " Cavitation number = " << 2 * (pinf - p_v) << endl;
	Information2 << " C_dest = " << C_dest << endl;
	Information2 << " C_prod = " << C_prod << endl;
	Information2 << " Angle of attack = " << AOA << endl;
	Information2 << "=======================" << endl;
	Information2 << " *** Informations about solver ***" << endl;
	Information2 << " Numerical Flux ID (1= Roe, 2=Lax, 3=AUSM+, 4=HLL, 5=HLLC, 6=LDFSS, 7=Average) = " << num_flag << endl;
	Information2 << " k2p = " << k2p << endl;
	Information2 << " k2rho = " << k2rho << endl;
	Information2 << " k4 = " << k4 << endl;
	Information2 << " High order method = " << HighOrderMethod << endl;
	Information2 << " eta (0=1st order and 1=2nd order) = " << eta << endl;
	Information2 << " Kappa (-1: fully upwind scheme, 0: Fromm's scheme, 1: centeral differencing, and 1/3: third order upwind-biased scheme) = " << Kappa << endl;
	Information2 << "=======================" << endl;


	system("pause");
}
