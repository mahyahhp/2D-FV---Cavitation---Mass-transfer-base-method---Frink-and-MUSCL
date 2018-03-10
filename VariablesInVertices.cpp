#include "Define_All_Includes_Header.h"


void  VariablesInVertices(int const Nv, int *No_NodeCell, int **NodeNode, int **NodeCell, double *Area, double *u,
	double *v, double *p, double *rho, double *alpha_l,
	double *u_node, double *v_node, double *p_node, double *rho_node, double *alpha_l_node)
{
	parallel_for(1, Nv + 1, [&](int i)
	{
		double sum_u = 0, sum_v = 0, sum_p = 0, sum_rho = 0, sum_alpha_l=0, sum_area=0;
		for (int j = 1; j <= No_NodeCell[i]; j++)
		{
			//sum_u = sum_u + u[NodeCell[i][j]][NodeNode[i][j]];
			sum_u       = sum_u       + u      [NodeCell[i][j]] * Area[j];
			sum_v       = sum_v       + v      [NodeCell[i][j]] * Area[j];
			sum_p       = sum_p       + p      [NodeCell[i][j]] * Area[j];
			sum_rho     = sum_rho     + rho    [NodeCell[i][j]] * Area[j];
			sum_alpha_l = sum_alpha_l + alpha_l[NodeCell[i][j]] * Area[j];

			sum_area = sum_area + Area[j];
		}
		//u_vertices[i] = sum_u / No_NodeCell[i];
		u_node[i]       = sum_u       / sum_area;
		v_node[i]       = sum_v       / sum_area;
		p_node[i]       = sum_p       / sum_area;
		rho_node[i]     = sum_rho     / sum_area;
		alpha_l_node[i] = sum_alpha_l / sum_area;
	});
}