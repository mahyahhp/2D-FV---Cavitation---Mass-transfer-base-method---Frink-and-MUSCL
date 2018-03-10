#include "Define_All_Includes_Header.h"
#include "Define_Arrays_Header.h"

void Jacobian_Matrices_NormalVectors(int const K, int **EToV,
	double *Vx, double const *Vy, double **LengFace,
	double *nx1, double *ny1, double *nx2, double *ny2, double *nx3, double *ny3)
{
	double *dr_dx, *dr_dy, *ds_dx, *ds_dy, *J, *J1, *J2, *J3;
	dr_dx = d1array(K + 1); dr_dy = d1array(K + 1); ds_dx = d1array(K + 1); ds_dy = d1array(K + 1);
	J = d1array(K + 1); J1 = d1array(K + 1); J2 = d1array(K + 1); J3 = d1array(K + 1);
		
	parallel_for(1, K + 1, [&](int j)
	{
		double dx_dr, dx_ds, dy_dr, dy_ds;
		dx_dr = -0.5*Vx[EToV[j][1]] + 0.5*Vx[EToV[j][2]];
			dx_ds = -0.5*Vx[EToV[j][1]] + 0.5*Vx[EToV[j][3]];
			dy_dr = -0.5*Vy[EToV[j][1]] + 0.5*Vy[EToV[j][2]];
			dy_ds = -0.5*Vy[EToV[j][1]] + 0.5*Vy[EToV[j][3]];
			J[j]= dx_dr* dy_ds - dx_ds* dy_dr;
			dr_dx[j] = dy_ds / J[j];
			dr_dy[j] = -dx_ds / J[j];
			ds_dx[j] = -dy_dr / J[j];
			ds_dy[j] = dx_dr / J[j];
	});

	parallel_for(1, K + 1, [&](int j)
	{
			//face 1
		double norm1, norm2, norm3;
			norm1 = sqrt(pow(ds_dx[j], 2) + pow(ds_dy[j], 2));
			nx1[j] = -ds_dx[j] / norm1;
			ny1[j] = -ds_dy[j] / norm1;
			//face 2
			norm2 = sqrt(pow(dr_dx[j] + ds_dx[j], 2) + pow(dr_dy[j] + ds_dy[j], 2));
			nx2[j] = (dr_dx[j] + ds_dx[j]) / norm2;
			ny2[j] = (dr_dy[j] + ds_dy[j]) / norm2;
			//face 3
			norm3 = sqrt(pow(dr_dx[j], 2) + pow(dr_dy[j], 2));
			nx3[j] = -dr_dx[j] / norm3;
			ny3[j] = -dr_dy[j] / norm3;

			//Jacobian computing
			J1[j] = LengFace[j][1] / 2;
			J2[j] = LengFace[j][2] / 2;
			J3[j] = LengFace[j][3] / 2;
	});
}
