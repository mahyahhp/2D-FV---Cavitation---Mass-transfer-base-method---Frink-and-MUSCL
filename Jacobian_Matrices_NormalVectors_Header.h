#ifndef _Jacobian_Matrices_NormalVectors_
#define _Jacobian_Matrices_NormalVectors_
extern void Jacobian_Matrices_NormalVectors(int const K, int **EToV,
	double *Vx, double const *Vy, double **LengFace,
	double *nx1, double *ny1, double *nx2, double *ny2, double *nx3, double *ny3);
#endif