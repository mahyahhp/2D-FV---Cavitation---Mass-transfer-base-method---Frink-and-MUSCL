#include "Define_All_Includes_Header.h"
#include "FluxVectors_Header.h"
#include "SolutionVector_Header.h"
double minval(double const a, double const b)
{
	double ans;
	ans = a;
	if (b<ans) ans = b;
	return ans;
}
double maxval(double const a, double const b)
{
	double ans;
	ans = a;
	if (b>ans) ans = b;
	return ans;
}
double M_plus_1(double const M)
{
	double ans;
	ans = 0.5*(M + abs(M));
	return ans;
}
double M_minus_1(double const M)
{
	double ans;
	ans = 0.5*(M - abs(M));
	return ans;
}
double M_plus_2(double const M)
{
	double ans;
	//if (abs(M) < 1)
	//{
	ans = 0.25*(M + 1)*(M + 1);
	//}
	//else
	//{
	//	ans = M_plus_1(M);
	//}
	return ans;
}
double M_minus_2(double const M)
{
	double ans;
	//if (abs(M) < 1)
	//{
	ans = -0.25*(M - 1)*(M - 1);
	//}
	//else
	//{
	//	ans = M_minus_1(M);
	//}
	return ans;
}
double M_plus_4(double const M, double const beta)
{
	double ans;
	if (abs(M) >= 1)
	{
		ans = M_plus_1(M);
	}
	else
	{
		ans = M_plus_2(M)*(1 - 16 * beta*M_minus_2(M));
	}
	return ans;
}
double M_minus_4(double const M, double const beta)
{
	double ans;
	if (abs(M) >= 1)
	{
		ans = M_minus_1(M);
	}
	else
	{
		ans = M_minus_2(M)*(1 + 16 * beta*M_plus_2(M));
	}
	return ans;
}
double p_plus_5(double const M, double const alpha)
{
	double ans;
	if (abs(M) >= 1)
	{
		ans = M_plus_1(M) / M;
	}
	else
	{
		ans = M_plus_2(M)*((+2 - M) - 16 * alpha*M*M_minus_2(M));
	}
	return ans;
}

double p_minus_5(double const M, double const alpha)
{
	double ans;
	if (abs(M) >= 1)
	{
		ans = M_minus_1(M) / M;
	}
	else
	{
		ans = M_minus_2(M)*((-2 - M) + 16 * alpha*M*M_plus_2(M));
	}
	return ans;
}


double sign(double const a)
{
	double result;
	if (a < 0)
	{
		result = -1;
	}
	else if (a > 0)
	{
		result = 1;
	}
	else
	{
		result = 0;
	}
	return result;
}


void WaveSpeed(double const u_mean, double const v_mean, double const nx, double const ny, double const beta,
	double &WS1L, double &WS2L, double &WS3L)
{
	double Vn = u_mean*nx + v_mean*ny;
	WS1L = Vn + sqrt(Vn*Vn + beta*beta);
	WS2L = Vn - sqrt(Vn*Vn + beta*beta);
	WS3L = Vn;
}


void MatrixInversion4x4(double const a11, double const a12, double const a13, double const a14,
	double const a21, double const a22, double const a23, double const a24,
	double const a31, double const a32, double const a33, double const a34,
	double const a41, double const a42, double const a43, double const a44, //inputs
	double &b11, double &b12, double &b13, double &b14,
	double &b21, double &b22, double &b23, double &b24,
	double &b31, double &b32, double &b33, double &b34,
	double &b41, double &b42, double &b43, double &b44)
{
	double det;
	det = a11*a22*a33*a44 - a11*a22*a34*a43 - a11*a23*a32*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 - a11*a24*a33*a42 - a12*a21*a33*a44 + a12*a21*a34*a43 + a12*a23*a31*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 + a12*a24*a33*a41 + a13*a21*a32*a44 - a13*a21*a34*a42 - a13*a22*a31*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 - a13*a24*a32*a41 - a14*a21*a32*a43 + a14*a21*a33*a42 + a14*a22*a31*a43 - a14*a22*a33*a41 - a14*a23*a31*a42 + a14*a23*a32*a41;
	b11 = (a22*a33*a44 - a22*a34*a43 - a23*a32*a44 + a23*a34*a42 + a24*a32*a43 - a24*a33*a42) / (det + 1e-14);
	b21 = (a21*a34*a43 - a21*a33*a44 + a23*a31*a44 - a23*a34*a41 - a24*a31*a43 + a24*a33*a41) / (det + 1e-14);
	b31 = (a21*a32*a44 - a21*a34*a42 - a22*a31*a44 + a22*a34*a41 + a24*a31*a42 - a24*a32*a41) / (det + 1e-14);
	b41 = (a21*a33*a42 - a21*a32*a43 + a22*a31*a43 - a22*a33*a41 - a23*a31*a42 + a23*a32*a41) / (det + 1e-14);

	b12 = (a12*a34*a43 - a12*a33*a44 + a13*a32*a44 - a13*a34*a42 - a14*a32*a43 + a14*a33*a42) / (det + 1e-14);
	b22 = (a11*a33*a44 - a11*a34*a43 - a13*a31*a44 + a13*a34*a41 + a14*a31*a43 - a14*a33*a41) / (det + 1e-14);
	b32 = (a11*a34*a42 - a11*a32*a44 + a12*a31*a44 - a12*a34*a41 - a14*a31*a42 + a14*a32*a41) / (det + 1e-14);
	b42 = (a11*a32*a43 - a11*a33*a42 - a12*a31*a43 + a12*a33*a41 + a13*a31*a42 - a13*a32*a41) / (det + 1e-14);

	b13 = (a12*a23*a44 - a12*a24*a43 - a13*a22*a44 + a13*a24*a42 + a14*a22*a43 - a14*a23*a42) / (det + 1e-14);
	b23 = (a11*a24*a43 - a11*a23*a44 + a13*a21*a44 - a13*a24*a41 - a14*a21*a43 + a14*a23*a41) / (det + 1e-14);
	b33 = (a11*a22*a44 - a11*a24*a42 - a12*a21*a44 + a12*a24*a41 + a14*a21*a42 - a14*a22*a41) / (det + 1e-14);
	b43 = (a11*a23*a42 - a11*a22*a43 + a12*a21*a43 - a12*a23*a41 - a13*a21*a42 + a13*a22*a41) / (det + 1e-14);

	b14 = (a12*a24*a33 - a12*a23*a34 + a13*a22*a34 - a13*a24*a32 - a14*a22*a33 + a14*a23*a32) / (det + 1e-14);
	b24 = (a11*a23*a34 - a11*a24*a33 - a13*a21*a34 + a13*a24*a31 + a14*a21*a33 - a14*a23*a31) / (det + 1e-14);
	b34 = (a11*a24*a32 - a11*a22*a34 + a12*a21*a34 - a12*a24*a31 - a14*a21*a32 + a14*a22*a31) / (det + 1e-14);
	b44 = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31) / (det + 1e-14);


}



void PreconditioningMatrix(double const u, double const v, double const rho, double const alpha_l, double const beta, double const drho_l,
	double & P_Matrix_inv11, double & P_Matrix_inv12, double & P_Matrix_inv13, double & P_Matrix_inv14,
	double & P_Matrix_inv21, double & P_Matrix_inv22, double & P_Matrix_inv23, double & P_Matrix_inv24,
	double & P_Matrix_inv31, double & P_Matrix_inv32, double & P_Matrix_inv33, double & P_Matrix_inv34,
	double & P_Matrix_inv41, double & P_Matrix_inv42, double & P_Matrix_inv43, double & P_Matrix_inv44)
{

	P_Matrix_inv11 = 1 / (rho * beta*beta);
	P_Matrix_inv12 = 0;
	P_Matrix_inv13 = 0;
	P_Matrix_inv14 = 0;

	P_Matrix_inv21 = 0;
	P_Matrix_inv22 = rho;
	P_Matrix_inv23 = 0;
	P_Matrix_inv24 = u * drho_l;

	P_Matrix_inv31 = 0;
	P_Matrix_inv32 = 0;
	P_Matrix_inv33 = rho;
	P_Matrix_inv34 = v * drho_l;

	P_Matrix_inv41 = alpha_l / (rho * beta*beta);
	P_Matrix_inv42 = 0;
	P_Matrix_inv43 = 0;
	P_Matrix_inv44 = 1;

}

void numerical_flux(int const num_flag, double const nx, double const ny, double const beta,
	double const uL, double const vL, double const rhoL, double const pL, double const alpha_lL,
	double const uR, double const vR, double const rhoR, double const pR, double const alpha_lR,
	double &num_flux_1, double &num_flux_2, double &num_flux_3, double &num_flux_4)

{
	if (num_flag == 1)
	{
		//***********//
		//*** Roe ***//
		//***********//

		double p_mean, u_mean, v_mean, rho_mean, alpha_l_mean;
		double Vn, VnL, VnR; //normal velocity
		double SoundSpeed;
		double Eigenvalue_1, Eigenvalue_2, Eigenvalue_3, Eigenvalue_4, WaveStrenghts_1, WaveStrenghts_2, WaveStrenghts_3, WaveStrenghts_4;
		double Eigenvector_11, Eigenvector_12, Eigenvector_13, Eigenvector_14,
			Eigenvector_21, Eigenvector_22, Eigenvector_23, Eigenvector_24,
			Eigenvector_31, Eigenvector_32, Eigenvector_33, Eigenvector_34,
			Eigenvector_41, Eigenvector_42, Eigenvector_43, Eigenvector_44;
		double det;
		double inv_Eigenvector_11, inv_Eigenvector_12, inv_Eigenvector_13, inv_Eigenvector_14,
			inv_Eigenvector_21, inv_Eigenvector_22, inv_Eigenvector_23, inv_Eigenvector_24,
			inv_Eigenvector_31, inv_Eigenvector_32, inv_Eigenvector_33, inv_Eigenvector_34,
			inv_Eigenvector_41, inv_Eigenvector_42, inv_Eigenvector_43, inv_Eigenvector_44;

		double dQ_1, dQ_2, dQ_3, dQ_4;
		double fL_1, fL_2, fL_3, fL_4, fR_1, fR_2, fR_3, fR_4;
		double D1_1, D2_1, D3_1, D4_1;
		double eps = 1e-14;
		double F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L;
		double F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R;

		//Define mean Roe values
		p_mean = (sqrt(rhoL)*pL + sqrt(rhoR)*pR) / (sqrt(rhoL) + sqrt(rhoR));
		u_mean = (sqrt(rhoL)*uL + sqrt(rhoR)*uR) / (sqrt(rhoL) + sqrt(rhoR));
		v_mean = (sqrt(rhoL)*vL + sqrt(rhoR)*vR) / (sqrt(rhoL) + sqrt(rhoR));
		rho_mean = sqrt(rhoL)*sqrt(rhoR);
		alpha_l_mean = (sqrt(rhoL)*alpha_lL + sqrt(rhoR)*alpha_lR) / (sqrt(rhoL) + sqrt(rhoR));

		//p_mean = (pL + pR) / (2);
		//u_mean = (uL + uR) / (2);
		//v_mean = (vL + vR) / (2);
		//rho_mean = (rhoL + rhoR) / (2);
		//alpha_l_mean = (alpha_lL + alpha_lR) / (2);

		Vn = u_mean*nx + v_mean*ny;
		VnL = uL*nx + vL*ny;
		VnR = uR*nx + vR*ny;

		////double rho_v = 0.01;




		////*** Drived by I ***//
		Eigenvalue_1 = Vn;
		Eigenvalue_2 = Vn + sqrt(beta*beta + Vn*Vn);
		Eigenvalue_3 = Vn - sqrt(beta*beta + Vn*Vn);
		Eigenvalue_4 = Vn;

		///Eigenvectors
		Eigenvector_11 = 0;
		Eigenvector_21 = -ny / (nx + eps);
		Eigenvector_31 = 1;
		Eigenvector_41 = 0;

		Eigenvector_13 = -(sqrt(beta *beta + Vn*Vn)*beta*beta * (-alpha_l_mean)) / (beta *beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) + v_mean*Vn + eps);
		Eigenvector_12 = -(sqrt(beta *beta + Vn*Vn)*beta*beta * (-rho_mean)) / (beta *beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) + v_mean*Vn + eps);
		Eigenvector_22 = (sqrt(beta * beta + Vn * Vn)*nx - Vn * nx + u_mean) / (sqrt(beta * beta + Vn * Vn)*ny - Vn * ny + v_mean + eps);
		Eigenvector_32 = 1;
		Eigenvector_42 = 0;

		Eigenvector_14 = -(sqrt(beta * beta + Vn * Vn)*beta * beta * (-alpha_l_mean)) / (-beta * beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) - v_mean*Vn + eps);
		Eigenvector_13 = -(sqrt(beta * beta + Vn * Vn)*beta * beta * (-rho_mean)) / (-beta * beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) - v_mean*Vn + eps);
		Eigenvector_23 = (sqrt(beta * beta + Vn * Vn)*nx + Vn * nx - u_mean) / (sqrt(beta * beta + Vn * Vn)*ny + Vn * ny - v_mean + eps);
		Eigenvector_33 = 1;
		Eigenvector_43 = 0;

		Eigenvector_14 = 0;
		Eigenvector_24 = 0;
		Eigenvector_34 = 0;
		Eigenvector_44 = 1;

		//Calculating WaveStrenghts
		//inverse of Eigenvectors

		MatrixInversion4x4(Eigenvector_11, Eigenvector_12, Eigenvector_13, Eigenvector_14,
			Eigenvector_21, Eigenvector_22, Eigenvector_23, Eigenvector_24,
			Eigenvector_31, Eigenvector_32, Eigenvector_33, Eigenvector_34,
			Eigenvector_41, Eigenvector_42, Eigenvector_43, Eigenvector_44, //inputs
			inv_Eigenvector_11, inv_Eigenvector_12, inv_Eigenvector_13, inv_Eigenvector_14,
			inv_Eigenvector_21, inv_Eigenvector_22, inv_Eigenvector_23, inv_Eigenvector_24,
			inv_Eigenvector_31, inv_Eigenvector_32, inv_Eigenvector_33, inv_Eigenvector_34,
			inv_Eigenvector_41, inv_Eigenvector_42, inv_Eigenvector_43, inv_Eigenvector_44);







		////*** Drived by azizollahi ***//
		//double C_bar = sqrt(beta*beta + Vn*Vn);
		//Eigenvalue_1 = Vn;
		//Eigenvalue_2 = Vn + C_bar;
		//Eigenvalue_3 = Vn - C_bar;
		//Eigenvalue_4 = Vn;

		///////Eigenvectors
		//Eigenvector_11 = 0;
		//Eigenvector_21 = (-2 * beta*beta*ny) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_31 = (2 * beta*beta*nx) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_41 = 0;

		//Eigenvector_12 = (C_bar*beta*beta*rho_mean) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_22 = (u_mean*(Vn + C_bar) + beta*beta*nx) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_32 = (v_mean*(Vn + C_bar) + beta*beta*ny) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_42 = 0;

		//Eigenvector_13 = (-C_bar*beta*beta*rho_mean) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_23 = (u_mean*(Vn - C_bar) + beta*beta*nx) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_33 = (v_mean*(Vn - C_bar) + beta*beta*ny) / (2 * beta*beta*C_bar*C_bar);
		//Eigenvector_43 = 0;

		//Eigenvector_14 = 0;
		//Eigenvector_24 = 0;
		//Eigenvector_34 = 0;
		//Eigenvector_44 = 1;

		///////inverse of Eigenvectors
		//inv_Eigenvector_11 = (u_mean*ny - v_mean*nx) / (rho_mean);
		//inv_Eigenvector_21 = (C_bar - Vn) / (rho_mean);
		//inv_Eigenvector_31 = (-C_bar - Vn) / (rho_mean);
		//inv_Eigenvector_41 = 0;

		//inv_Eigenvector_12 = -v_mean*Vn - beta*beta*ny;
		//inv_Eigenvector_22 = beta*beta*nx;
		//inv_Eigenvector_32 = beta*beta*nx;
		//inv_Eigenvector_42 = 0;

		//inv_Eigenvector_13 = u_mean*Vn + beta*beta*nx;
		//inv_Eigenvector_23 = beta*beta*ny;
		//inv_Eigenvector_33 = beta*beta*ny;
		//inv_Eigenvector_43 = 0;

		//inv_Eigenvector_14 = 0;
		//inv_Eigenvector_24 = 0;
		//inv_Eigenvector_34 = 0;
		//inv_Eigenvector_44 = 2 * beta*beta*C_bar*C_bar;








		//primitive values differences
		dQ_1 = (pR - pL);
		dQ_2 = (uR - uL);
		dQ_3 = (vR - vL);
		dQ_4 = (alpha_lR - alpha_lL);
		//dQ_1 = (pR - pL)*minval(sqrt(u_mean*u_mean + v_mean*v_mean) / C_bar, 1);
		//dQ_2 = (uR - uL)*minval(sqrt(u_mean*u_mean + v_mean*v_mean) / C_bar, 1);
		//dQ_3 = (vR - vL)*minval(sqrt(u_mean*u_mean + v_mean*v_mean) / C_bar, 1);
		//dQ_4 = (alpha_lR - alpha_lL)*minval(sqrt(u_mean*u_mean + v_mean*v_mean) / C_bar, 1);

		//Fluxes
		FluxVectors(pL, uL, vL, rhoL, alpha_lL,  //inputs
			F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L);

		FluxVectors(pR, uR, vR, rhoR, alpha_lR,  //inputs
			F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R);

		fL_1 = F1L*nx + G1L*ny;
		fL_2 = F2L*nx + G2L*ny;
		fL_3 = F3L*nx + G3L*ny;
		fL_4 = F4L*nx + G4L*ny;

		fR_1 = F1R*nx + G1R*ny;
		fR_2 = F2R*nx + G2R*ny;
		fR_3 = F3R*nx + G3R*ny;
		fR_4 = F4R*nx + G4R*ny;

		WaveStrenghts_1 = inv_Eigenvector_11*dQ_1 + inv_Eigenvector_12*dQ_2 + inv_Eigenvector_13*dQ_3 + inv_Eigenvector_14*dQ_4;
		WaveStrenghts_2 = inv_Eigenvector_21*dQ_1 + inv_Eigenvector_22*dQ_2 + inv_Eigenvector_23*dQ_3 + inv_Eigenvector_24*dQ_4;
		WaveStrenghts_3 = inv_Eigenvector_31*dQ_1 + inv_Eigenvector_32*dQ_2 + inv_Eigenvector_33*dQ_3 + inv_Eigenvector_34*dQ_4;
		WaveStrenghts_4 = inv_Eigenvector_41*dQ_1 + inv_Eigenvector_42*dQ_2 + inv_Eigenvector_43*dQ_3 + inv_Eigenvector_44*dQ_4;

		//Preconditioning matrix
		double P11, P12, P13, P14, P21, P22, P23, P24, P31, P32, P33, P34, P41, P42, P43, P44, D1, D2, D3, D4;
		PreconditioningMatrix(u_mean, v_mean, rho_mean, alpha_l_mean, beta, 1 - 0.01,
			P11, P12, P13, P14,
			P21, P22, P23, P24,
			P31, P32, P33, P34,
			P41, P42, P43, P44);
		//PreconditioningMatrix(uL, vL, rhoL, alpha_lL, beta, 1-0.01,
		//	P11, P12, P13, P14,
		//	P21, P22, P23, P24,
		//	P31, P32, P33, P34,
		//	P41, P42, P43, P44);
		//PreconditioningMatrix(uR, vR, rhoR, alpha_lR, beta, 1 - 0.01,
		//	P11, P12, P13, P14,
		//	P21, P22, P23, P24,
		//	P31, P32, P33, P34,
		//	P41, P42, P43, P44);


		//Calculating Roe's dissipation and numerical fluxes
		D1_1 = Eigenvector_11 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_12 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_13 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_14 * abs(Eigenvalue_4)*WaveStrenghts_4;
		D2_1 = Eigenvector_21 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_22 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_23 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_24 * abs(Eigenvalue_4)*WaveStrenghts_4;
		D3_1 = Eigenvector_31 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_32 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_33 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_34 * abs(Eigenvalue_4)*WaveStrenghts_4;
		D4_1 = Eigenvector_41 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_42 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_43 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_44 * abs(Eigenvalue_4)*WaveStrenghts_4;

		//*** Classical Roe Scheme ***
		D1 = D1_1*P11 + D2_1*P12 + D3_1*P13 + D4_1*P14;
		D2 = D1_1*P21 + D2_1*P22 + D3_1*P23 + D4_1*P24;
		D3 = D1_1*P31 + D2_1*P32 + D3_1*P33 + D4_1*P34;
		D4 = D1_1*P41 + D2_1*P42 + D3_1*P43 + D4_1*P44;

		////*** Low-Speed-Roe scheme ***//
		//D1_1 = abs(Vn)*dQ_1;
		//D2_1 = abs(Vn)*dQ_2;
		//D3_1 = abs(Vn)*dQ_3;
		//D4_1 = abs(Vn)*dQ_4;

		//D1 = D1_1*P11 + D2_1*P12 + D3_1*P13 + D4_1*P14;
		//D2 = D1_1*P21 + D2_1*P22 + D3_1*P23 + D4_1*P24;
		//D3 = D1_1*P31 + D2_1*P32 + D3_1*P33 + D4_1*P34;
		//D4 = D1_1*P41 + D2_1*P42 + D3_1*P43 + D4_1*P44;

		num_flux_1 = 0.5 * (fL_1 + fR_1) - 0.5 * D1;
		num_flux_2 = 0.5 * (fL_2 + fR_2) - 0.5 * D2;
		num_flux_3 = 0.5 * (fL_3 + fR_3) - 0.5 * D3;
		num_flux_4 = 0.5 * (fL_4 + fR_4) - 0.5 * D4;



		////************************************//
		////*** Roe all-speed numerical flux ***//
		////************************************//
		//double c2 = 0.05; //should be larger than a threshold value 0.04
		//double Uc = 0.5*(VnL + VnR) - c2 / (1 * 1)*(pR - pL);
		//double C_bar = sqrt(beta*beta + Vn*Vn);

		//double Mach = sqrt(Vn*Vn)/C_bar;
		//double f_M = minval(Mach*sqrt(4 + (1 - Mach*Mach)*(1 - Mach*Mach)) / (1 + Mach*Mach),1);
		//f_M = 0.0;

		//Eigenvalue_1 = Vn;
		//Eigenvalue_2 = Vn + f_M*C_bar;
		//Eigenvalue_3 = Vn - f_M*C_bar;
		//Eigenvalue_4 = Vn;

		///////Eigenvectors

		////////Eigenvalues # 5
		//////Eigenvalue_1 = Vn;
		//////Eigenvalue_2 = Vn + sqrt(beta*beta + Vn*Vn);
		//////Eigenvalue_3 = Vn - sqrt(beta*beta + Vn*Vn);
		//////Eigenvalue_4 = Vn;

		/////Eigenvectors
		//Eigenvector_11 = 0;
		//Eigenvector_21 = -ny / (nx + eps);
		//Eigenvector_31 = 1;
		//Eigenvector_41 = 0;

		//Eigenvector_13 = -(sqrt(beta *beta + Vn*Vn)*beta*beta * (-alpha_l_mean)) / (beta *beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) + v_mean*Vn);
		//Eigenvector_12 = -(sqrt(beta *beta + Vn*Vn)*beta*beta * (-rho_mean)) / (beta *beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) + v_mean*Vn);
		//Eigenvector_22 = (sqrt(beta * beta + Vn * Vn)*nx - Vn * nx + u_mean) / (sqrt(beta * beta + Vn * Vn)*ny - Vn * ny + v_mean);
		//Eigenvector_32 = 1;
		//Eigenvector_42 = 0;

		//Eigenvector_14 = -(sqrt(beta * beta + Vn * Vn)*beta * beta * (-alpha_l_mean)) / (-beta * beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) - v_mean*Vn);
		//Eigenvector_13 = -(sqrt(beta * beta + Vn * Vn)*beta * beta * (-rho_mean)) / (-beta * beta * ny + v_mean*sqrt(beta * beta + Vn * Vn) - v_mean*Vn);
		//Eigenvector_23 = (sqrt(beta * beta + Vn * Vn)*nx + Vn * nx - u_mean) / (sqrt(beta * beta + Vn * Vn)*ny + Vn * ny - v_mean);
		//Eigenvector_33 = 1;
		//Eigenvector_43 = 0;

		//Eigenvector_14 = 0;
		//Eigenvector_24 = 0;
		//Eigenvector_34 = 0;
		//Eigenvector_44 = 1;

		////Eigenvector_11 = 0;
		////Eigenvector_21 = (-2 * beta*beta*ny) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_31 = (2 * beta*beta*nx) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_41 = 0;

		////Eigenvector_12 = (C_bar*beta*beta*rho_mean) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_22 = (u_mean*(Vn + C_bar) + beta*beta*nx) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_32 = (v_mean*(Vn + C_bar) + beta*beta*ny) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_42 = 0;

		////Eigenvector_13 = (-C_bar*beta*beta*rho_mean) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_23 = (u_mean*(Vn - C_bar) + beta*beta*nx) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_33 = (v_mean*(Vn - C_bar) + beta*beta*ny) / (2 * beta*beta*C_bar*C_bar);
		////Eigenvector_43 = 0;

		////Eigenvector_14 = 0;
		////Eigenvector_24 = 0;
		////Eigenvector_34 = 0;
		////Eigenvector_44 = 1;


		////primitive values differences
		//dQ_1 = (pR - pL);
		//dQ_2 = (uR - uL);
		//dQ_3 = (vR - vL);
		//dQ_4 = (alpha_lR - alpha_lL);

		////Fluxes
		//FluxVectors(pL, uL, vL, rhoL, alpha_lL,  //inputs
		//	F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L);

		//FluxVectors(pR, uR, vR, rhoR, alpha_lR,  //inputs
		//	F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R);

		//fL_1 = F1L*nx + G1L*ny;
		//fL_2 = F2L*nx + G2L*ny;
		//fL_3 = F3L*nx + G3L*ny;
		//fL_4 = F4L*nx + G4L*ny;

		//fR_1 = F1R*nx + G1R*ny;
		//fR_2 = F2R*nx + G2R*ny;
		//fR_3 = F3R*nx + G3R*ny;
		//fR_4 = F4R*nx + G4R*ny;

		//double fp_1, fp_2, fp_3, fp_4;
		//fp_1 = Uc * 1;
		//fp_2 = Uc * rho_mean*u_mean + p_mean*nx;
		//fp_3 = Uc * rho_mean*v_mean + p_mean*ny;
		//fp_4 = Uc * alpha_l_mean;

		////Calculating WaveStrenghts
		////inverse of Eigenvectors

		//MatrixInversion4x4(Eigenvector_11, Eigenvector_12, Eigenvector_13, Eigenvector_14,
		//	Eigenvector_21, Eigenvector_22, Eigenvector_23, Eigenvector_24,
		//	Eigenvector_31, Eigenvector_32, Eigenvector_33, Eigenvector_34,
		//	Eigenvector_41, Eigenvector_42, Eigenvector_43, Eigenvector_44, //inputs
		//	inv_Eigenvector_11, inv_Eigenvector_12, inv_Eigenvector_13, inv_Eigenvector_14,
		//	inv_Eigenvector_21, inv_Eigenvector_22, inv_Eigenvector_23, inv_Eigenvector_24,
		//	inv_Eigenvector_31, inv_Eigenvector_32, inv_Eigenvector_33, inv_Eigenvector_34,
		//	inv_Eigenvector_41, inv_Eigenvector_42, inv_Eigenvector_43, inv_Eigenvector_44);



		/////////inverse of Eigenvectors
		////inv_Eigenvector_11 = (u_mean*ny - v_mean*nx) / (rho_mean);
		////inv_Eigenvector_21 = (C_bar - Vn) / (rho_mean);
		////inv_Eigenvector_31 = (-C_bar - Vn) / (rho_mean);
		////inv_Eigenvector_41 = 0;

		////inv_Eigenvector_12 = -v_mean*Vn - beta*beta*ny;
		////inv_Eigenvector_22 = beta*beta*nx;
		////inv_Eigenvector_32 = beta*beta*nx;
		////inv_Eigenvector_42 = 0;

		////inv_Eigenvector_13 = u_mean*Vn + beta*beta*nx;
		////inv_Eigenvector_23 = beta*beta*ny;
		////inv_Eigenvector_33 = beta*beta*ny;
		////inv_Eigenvector_43 = 0;

		////inv_Eigenvector_14 = 0;
		////inv_Eigenvector_24 = 0;
		////inv_Eigenvector_34 = 0;
		////inv_Eigenvector_44 = 2 * beta*beta*C_bar*C_bar;



		//WaveStrenghts_1 = inv_Eigenvector_11*dQ_1 + inv_Eigenvector_12*dQ_2 + inv_Eigenvector_13*dQ_3 + inv_Eigenvector_14*dQ_4;
		//WaveStrenghts_2 = inv_Eigenvector_21*dQ_1 + inv_Eigenvector_22*dQ_2 + inv_Eigenvector_23*dQ_3 + inv_Eigenvector_24*dQ_4;
		//WaveStrenghts_3 = inv_Eigenvector_31*dQ_1 + inv_Eigenvector_32*dQ_2 + inv_Eigenvector_33*dQ_3 + inv_Eigenvector_34*dQ_4;
		//WaveStrenghts_4 = inv_Eigenvector_41*dQ_1 + inv_Eigenvector_42*dQ_2 + inv_Eigenvector_43*dQ_3 + inv_Eigenvector_44*dQ_4;

		////Preconditioning matrix
		//double P11, P12, P13, P14, P21, P22, P23, P24, P31, P32, P33, P34, P41, P42, P43, P44, D1, D2, D3, D4;
		////PreconditioningMatrix(uL, vL, rhoL, alpha_lL, beta, 1-0.01,
		////	P11, P12, P13, P14,
		////	P21, P22, P23, P24,
		////	P31, P32, P33, P34,
		////	P41, P42, P43, P44);
		////PreconditioningMatrix(uR, vR, rhoR, alpha_lR, beta, 1 - 0.01,
		////	P11, P12, P13, P14,
		////	P21, P22, P23, P24,
		////	P31, P32, P33, P34,
		////	P41, P42, P43, P44);
		//PreconditioningMatrix(u_mean, v_mean, rho_mean, alpha_l_mean, beta, 1 - 0.01,
		//	P11, P12, P13, P14,
		//	P21, P22, P23, P24,
		//	P31, P32, P33, P34,
		//	P41, P42, P43, P44);

		////Calculating Roe's dissipation and numerical fluxes
		//D1_1 = Eigenvector_11 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_12 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_13 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_14 * abs(Eigenvalue_4)*WaveStrenghts_4;
		//D2_1 = Eigenvector_21 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_22 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_23 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_24 * abs(Eigenvalue_4)*WaveStrenghts_4;
		//D3_1 = Eigenvector_31 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_32 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_33 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_34 * abs(Eigenvalue_4)*WaveStrenghts_4;
		//D4_1 = Eigenvector_41 * abs(Eigenvalue_1)*WaveStrenghts_1 + Eigenvector_42 * abs(Eigenvalue_2)*WaveStrenghts_2 + Eigenvector_43 * abs(Eigenvalue_3)*WaveStrenghts_3 + Eigenvector_44 * abs(Eigenvalue_4)*WaveStrenghts_4;

		//D1 = D1_1*P11 + D2_1*P12 + D3_1*P13 + D4_1*P14;
		//D2 = D1_1*P21 + D2_1*P22 + D3_1*P23 + D4_1*P24;
		//D3 = D1_1*P31 + D2_1*P32 + D3_1*P33 + D4_1*P34;
		//D4 = D1_1*P41 + D2_1*P42 + D3_1*P43 + D4_1*P44;

		//num_flux_1 = 0.5 * f_M * (fL_1 + fR_1) + (1 - f_M)*fp_1 - 0.5 * D1;
		//num_flux_2 = 0.5 * f_M * (fL_2 + fR_2) + (1 - f_M)*fp_2 - 0.5 * D2;
		//num_flux_3 = 0.5 * f_M * (fL_3 + fR_3) + (1 - f_M)*fp_3 - 0.5 * D3;
		//num_flux_4 = 0.5 * f_M * (fL_4 + fR_4) + (1 - f_M)*fp_4 - 0.5 * D4;




	}
	else if (num_flag == 2)
	{
		//******************************//
		//*** Lax-Friedrichs/Rusonov ***//
		//******************************//

		double F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L;
		double F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R;
		double Q1L, Q2L, Q3L, Q4L, Q1R, Q2R, Q3R, Q4R;
		double fL1, fL2, fL3, fL4, fR1, fR2, fR3, fR4;
		double dQ1, dQ2, dQ3, dQ4, lambda;
		double WS1L, WS2L, WS3L, WS4L, WS1R, WS2R, WS3R, WS4R;

		///////////////////////

		FluxVectors(pL, uL, vL, rhoL, alpha_lL,  //inputs
			F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L);

		FluxVectors(pR, uR, vR, rhoR, alpha_lR,  //inputs
			F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R);


		SolutionVector(pL, uL, vL, alpha_lL,  //inputs
			Q1L, Q2L, Q3L, Q4L);

		SolutionVector(pR, uR, vR, alpha_lR,  //inputs
			Q1R, Q2R, Q3R, Q4R);


		dQ1 = Q1L - Q1R;
		dQ2 = Q2L - Q2R;
		dQ3 = Q3L - Q3R;
		dQ4 = Q4L - Q4R;

		////wave speed #1
		//WaveSpeed(uL, vL, nx, ny, beta, WS1L, WS2L, WS3L);
		//WaveSpeed(uR, vR, nx, ny, beta, WS1R, WS2R, WS3R);

		//lambda = maxval(abs(WS1L), abs(WS1R));
		//lambda = maxval(lambda, abs(WS2L));
		//lambda = maxval(lambda, abs(WS2R));
		//lambda = maxval(lambda, abs(WS3L));
		//lambda = maxval(lambda, abs(WS3R));
		  
		double coef = 0.2; 
		WS1L = sqrt(uL*uL + vL*vL) + sqrt(uL*uL + vL*vL + beta*beta);
		WS1R = sqrt(uR*uR + vR*vR) + sqrt(uR*uR + vR*vR + beta*beta);
		lambda = maxval(abs(WS1L), abs(WS1R));

		////wave speed #2
		//double VnL = uL*nx + vL*ny;
		//double VnR = uR*nx + vR*ny;
		//WS1L = abs(rhoL*VnL);
		//WS2L = abs(rhoL*VnL + sqrt(VnL*VnL*rhoL*rhoL + 1));
		//WS3L = abs(rhoL*VnL - sqrt(VnL*VnL*rhoL*rhoL + 1));
		//WS4L = VnL;

		//WS1R = abs(rhoR*VnR);
		//WS2R = abs(rhoR*VnR + sqrt(VnR*VnR*rhoR*rhoR + 1));
		//WS3R = abs(rhoR*VnR - sqrt(VnR*VnR*rhoR*rhoR + 1));
		//WS4R = VnR;

		//lambda = maxval(WS1L, WS1R);
		//lambda = maxval(lambda, WS2L);
		//lambda = maxval(lambda, WS2R);
		//lambda = maxval(lambda, WS3L);
		//lambda = maxval(lambda, WS3R);
		//lambda = maxval(lambda, WS4L);
		//lambda = maxval(lambda, WS4R);


		fL1 = F1L*nx + G1L*ny;
		fL2 = F2L*nx + G2L*ny;
		fL3 = F3L*nx + G3L*ny;
		fL4 = F4L*nx + G4L*ny;

		fR1 = F1R*nx + G1R*ny;
		fR2 = F2R*nx + G2R*ny;
		fR3 = F3R*nx + G3R*ny;
		fR4 = F4R*nx + G4R*ny;

		//num_flux_1 = 0.5 * (fL1 + fR1) + 0.5*lambda*(dQ1);
		//num_flux_2 = 0.5 * (fL2 + fR2) + 0.5*lambda*(dQ2);
		//num_flux_3 = 0.5 * (fL3 + fR3) + 0.5*lambda*(dQ3);
		//num_flux_4 = 0.5 * (fL4 + fR4) + 0.5*lambda*(dQ4);

		double P11, P12, P13, P14, P21, P22, P23, P24, P31, P32, P33, P34, P41, P42, P43, P44, dQ1_, dQ2_, dQ3_, dQ4_;
		PreconditioningMatrix(uL, vL, rhoL, alpha_lL, beta, 1 - 0.01,
			P11, P12, P13, P14,
			P21, P22, P23, P24,
			P31, P32, P33, P34,
			P41, P42, P43, P44);

		dQ1_ = dQ1*P11 + dQ2*P12 + dQ3*P13 + dQ4*P14;
		dQ2_ = dQ1*P21 + dQ2*P22 + dQ3*P23 + dQ4*P24;
		dQ3_ = dQ1*P31 + dQ2*P32 + dQ3*P33 + dQ4*P34;
		dQ4_ = dQ1*P41 + dQ2*P42 + dQ3*P43 + dQ4*P44;

		num_flux_1 = 0.5 * (fL1 + fR1) + 0.5*coef*lambda*(dQ1_);
		num_flux_2 = 0.5 * (fL2 + fR2) + 0.5*coef*lambda*(dQ2_);
		num_flux_3 = 0.5 * (fL3 + fR3) + 0.5*coef*lambda*(dQ3_);
		num_flux_4 = 0.5 * (fL4 + fR4) + 0.5*coef*lambda*(dQ4_);

	}
	else if (num_flag == 3)
	{
		//	//*************//
		//	//*** AUSM Family ***//
		//	//*************//


		////Schemes
		//double beta, alpha, c_bar, VnL, VnR, ML, MR, M_bar, p_bar;
		//double rhosai1, rhosai2, rhosai3, rhosai4, temp;
		//double aaR, aaL, rrhoR, rrhoL, ppR, ppL, HHR, HHL, uuR, uuL, vvR, vvL;

		////aaR=aL;
		////aaL=aR;
		////rrhoR=rhoL;
		////rrhoL=rhoR;
		////ppR=pL;
		////ppL=pR;
		////HHR=HL;
		////HHL=HR;
		////uuR=uL;
		////uuL=uR;
		////vvR=vL;
		////vvL=vR;

		//aaR=aR;
		//aaL=aL;
		//rrhoR=rhoL;
		//rrhoL=rhoR;
		//ppR=pL;
		//ppL=pR;
		//HHR=HL;
		//HHL=HL;
		//uuR=uR;
		//uuL=uL;
		//vvR=vL;
		//vvL=vR;

		////AUSM
		//beta=0;
		//alpha=0;
		////AUSM +
		//beta=1.0/8.0;
		//alpha=3.0/16.0;

		//c_bar=0.5*(aaL+aaR);
		//VnL = uuL*nx + vvL*ny;
		//VnR = uuR*nx + vvR*ny;
		//ML=VnL/c_bar;
		//MR=VnR/c_bar;
		//M_bar=M_plus_4(ML,beta)+M_minus_4(MR,beta);
		//p_bar=ppL*p_plus_5(ML,alpha)+ppR*p_minus_5(MR,alpha);
		//
		//if(M_bar>=0)
		//{
		//	rhosai1=rrhoL*1;
		//	rhosai2=rrhoL*uuL;
		//	rhosai3=rrhoL*vvL;
		//	rhosai4=rrhoL*HHL;
		//}
		//else
		//{
		//	rhosai1=rrhoR*1;
		//	rhosai2=rrhoR*uuR;
		//	rhosai3=rrhoR*vvR;
		//	rhosai4=rrhoR*HHR;
		//}

		//num_flux_1=(rhoL*uL)*nx+(rhoL*vL)*ny;
		//num_flux_2=(rhoL*uL*uL+pL)*nx+(rhoL*uL*vL)*ny;
		//num_flux_3=(rhoL*uL*vL)*nx+(rhoL*vL*vL+pL)*ny;
		//num_flux_5=(rhoL*uL*HL)*nx+(rhoL*vL*HL)*ny;

		//num_flux_1 = c_bar*M_bar*rhosai1 + p_bar*0;
		//num_flux_2 = c_bar*M_bar*rhosai2 + p_bar*nx;
		//num_flux_3 = c_bar*M_bar*rhosai3 + p_bar*ny;
		//num_flux_5 = c_bar*M_bar*rhosai4 + p_bar*0;


		////AUSM+-up Scheme
		//double beta, alpha, c_bar, VnL, VnR, ML, MR, M_bar, p_bar;
		//double rhosai1, rhosai2, rhosai3, rhosai4;
		//double sigma, Kp, M_co, rho_bar, M_bar_2, M0, fa, M_bar_p, Ku, p_bar_u;

		//c_bar=0.5*(aL+aR);
		//VnL = uL*nx + vL*ny;
		//VnR = uR*nx + vR*ny;
		//ML=VnL/c_bar;
		//MR=VnR/c_bar;

		//sigma=1;
		//Kp=0.25;
		//beta=1.0/8.0;
		//M_co=1e-4;
		//rho_bar=0.5*(rhoL+rhoR);
		//M_bar_2=0.5*(ML*ML+MR*MR);
		//M_bar=M_plus_4(ML,beta)+M_minus_4(MR,beta);

		//alpha=3.0/16.0;
		//Ku=0.75;
		//p_bar_u=-2*Ku*rho_bar*c_bar*c_bar*p_plus_5(ML,alpha)*p_minus_5(MR,alpha)*(MR-ML);
		//p_bar=pL*p_plus_5(ML,alpha)+pR*p_minus_5(MR,alpha)+p_bar_u;
		//
		//if(M_bar>=0)
		//{
		//	rhosai1=rhoL*1;
		//	rhosai2=rhoL*uL;
		//	rhosai3=rhoL*vL;
		//	rhosai4=rhoL*HL;
		//}
		//else
		//{
		//	rhosai1=rhoR*1;
		//	rhosai2=rhoR*uR;
		//	rhosai3=rhoR*vR;
		//	rhosai4=rhoR*HR;
		//}

		//num_flux_1 = c_bar*M_bar*rhosai1 + p_bar*0;
		//num_flux_2 = c_bar*M_bar*rhosai2 + p_bar*nx;
		//num_flux_3 = c_bar*M_bar*rhosai3 + p_bar*ny;
		//num_flux_5 = c_bar*M_bar*rhosai4 + p_bar*0;





		////*** AUSM+-up Scheme ***//
		//double beta, alpha, c_bar, VnL, VnR, ML, MR, M_bar, p_bar;
		//double rhosai1, rhosai2, rhosai3, rhosai4;
		//double sigma, Kp, M_co, rho_bar, M_bar_2, M0, fa, M_bar_p, Ku, p_bar_u, eL, eR;

		//
		//beta = 1.0 / 8.0; //1.0/8.0
		//alpha = 3.0 / 16.0; //3.0/16.0

		//M_co = 0.01; //1e-4
		//Ku = 0.75; //0.75
		//sigma = 1; //1
		//Kp = 0.25; //0.25

		//VnL = uL*nx + vL*ny;
		//VnR = uR*nx + vR*ny;

		////double aL = VnL + sqrt(VnL*VnL + beta*beta);
		////double aR = VnR + sqrt(VnR*VnR + beta*beta);

		//double aL = sqrt(uL*uL + vL*vL) + sqrt(uL*uL + vL*vL + beta*beta);
		//double aR = sqrt(uR*uR + vR*vR) + sqrt(uR*uR + vR*vR + beta*beta);

		//c_bar = 0.5*(aL + aR);
		////c_bar = sqrt(aL*aR);
		//rho_bar = 0.5*(rhoL + rhoR);

		//ML = VnL / c_bar;
		//MR = VnR / c_bar;

		//M_bar_2 = 0.5*(ML*ML + MR*MR);
		//double Mp = -Kp*maxval(1 - sigma*M_bar_2, 0)*(pR - pL) / (rho_bar*c_bar*c_bar);
		//M_bar = M_plus_4(ML, beta) + M_minus_4(MR, beta) + Mp;
		////M_bar = M_plus_4(ML, beta) + M_minus_4(MR, beta);

		//p_bar_u = -2 * Ku*rho_bar*c_bar*c_bar*p_plus_5(ML, alpha)*p_minus_5(MR, alpha)*(MR - ML);
		//p_bar = pL*p_plus_5(ML, alpha) + pR*p_minus_5(MR, alpha) + p_bar_u;

		//if (M_bar >= 0)
		//{
		//	rhosai1 = 1;
		//	rhosai2 = rhoL*uL;
		//	rhosai3 = rhoL*vL;
		//	rhosai4 = alpha_lL;
		//}
		//else
		//{
		//	rhosai1 = 1;
		//	rhosai2 = rhoR*uR;
		//	rhosai3 = rhoR*vR;
		//	rhosai4 = alpha_lR;
		//}

		//double F_AUSM_1, F_AUSM_2, F_AUSM_3, F_AUSM_4, tx, ty;
		//F_AUSM_1 = c_bar*M_bar*rhosai1 + p_bar * 0;
		//F_AUSM_2 = c_bar*M_bar*rhosai2 + p_bar*nx;
		//F_AUSM_3 = c_bar*M_bar*rhosai3 + p_bar*ny;
		//F_AUSM_4 = c_bar*M_bar*rhosai4 + p_bar * 0;

		//num_flux_1 = F_AUSM_1;
		//num_flux_2 = F_AUSM_2;
		//num_flux_3 = F_AUSM_3;
		//num_flux_4 = F_AUSM_4;


		////tx = -ny;
		////ty = +nx;
		////num_flux_1 = F_AUSM_1;
		////num_flux_2 = nx*F_AUSM_2 + tx*F_AUSM_3;
		////num_flux_3 = ny*F_AUSM_2 + ty*F_AUSM_3;
		////num_flux_5 = F_AUSM_4;
		////num_flux_5 = F_AUSM_4;







		//AUSM+-up all speed Scheme
		double alpha, c_bar, VnL, VnR, ML, MR, M_bar, p_bar;
		double rhosai1, rhosai2, rhosai3, rhosai4;
		double sigma, Kp, M_co, rho_bar, M_bar_2, M0, fa, M_bar_p, Ku, p_bar_u, eL, eR;

		
		double beta_AUSM = 1.0 / 8.0; //1.0/8.0
		alpha = 3.0 / 16.0; //3.0/16.0

		M_co = 0.2; //1e-4
		Ku = 0.75; //0.75
		sigma = 1; //1
		Kp = 0.25; //0.25

		VnL = uL*nx + vL*ny;
		VnR = uR*nx + vR*ny;

		//double aL = VnL + sqrt(VnL*VnL + beta*beta);
		//double aR = VnR + sqrt(VnR*VnR + beta*beta);

		double aL = sqrt(uL*uL + vL*vL + beta*beta);
		double aR = sqrt(uR*uR + vR*vR + beta*beta);

		c_bar = 0.5*(aL + aR);
		rho_bar = 0.5*(rhoL + rhoR);

		ML = VnL / c_bar;
		MR = VnR / c_bar;

		M_bar_2 = 0.5*(ML*ML + MR*MR);
		M0 = minval(1, maxval(M_bar_2, M_co*M_co));
		fa = M0*(2 - M0);
		M_bar_p = -Kp / fa*maxval(1 - sigma*M_bar_2*M_bar_2, 0)*(pR - pL) / rho_bar / c_bar / c_bar;
		M_bar = M_plus_4(ML, beta_AUSM) + M_minus_4(MR, beta_AUSM) + M_bar_p;

		alpha = alpha*(-4 + 5 * fa*fa);
		p_bar_u = -2 * fa*Ku*rho_bar*c_bar*c_bar*p_plus_5(ML, alpha)*p_minus_5(MR, alpha)*(MR - ML);
		p_bar = pL*p_plus_5(ML, alpha) + pR*p_minus_5(MR, alpha) + p_bar_u;


		if (M_bar >= 0)
		{
			rhosai1 = 1;
			rhosai2 = rhoL*uL;
			rhosai3 = rhoL*vL;
			rhosai4 = alpha_lL;
		}
		else
		{
			rhosai1 = 1;
			rhosai2 = rhoR*uR;
			rhosai3 = rhoR*vR;
			rhosai4 = alpha_lR;
		}

		double F_AUSM_1, F_AUSM_2, F_AUSM_3, F_AUSM_4, tx, ty;
		F_AUSM_1 = c_bar*M_bar*rhosai1 + p_bar * 0;
		F_AUSM_2 = c_bar*M_bar*rhosai2 + p_bar*nx;
		F_AUSM_3 = c_bar*M_bar*rhosai3 + p_bar*ny;
		F_AUSM_4 = c_bar*M_bar*rhosai4 + p_bar * 0;

		num_flux_1 = F_AUSM_1;
		num_flux_2 = F_AUSM_2;
		num_flux_3 = F_AUSM_3;
		num_flux_4 = F_AUSM_4;
	}
	else if (num_flag == 4)
	{
		////double tx, ty;
		////double UL_1, UL_2, UL_3, UL_4, UR_1, UR_2, UR_3, UR_4;
		////double U_hat_L_1, U_hat_L_2, U_hat_L_3, U_hat_L_4, U_hat_R_1, U_hat_R_2, U_hat_R_3, U_hat_R_4;
		////double rho_hat_L, u_hat_L, v_hat_L, E_hat_L, p_hat_L, H_hat_L, rho_hat_R, u_hat_R, v_hat_R, E_hat_R, p_hat_R, H_hat_R;
		////double Flux_hat_L_1, Flux_hat_L_2, Flux_hat_L_3, Flux_hat_L_4, Flux_hat_R_1, Flux_hat_R_2, Flux_hat_R_3, Flux_hat_R_4;
		////double SL, SR;
		////double F_HLL_1, F_HLL_2, F_HLL_3, F_HLL_4;

		////tx=-ny;
		////ty=+nx;

		////UL_1=rhoL;
		////UL_2=rhoL*uL;
		////UL_3=rhoL*vL;
		////UL_4=EL;

		////UR_1=rhoR;
		////UR_2=rhoR*uR;
		////UR_3=rhoR*vR;
		////UR_4=ER;

		//////  [1   0   0  0]
		//////T=[0  nx  ny  0]
		//////  [0  tx  ty  0]
		//////  [0   0   0  1]
		//////       [1   0   0  0]
		//////inv(T)=[0  nx  tx  0]
		//////       [0  ny  ty  0]
		//////       [0   0   0  1]

		//////U_hat=T*U
		////U_hat_L_1=UL_1;
		////U_hat_L_2=nx*UL_2+ny*UL_3;
		////U_hat_L_3=tx*UL_2+ty*UL_3;
		////U_hat_L_4=UL_4;

		////U_hat_R_1=UR_1;
		////U_hat_R_2=nx*UR_2+ny*UR_3;
		////U_hat_R_3=tx*UR_2+ty*UR_3;
		////U_hat_R_4=UR_4;

		////rho_hat_L=U_hat_L_1;
		////u_hat_L=U_hat_L_2/U_hat_L_1;
		////v_hat_L=U_hat_L_3/U_hat_L_1;
		////E_hat_L=U_hat_L_4;
		////p_hat_L=pL;
		////H_hat_L=(E_hat_L+p_hat_L)/rho_hat_L;

		////rho_hat_R=U_hat_R_1;
		////u_hat_R=U_hat_R_2/U_hat_R_1;
		////v_hat_R=U_hat_R_3/U_hat_R_1;
		////E_hat_R=U_hat_R_4;
		////p_hat_R=pR;
		////H_hat_R=(E_hat_R+p_hat_R)/rho_hat_R;


		////Flux_hat_L_1=rho_hat_L*u_hat_L;
		////Flux_hat_L_2=rho_hat_L*u_hat_L*u_hat_L+p_hat_L;
		////Flux_hat_L_3=rho_hat_L*u_hat_L*v_hat_L;
		////Flux_hat_L_4=u_hat_L*(E_hat_L+p_hat_L);

		////Flux_hat_R_1=rho_hat_R*u_hat_R;
		////Flux_hat_R_2=rho_hat_R*u_hat_R*u_hat_R+p_hat_R;
		////Flux_hat_R_3=rho_hat_R*u_hat_R*v_hat_R;
		////Flux_hat_R_4=u_hat_R*(E_hat_R+p_hat_R);

		//////Always rho_hat=rho & E_hat=E & since u_mean^2+v_mean^2=u_hat^2+v_hat^2 --> p_hat=p and therefore a_hat=a
		////SL=minval(u_hat_L-aL, u_hat_R-aR);
		////SR=maxval(u_hat_L+aL, u_hat_R+aR);
		//////S_star=( p_hat_R-p_hat_L+rho_hat_L*u_hat_L*(SL-u_hat_L)-rho_hat_R*u_hat_R*(SR-u_hat_R) )/( rho_hat_L*(SL-u_hat_L)-rho_hat_R*(SR-u_hat_R) );

		////if (SL>=0)
		////{
		////	F_HLL_1=Flux_hat_L_1;
		////	F_HLL_2=Flux_hat_L_2;
		////	F_HLL_3=Flux_hat_L_3;
		////	F_HLL_4=Flux_hat_L_4;
		////}
		////else if (SR<=0)
		////{
		////	F_HLL_1=Flux_hat_R_1;
		////	F_HLL_2=Flux_hat_R_2;
		////	F_HLL_3=Flux_hat_R_3;
		////	F_HLL_4=Flux_hat_R_4;
		////}
		////else
		////{
		////	F_HLL_1=( SR*Flux_hat_L_1-SL*Flux_hat_R_1+SL*SR*(U_hat_R_1-U_hat_L_1) )/( SR-SL );
		////	F_HLL_2=( SR*Flux_hat_L_2-SL*Flux_hat_R_2+SL*SR*(U_hat_R_2-U_hat_L_2) )/( SR-SL );
		////	F_HLL_3=( SR*Flux_hat_L_3-SL*Flux_hat_R_3+SL*SR*(U_hat_R_3-U_hat_L_3) )/( SR-SL );
		////	F_HLL_4=( SR*Flux_hat_L_4-SL*Flux_hat_R_4+SL*SR*(U_hat_R_4-U_hat_L_4) )/( SR-SL );
		////}
		////num_flux_1 = F_HLL_1;
		////num_flux_2 = nx*F_HLL_2 + tx*F_HLL_3;
		////num_flux_3 = ny*F_HLL_2 + ty*F_HLL_3;
		////num_flux_5 = F_HLL_4;
		////num_flux_5 = F_HLL_4;




		//double tx, ty;
		//double UL_1, UL_2, UL_3, UL_4, UR_1, UR_2, UR_3, UR_4;
		//double U_hat_L_1, U_hat_L_2, U_hat_L_3, U_hat_L_4, U_hat_R_1, U_hat_R_2, U_hat_R_3, U_hat_R_4;
		//double rho_hat_L, u_hat_L, v_hat_L, E_hat_L, p_hat_L, H_hat_L, rho_hat_R, u_hat_R, v_hat_R, E_hat_R, p_hat_R, H_hat_R;
		//double Flux_hat_L_1, Flux_hat_L_2, Flux_hat_L_3, Flux_hat_L_4, Flux_hat_R_1, Flux_hat_R_2, Flux_hat_R_3, Flux_hat_R_4;
		//double SL, SR;
		//double F_HLL_1, F_HLL_2, F_HLL_3, F_HLL_4;

		//tx = -ny;
		//ty = +nx;

		//UL_1 = rhoL;
		//UL_2 = rhoL*uL;
		//UL_3 = rhoL*vL;
		//UL_4 = EL;

		//UR_1 = rhoR;
		//UR_2 = rhoR*uR;
		//UR_3 = rhoR*vR;
		//UR_4 = ER;

		////  [1   0   0  0]
		////T=[0  nx  ny  0]
		////  [0  tx  ty  0]
		////  [0   0   0  1]
		////       [1   0   0  0]
		////inv(T)=[0  nx  tx  0]
		////       [0  ny  ty  0]
		////       [0   0   0  1]

		////U_hat=T*U
		//U_hat_L_1 = UL_1;
		//U_hat_L_2 = nx*UL_2 + ny*UL_3;
		//U_hat_L_3 = tx*UL_2 + ty*UL_3;
		//U_hat_L_4 = UL_4;

		//U_hat_R_1 = UR_1;
		//U_hat_R_2 = nx*UR_2 + ny*UR_3;
		//U_hat_R_3 = tx*UR_2 + ty*UR_3;
		//U_hat_R_4 = UR_4;

		//rho_hat_L = U_hat_L_1;
		//u_hat_L = U_hat_L_2 / U_hat_L_1;
		//v_hat_L = U_hat_L_3 / U_hat_L_1;
		//E_hat_L = U_hat_L_4;
		//p_hat_L = pL;
		//H_hat_L = (E_hat_L + p_hat_L) / rho_hat_L;

		//rho_hat_R = U_hat_R_1;
		//u_hat_R = U_hat_R_2 / U_hat_R_1;
		//v_hat_R = U_hat_R_3 / U_hat_R_1;
		//E_hat_R = U_hat_R_4;
		//p_hat_R = pR;
		//H_hat_R = (E_hat_R + p_hat_R) / rho_hat_R;


		//Flux_hat_L_1 = rho_hat_L*u_hat_L;
		//Flux_hat_L_2 = rho_hat_L*u_hat_L*u_hat_L + p_hat_L;
		//Flux_hat_L_3 = rho_hat_L*u_hat_L*v_hat_L;
		//Flux_hat_L_4 = u_hat_L*(E_hat_L + p_hat_L);

		//Flux_hat_R_1 = rho_hat_R*u_hat_R;
		//Flux_hat_R_2 = rho_hat_R*u_hat_R*u_hat_R + p_hat_R;
		//Flux_hat_R_3 = rho_hat_R*u_hat_R*v_hat_R;
		//Flux_hat_R_4 = u_hat_R*(E_hat_R + p_hat_R);

		//double rho, u_mean, v_mean, H, c2, c, t1, t2, t3;
		//rho=rho_hat_L*rho_hat_R;
		//u_mean=(rho_hat_L*u_hat_L+rho_hat_R*u_hat_R)/(rho_hat_L+rho_hat_R);
		//v_mean=(rho_hat_L*v_hat_L+rho_hat_R*v_hat_R)/(rho_hat_L+rho_hat_R);
		//H=(rho_hat_L*H_hat_L+rho_hat_R*H_hat_R)/(rho_hat_L+rho_hat_R);
		//c=(rho_hat_L*aL+rho_hat_R*aR)/(rho_hat_L+rho_hat_R);
		//c2=c*c;

		//SL=minval(u_hat_L-aL, u_mean-c);
		//SR=maxval(u_hat_R+aR, u_mean+c);
		//t1 = (minval(SR,0)-minval(0,SL))/(SR-SL);
		//t2 = 1-t1;
		//t3 = (SR*abs(SL)-SL*abs(SR))/(2*(SR-SL));

		//F_HLL_1 = t1*Flux_hat_R_1 + t2*Flux_hat_L_1 - t3*(U_hat_R_1-U_hat_L_1);
		//F_HLL_2 = t1*Flux_hat_R_2 + t2*Flux_hat_L_2 - t3*(U_hat_R_2-U_hat_L_2);
		//F_HLL_3 = t1*Flux_hat_R_3 + t2*Flux_hat_L_3 - t3*(U_hat_R_3-U_hat_L_3);
		//F_HLL_4 = t1*Flux_hat_R_4 + t2*Flux_hat_L_4 - t3*(U_hat_R_4-U_hat_L_4);

		//num_flux_1 = F_HLL_1;
		//num_flux_2 = nx*F_HLL_2+tx*F_HLL_3;
		//num_flux_3 = ny*F_HLL_2+ty*F_HLL_3;
		//num_flux_5 = F_HLL_4;
		//num_flux_5 = F_HLL_4;


	}




	else if (num_flag == 5)
	{
		//double tx, ty;
		////double UL_1, UL_2, UL_3, UL_4, UR_1, UR_2, UR_3, UR_4;
		////double U_hat_L_1, U_hat_L_2, U_hat_L_3, U_hat_L_4, U_hat_R_1, U_hat_R_2, U_hat_R_3, U_hat_R_4;
		////double rho_hat_L, u_hat_L, v_hat_L, E_hat_L, p_hat_L, H_hat_L, rho_hat_R, u_hat_R, v_hat_R, E_hat_R, p_hat_R, H_hat_R;
		////double Flux_hat_L_1, Flux_hat_L_2, Flux_hat_L_3, Flux_hat_L_4, Flux_hat_R_1, Flux_hat_R_2, Flux_hat_R_3, Flux_hat_R_4;
		////double SL, SR;
		////double F_HLLC_1, F_HLLC_2, F_HLLC_3, F_HLLC_4;
		////double S_star, p_star;
		////double U_hat_star_L_1, U_hat_star_L_2, U_hat_star_L_3, U_hat_star_L_4;
		////double U_hat_star_R_1, U_hat_star_R_2, U_hat_star_R_3, U_hat_star_R_4;

		////tx=-ny;
		////ty=+nx;

		////UL_1=rhoL;
		////UL_2=rhoL*uL;
		////UL_3=rhoL*vL;
		////UL_4=EL;

		////UR_1=rhoR;
		////UR_2=rhoR*uR;
		////UR_3=rhoR*vR;
		////UR_4=ER;

		//////  [1   0   0  0]
		//////T=[0  nx  ny  0]
		//////  [0  tx  ty  0]
		//////  [0   0   0  1]
		//////       [1   0   0  0]
		//////inv(T)=[0  nx  tx  0]
		//////       [0  ny  ty  0]
		//////       [0   0   0  1]

		//////U_hat=T*U
		////U_hat_L_1=UL_1;
		////U_hat_L_2=nx*UL_2+ny*UL_3;
		////U_hat_L_3=tx*UL_2+ty*UL_3;
		////U_hat_L_4=UL_4;

		////U_hat_R_1=UR_1;
		////U_hat_R_2=nx*UR_2+ny*UR_3;
		////U_hat_R_3=tx*UR_2+ty*UR_3;
		////U_hat_R_4=UR_4;

		////rho_hat_L=U_hat_L_1;
		////u_hat_L=U_hat_L_2/U_hat_L_1;
		////v_hat_L=U_hat_L_3/U_hat_L_1;
		////E_hat_L=U_hat_L_4;
		////p_hat_L=pL;
		////H_hat_L=(E_hat_L+p_hat_L)/rho_hat_L;

		////rho_hat_R=U_hat_R_1;
		////u_hat_R=U_hat_R_2/U_hat_R_1;
		////v_hat_R=U_hat_R_3/U_hat_R_1;
		////E_hat_R=U_hat_R_4;
		////p_hat_R=pR;
		////H_hat_R=(E_hat_R+p_hat_R)/rho_hat_R;


		////Flux_hat_L_1=rho_hat_L*u_hat_L;
		////Flux_hat_L_2=rho_hat_L*u_hat_L*u_hat_L+p_hat_L;
		////Flux_hat_L_3=rho_hat_L*u_hat_L*v_hat_L;
		////Flux_hat_L_4=u_hat_L*(E_hat_L+p_hat_L);

		////Flux_hat_R_1=rho_hat_R*u_hat_R;
		////Flux_hat_R_2=rho_hat_R*u_hat_R*u_hat_R+p_hat_R;
		////Flux_hat_R_3=rho_hat_R*u_hat_R*v_hat_R;
		////Flux_hat_R_4=u_hat_R*(E_hat_R+p_hat_R);

		//////Always rho_hat=rho & E_hat=E & since u_mean^2+v_mean^2=u_hat^2+v_hat^2 --> p_hat=p and therefore a_hat=a
		////SL=minval(u_hat_L-aL, u_hat_R-aR);
		////SR=maxval(u_hat_R+aR, u_hat_R+aR);
		////S_star=( p_hat_R-p_hat_L+rho_hat_L*u_hat_L*(SL-u_hat_L)-rho_hat_R*u_hat_R*(SR-u_hat_R) )/( rho_hat_L*(SL-u_hat_L)-rho_hat_R*(SR-u_hat_R) );
		////p_star=rho_hat_L*(u_hat_L-SL)*(u_hat_L-S_star)+pL;



		////U_hat_star_L_1=rho_hat_L*(SL-u_hat_L)/(SL-S_star)*(1)+(1/(SL-S_star))*(0);
		////U_hat_star_L_2=rho_hat_L*(SL-u_hat_L)/(SL-S_star)*(u_hat_L)+(1/(SL-S_star))*(p_star-pL);
		////U_hat_star_L_3=rho_hat_L*(SL-u_hat_L)/(SL-S_star)*(v_hat_L)+(1/(SL-S_star))*(0);
		////U_hat_star_L_4=rho_hat_L*(SL-u_hat_L)/(SL-S_star)*(E_hat_L/rho_hat_L)+(1/(SL-S_star))*(p_star*S_star-pL*u_hat_L);

		////U_hat_star_R_1=rho_hat_R*(SL-u_hat_R)/(SL-S_star)*(1)+(1/(SL-S_star))*(0);
		////U_hat_star_R_2=rho_hat_R*(SL-u_hat_R)/(SL-S_star)*(u_hat_R)+(1/(SL-S_star))*(p_star-pL);
		////U_hat_star_R_3=rho_hat_R*(SL-u_hat_R)/(SL-S_star)*(v_hat_R)+(1/(SL-S_star))*(0);
		////U_hat_star_R_4=rho_hat_R*(SL-u_hat_R)/(SL-S_star)*(E_hat_R/rho_hat_L)+(1/(SL-S_star))*(p_star*S_star-pL*u_hat_R);

		////if (SL>=0)
		////{
		////	F_HLLC_1=Flux_hat_L_1;
		////	F_HLLC_2=Flux_hat_L_2;
		////	F_HLLC_3=Flux_hat_L_3;
		////	F_HLLC_4=Flux_hat_L_4;
		////}
		////else if (SR<0)
		////{
		////	F_HLLC_1=Flux_hat_R_1;
		////	F_HLLC_2=Flux_hat_R_2;
		////	F_HLLC_3=Flux_hat_R_3;
		////	F_HLLC_4=Flux_hat_R_4;
		////}
		////else if (SL<0 && S_star>=0)
		////{
		////	F_HLLC_1=Flux_hat_L_1+SL*(U_hat_star_L_1-U_hat_L_1);
		////	F_HLLC_2=Flux_hat_L_2+SL*(U_hat_star_L_2-U_hat_L_2);
		////	F_HLLC_3=Flux_hat_L_3+SL*(U_hat_star_L_3-U_hat_L_3);
		////	F_HLLC_4=Flux_hat_L_4+SL*(U_hat_star_L_4-U_hat_L_4);
		////}
		////else
		////{
		////	F_HLLC_1=Flux_hat_R_1+SL*(U_hat_star_R_1-U_hat_R_1);
		////	F_HLLC_2=Flux_hat_R_2+SL*(U_hat_star_R_2-U_hat_R_2);
		////	F_HLLC_3=Flux_hat_R_3+SL*(U_hat_star_R_3-U_hat_R_3);
		////	F_HLLC_4=Flux_hat_R_4+SL*(U_hat_star_R_4-U_hat_R_4);
		////}

		////num_flux_1 = F_HLLC_1;
		////num_flux_2 = nx*F_HLLC_2+tx*F_HLLC_3;
		////num_flux_3 = ny*F_HLLC_2+ty*F_HLLC_3;
		////num_flux_5 = F_HLLC_4;
		////num_flux_5 = F_HLLC_4;
		//




		////From paper: ON THE CHOICE OF WAVESPEEDS FOR THE HLLC RIEMANN SOLVER
		//double eL, eR;
		//double qL, qR, cL, cR;
		//double Rrho, u_tilde, v_tilde, c_tilde, H_tilde, q_tilde, SL, SR, SM;
		//double rho_star_L, p_star, rho_u_star_L, rho_v_star_L, E_star_L, q_star_L;
		//double rho_star_R, rho_u_star_R, rho_v_star_R, E_star_R, q_star_R;
		//double U_star_L_1, U_star_L_2, U_star_L_3, U_star_L_4, U_star_R_1, U_star_R_2, U_star_R_3, U_star_R_4;
		//double F_L_1, F_L_2, F_L_3, F_L_4, F_R_1, F_R_2, F_R_3, F_R_4;
		//double F_star_L_1, F_star_L_2, F_star_L_3, F_star_L_4, F_star_R_1, F_star_R_2, F_star_R_3, F_star_R_4, e_star_L, e_star_R;

		//eL = (EL - 0.5*rhoL * (uL * uL + vL * vL)) / rhoL;
		//eR = (ER - 0.5*rhoR * (uR * uR + vR * vR)) / rhoR;

		//qL = uL*nx + vL*ny; //Eq. 19
		//qR = uR*nx + vR*ny;
		//cL = aL;
		//cR = aR;

		//Rrho = sqrt(rhoR / rhoL); //Eq. 51
		//u_tilde = (uL + uR*Rrho) / (1 + Rrho);
		//v_tilde = (vL + vR*Rrho) / (1 + Rrho);
		//c_tilde = (aL + aR*Rrho) / (1 + Rrho);
		//H_tilde = (HL + HR*Rrho) / (1 + Rrho);
		//q_tilde = u_tilde*nx + v_tilde*ny;
		//SL = minval(qL - cL, q_tilde - c_tilde);
		//SR = maxval(qR + cR, q_tilde + c_tilde);


		//SM = (rhoR*qR*(SR - qR) - rhoL*qL*(SL - qL) + pL - pR) / (rhoR*(SR - qR) - rhoL*(SL - qL)); //Eqs. 35 to 40


		//rho_star_L = rhoL*(SL - qL) / (SL - SM);
		//p_star=rhoL*(qL-SL)*(qL-SM)+pL;
		//rho_u_star_L = ((SL - qL)*rhoL*uL + (p_star - pL)*nx) / (SL - SM);
		//rho_v_star_L = ((SL - qL)*rhoL*vL + (p_star - pL)*ny) / (SL - SM);
		//E_star_L = ((SL - qL)*EL -pL*qL+p_star*SM) / (SL - SM);
		//q_star_L = SM; //Eq. 33
		//e_star_L = (E_star_L - 0.5*rho_star_L * ((rho_u_star_L / rho_star_L)* (rho_u_star_L / rho_star_L) + (rho_v_star_L / rho_star_L) * (rho_v_star_L / rho_star_L))) / rho_star_L;

		//rho_star_R = rhoR*(SR - qR) / (SR - SM);
		//p_star = rhoR*(qR - SR)*(qR - SM) + pR;
		//rho_u_star_R = ((SR - qR)*rhoR*uR + (p_star - pR)*nx) / (SR - SM);
		//rho_v_star_R = ((SR - qR)*rhoR*vR + (p_star - pR)*ny) / (SR - SM);
		//E_star_R = ((SR - qR)*ER - pR*qR + p_star*SM) / (SR - SM);
		//q_star_R = SM;
		//e_star_R = (E_star_R - 0.5*rho_star_R * ((rho_u_star_R / rho_star_R)* (rho_u_star_R / rho_star_R) + (rho_v_star_R / rho_star_R) * (rho_v_star_R / rho_star_R))) / rho_star_R;


		//U_star_L_1 = rho_star_L; //Eq. (32)
		//U_star_L_2 = rho_u_star_L;
		//U_star_L_3 = rho_v_star_L;
		////U_star_L_4 = E_star_L;
		//U_star_L_4 = beta*E_star_L + (1 - beta)*(rho_star_L*e_star_L);

		//U_star_R_1 = rho_star_R;
		//U_star_R_2 = rho_u_star_R;
		//U_star_R_3 = rho_v_star_R;
		////U_star_R_4 = E_star_R;
		//U_star_R_4 = beta*E_star_R + (1 - beta)*(rho_star_R*e_star_R);

		////F1L = rhoL * uL;
		////F2L = rhoL * uL * uL + pL;
		////F3L = rhoL * uL * vL;
		////F4L = rhoL * uL * eL + beta*uL * (0.5*rhoL * (uL * uL + vL * vL) + pL);
		////G1L = rhoL * vL;
		////G2L = rhoL * vL * uL;
		////G3L = rhoL * vL * vL + pL;
		////G4L = rhoL * vL * eL + beta*vL * (0.5*rhoL * (uL * uL + vL * vL) + pL);

		////F1R = rhoR * uR;
		////F2R = rhoR * uR * uR + pR;
		////F3R = rhoR * uR * vR;
		////F4R = rhoR * uR * eR + beta*uR * (0.5*rhoR * (uR * uR + vR * vR) + pR);
		////G1R = rhoR * vR;
		////G2R = rhoR * vR * uR;
		////G3R = rhoR * vR * vR + pR;
		////G4R = rhoR * vR * eR + beta*vR * (0.5*rhoR * (uR * uR + vR * vR) + pR);

		//F_L_1 = rhoL*qL; //Eq. (32)
		//F_L_2 = rhoL*uL*qL+pL*nx;
		//F_L_3 = rhoL*vL*qL + pL*ny;
		////F_L_4 = (EL+pL)*qL;
		//F_L_4 = rhoL * qL * eL + beta*qL * (0.5*rhoL * (uL * uL + vL * vL) + pL);

		//F_R_1 = rhoR*qR; //Eq. (32)
		//F_R_2 = rhoR*uR*qR + pR*nx;
		//F_R_3 = rhoR*vR*qR + pR*ny;
		////F_R_4 = (ER + pR)*qR;
		//F_R_4 = rhoR * qR * eR + beta*qR * (0.5*rhoR * (uR * uR + vR * vR) + pR);

		//F_star_L_1 = rho_star_L*q_star_L;
		//F_star_L_2 = rho_u_star_L*q_star_L + p_star*nx;
		//F_star_L_3 = rho_v_star_L*q_star_L + p_star*ny;
		////F_star_L_4 = (E_star_L + p_star)*q_star_L;
		//F_star_L_4 = rho_star_L * q_star_L * e_star_L + beta*q_star_L * (
		//	0.5*rho_star_L * ((rho_u_star_L / rho_star_L)* (rho_u_star_L / rho_star_L) + (rho_v_star_L / rho_star_L) * (rho_v_star_L / rho_star_L)) + p_star);


		//F_star_R_1 = rho_star_R*q_star_R;
		//F_star_R_2 = rho_u_star_R*q_star_R + p_star*nx;
		//F_star_R_3 = rho_v_star_R*q_star_R + p_star*ny;
		////F_star_R_4 = (E_star_R + p_star)*q_star_R;
		//F_star_R_4 = rho_star_R * q_star_R * e_star_R + beta*q_star_R * (
		//	0.5*rho_star_R * ((rho_u_star_R / rho_star_R)* (rho_u_star_R / rho_star_R) + (rho_v_star_R / rho_star_R) * (rho_v_star_R / rho_star_R)) + p_star);

		//if (SL>0) //Eq. 26
		//{
		//	num_flux_1 = F_L_1;
		//	num_flux_2 = F_L_2;
		//	num_flux_3 = F_L_3;
		//	num_flux_5 = F_L_4;
		//}
		//else if (SR<0)
		//{
		//	num_flux_1 = F_R_1;
		//	num_flux_2 = F_R_2;
		//	num_flux_3 = F_R_3;
		//	num_flux_5 = F_R_4;
		//}
		//else if (SL<=0 && SM>0)
		//{
		//	num_flux_1 = F_star_L_1;
		//	num_flux_2 = F_star_L_2;
		//	num_flux_3 = F_star_L_3;
		//	num_flux_5 = F_star_L_4;
		//}
		//else
		//{
		//	num_flux_1 = F_star_R_1;
		//	num_flux_2 = F_star_R_2;
		//	num_flux_3 = F_star_R_3;
		//	num_flux_5 = F_star_R_4;
		//}



	}



	else if (num_flag == 6)
	{
		//	//*************//
		//	//*** LDFSS ***//
		//	//*************//


		double VnL = uL*nx + vL*ny;
		double VnR = uR*nx + vR*ny;

		//double aL = VnL + sqrt(VnL*VnL + beta*beta);
		//double aR = VnR + sqrt(VnR*VnR + beta*beta);

		double aL = sqrt(uL*uL + vL*vL) + sqrt(uL*uL + vL*vL + beta*beta);
		double aR = sqrt(uR*uR + vR*vR) + sqrt(uR*uR + vR*vR + beta*beta);

		double a_mid = 0.5*(aL + aR);
		double V_mid = 0.5*(VnL + VnR);

		//double v_ref_mid = sqrt(minval(a_mid*a_mid, maxval(V_mid*V_mid, 15 * 15)));
		double v_ref_mid = a_mid;
		double M_ref_mid = v_ref_mid / a_mid;

		double a_tilde = sqrt(((1 - M_ref_mid*M_ref_mid)*(1 - M_ref_mid*M_ref_mid)*V_mid*V_mid + 4 * v_ref_mid*v_ref_mid)) / (1 + M_ref_mid*M_ref_mid);

		double ML = VnL / a_tilde;
		double MR = VnR / a_tilde;

		double alpha_plus_L = 0.5*(1 + sign(ML));
		double alpha_minus_L = 0.5*(1 - sign(ML));

		double alpha_plus_R = 0.5*(1 + sign(MR));
		double alpha_minus_R = 0.5*(1 - sign(MR));

		double beta_L = -maxval(0, 1 - floor(abs(ML)));
		double beta_R = -maxval(0, 1 - floor(abs(MR)));

		double M_plus = alpha_plus_L*(1 + beta_L)*ML - beta_L*M_plus_2(ML);
		double M_minus = alpha_minus_R*(1 + beta_R)*MR - beta_R*M_minus_2(MR);

		double M_mid = 0.5*(M_plus - alpha_plus_L*ML - M_minus + alpha_minus_R*MR);

		double rho_mid = 0.5*(rhoL + rhoR);

		double p_plus = alpha_plus_L*(1 + beta_L) - beta_L*0.5*(1 + ML);
		double p_minus = alpha_minus_R*(1 + beta_R) - beta_R*0.5*(1 - MR);

		double U_plus = a_tilde*(M_plus - M_mid*(1 - (pL - pR) / (2 * rhoL*v_ref_mid*v_ref_mid)));
		double U_minus = a_tilde*(M_minus + M_mid*(1 + (pL - pR) / (2 * rhoR*v_ref_mid*v_ref_mid)));

		double p_mid = 0.5*(pL + pR) + 0.5*(p_plus - p_minus)*(pL - pR) + rho_mid*v_ref_mid*v_ref_mid*(p_plus + p_minus - 1);


		num_flux_1 = U_plus*(1) + U_minus*(1) + p_mid*(0);
		num_flux_2 = U_plus*(rhoL*uL) + U_minus*(rhoR*uR) + p_mid*(nx);
		num_flux_3 = U_plus*(rhoL*vL) + U_minus*(rhoR*vR) + p_mid*(ny);
		num_flux_4 = U_plus*(alpha_lL)+U_minus*(alpha_lR)+p_mid*(0);

	}
	else if (num_flag == 7)
	{
		//******************************//
		//*** Lax-Friedrichs/Rusonov ***//
		//******************************//

		double F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L;
		double F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R;
		double fL1, fL2, fL3, fL4, fR1, fR2, fR3, fR4;

		///////////////////////

		FluxVectors(pL, uL, vL, rhoL, alpha_lL,  //inputs
			F1L, F2L, F3L, F4L, G1L, G2L, G3L, G4L);

		FluxVectors(pR, uR, vR, rhoR, alpha_lR,  //inputs
			F1R, F2R, F3R, F4R, G1R, G2R, G3R, G4R);


		fL1 = F1L*nx + G1L*ny;
		fL2 = F2L*nx + G2L*ny;
		fL3 = F3L*nx + G3L*ny;
		fL4 = F4L*nx + G4L*ny;

		fR1 = F1R*nx + G1R*ny;
		fR2 = F2R*nx + G2R*ny;
		fR3 = F3R*nx + G3R*ny;
		fR4 = F4R*nx + G4R*ny;

		num_flux_1 = 0.5 * (fL1 + fR1);
		num_flux_2 = 0.5 * (fL2 + fR2);
		num_flux_3 = 0.5 * (fL3 + fR3);
		num_flux_4 = 0.5 * (fL4 + fR4);

	}



}
