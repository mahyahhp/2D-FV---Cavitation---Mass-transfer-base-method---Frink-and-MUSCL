

void FluxVectors(double const p_temp, double const u_temp, double const v_temp, double const rho_temp, double const alpha_l_temp,
	double &F1, double &F2, double &F3, double &F4, double &G1, double &G2, double &G3, double &G4)
{
	F1 = u_temp;
	F2 = rho_temp * u_temp * u_temp + p_temp;
	F3 = rho_temp * u_temp * v_temp;
	F4 = alpha_l_temp*u_temp;

	G1 = v_temp;
	G2 = rho_temp * u_temp * v_temp;
	G3 = rho_temp * v_temp * v_temp + p_temp;
	G4 = alpha_l_temp*v_temp;
}