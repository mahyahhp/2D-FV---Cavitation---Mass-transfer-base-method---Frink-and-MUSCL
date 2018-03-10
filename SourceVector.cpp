void SourceVector(double const m_minus, double const m_plus, double const rho_v, //inputs
	double &S1, double &S2, double &S3, double &S4)
{
	S1 = 0.999*S1 + 0.001*(m_plus + m_minus)*(1 - 1 / rho_v);
	S2 = 0.999*S2 + 0.001*(0);
	S3 = 0.999*S3 + 0.001 * (0);
	S4 = 0.999*S4 + 0.001*(m_plus + m_minus)*(1);

	//S1 = (m_plus + m_minus)*(1 - 1 / rho_v);
	//S2 = (0);
	//S3 = (0);
	//S4 = (m_plus + m_minus)*(1);
}