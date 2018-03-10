#include "Define_All_includes_Header.h"

void x_and_y(int const Np, int const K, double const *r, double const *s, int **EToV, double const *Vx, double const * Vy, double ** x, double ** y)
{
	parallel_for(1, K + 1, [&](int j)
	{
		for (int i = 1; i <= Np; i++)
		{
			x[j][i] = 0.5*(-(r[i] + s[i])*Vx[EToV[j][1]] + (1 + r[i])*Vx[EToV[j][2]] + (1 + s[i])*Vx[EToV[j][3]]);
			y[j][i] = 0.5*(-(r[i] + s[i])*Vy[EToV[j][1]] + (1 + r[i])*Vy[EToV[j][2]] + (1 + s[i])*Vy[EToV[j][3]]);
			//cout << i << "   " << j << "   " << x[i][j] << "   " << y[i][j] << endl;
		}
	});
	//system("pause");
}