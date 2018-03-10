#include "Define_All_Includes_Header.h"

//double *d1array(int a){
//	double *result = new double[a];
//	return result;
//}
//double **d2array(int a, int b){
//	double **result;
//	result = new double*[a];
//	for (int i = 0; i <= a - 1; i++)
//		result[i] = new double[b];
//	return result;
//}
//double ***d3array(int a, int b, int c){
//	double ***result;
//	result = new double**[a];
//	for (int i = 0; i <= a - 1; i++){
//		result[i] = new double*[b];
//		for (int j = 0; j <= b - 1; j++){
//			result[i][j] = new double[c];
//		}
//	}
//	return result;
//}
//int *i1array(int a){
//	int *result = new int[a];
//	return result;
//}
//int **i2array(int a, int b){
//	int **result;
//	result = new int*[a];
//	for (int i = 0; i <= a - 1; i++)
//		result[i] = new int[b];
//	return result;
//}
//int ***i3array(int a, int b, int c){
//	int ***result;
//	result = new int**[a];
//	for (int i = 0; i <= a - 1; i++){
//		result[i] = new int*[b];
//		for (int j = 0; j <= b - 1; j++){
//			result[i][j] = new int[c];
//		}
//	}
//	return result;
//}






double *d1array(int a){
	double *result = new double[a];
	return result;
}
double **d2array(int a, int b){
	double **result;
	result = new double*[a];
	parallel_for(0, a, [&](int i)
	{
		result[i] = new double[b];
	});
	return result;
}
double ***d3array(int a, int b, int c){
	double ***result;
	result = new double**[a];
	parallel_for(0, a, [&](int i)
	{
		result[i] = new double*[b];
		for (int j = 0; j <= b - 1; j++){
			result[i][j] = new double[c];
		}
	});
	return result;
}
int *i1array(int a){
	int *result = new int[a];
	return result;
}
int **i2array(int a, int b){
	int **result;
	result = new int*[a];
	parallel_for(0, a, [&](int i)
	{
		result[i] = new int[b];
	});
	return result;
}
int ***i3array(int a, int b, int c){
	int ***result;
	result = new int**[a];
	parallel_for(0, a, [&](int i)
	{
		result[i] = new int*[b];
		for (int j = 0; j <= b - 1; j++){
			result[i][j] = new int[c];
		}
	});
	return result;
}