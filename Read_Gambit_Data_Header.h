
#ifndef _Read_Gambit_
#define _Read_Gambit_
extern void read_data_control_info(int & Nnod, int & Nelem, int & NBCs);
extern void read_main_data(double *x, double *y, int ** Element, double * Area, double **Length_Face, 
	int *BC_NumElem, int **BC_Type_Elem, int **BC_Type_Edge, int *No_CodeCell, int **NodeCell, int **NodeNode,
	int *BC_Elem_ID, int *BC_Edge_Counter, int **BC_Edge, int **BC_Edge_type);
#endif