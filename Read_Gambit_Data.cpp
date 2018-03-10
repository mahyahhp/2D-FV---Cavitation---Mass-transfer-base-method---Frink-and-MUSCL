#include "Define_All_Includes_Header.h"
#include "Define_Arrays_Header.h"

ifstream file;

void read_data_control_info(int & Nnod, int & Nelem, int &NBCs)
//(int NUMNP,int NELEM)
{
	file.open("input.neu");
	//******* CONTROL INFO *******//
	char A[10][80];

	file >> A[1] >> A[2] >> A[3];
	file >> A[1] >> A[2] >> A[3] >> A[4];
	file >> A[1];
	file >> A[1] >> A[2] >> A[3] >> A[4];
	file >> A[1] >> A[2] >> A[3] >> A[4] >> A[5] >> A[6];

	int NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL;
	file >> NUMNP >> NELEM >> NGRPS >> NBSETS >> NDFCD >> NDFVL;
	file >> A[1];

	Nnod = NUMNP;
	Nelem = NELEM;
	NBCs = NBSETS;
	//****************************//
	file.close();
}
void read_main_data(double *x, double *y, int ** Element, double * Area, double **Length_Face, 
	int *BC_NumElem, int **BC_Type_Elem, int **BC_Type_Edge, int *No_NodeCell, int **NodeCell, int **NodeNode,
	int *BC_Elem_ID, int *BC_Edge_Counter, int **BC_Edge, int **BC_Edge_type)
{
	file.open("input.neu");
	//******* CONTROL INFO *******//
	char A[10][80];

	file >> A[1] >> A[2] >> A[3];
	file >> A[1] >> A[2] >> A[3] >> A[4];
	file >> A[1];
	file >> A[1] >> A[2] >> A[3] >> A[4];
	file >> A[1] >> A[2] >> A[3] >> A[4] >> A[5] >> A[6];

	int NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL;
	file >> NUMNP >> NELEM >> NGRPS >> NBSETS >> NDFCD >> NDFVL;
	file >> A[1];
	//****************************//

	//******* NODAL COORDINATE *******//
	file >> A[1] >> A[2] >> A[3];


	for (int j = 1; j <= NUMNP; j++){
		int i;
		file >> i; //shomareye noghte
		file >> x[j]; //mokhtasate x noghteh
		file >> y[j];
	}
	file >> A[1];
	//********************************//

	//******* ELEMENTS/CELLS *******//
	file >> A[1] >> A[2];

	int *element_type_code, *NDF, **IE;
	double *Area_element;
	element_type_code = i1array(NELEM + 1); NDF = i1array(NELEM + 1); IE = i2array(NELEM + 1, 20 + 1); Area_element = d1array(NELEM + 1);
	int k;

	for (int j = 1; j <= NELEM; j++)
	{
		int i;
		file >> i; //shomareye eleman
		file >> element_type_code[i]; //codi ke noe eleman ra moshakhas mikonad. masalan 3 yani triangle
		file >> NDF[i]; //number of degree of freedom.
		for (k = 1; k <= NDF[i]; k++)
		{
			file >> IE[i][k];
			Element[i][k] = IE[i][k]; //@@@@@@@@@@@@@@@@
			//cout << i << "   " << k << "   " << Element[i][k] << endl;
		}
	}

	parallel_for(1, NELEM + 1, [&](int i)
	{
		for (int k = 1; k <= NDF[i]; k++)
		{
			No_NodeCell[Element[i][k]] = 0;
			for (int j = 1; j <= NELEM; j++)
			{
				for (int m = 1; m <= NDF[i]; m++)
				{
					//if (x[Element[i][k]] == x[Element[j][m]] & y[Element[i][k]] == y[Element[j][m]] & j!=i)
					if (x[Element[i][k]] == x[Element[j][m]] & y[Element[i][k]] == y[Element[j][m]])
					{
						No_NodeCell[Element[i][k]] = No_NodeCell[Element[i][k]] + 1; //Noumber of Connected Elements to a Node
						NodeCell[Element[i][k]][No_NodeCell[Element[i][k]]] = j; //شماره الملنی که به گره متصل است elemeny number that is connect through a node
						NodeNode[Element[i][k]][No_NodeCell[Element[i][k]]] = m;
						//cout << "Element = " << i << " , Node = " << k << " , Node number = " << Element[i][k] <<
						//	" , NodeCell = " << NodeCell[Element[i][k]][No_NodeCell[k][i]] << " , NodeNode = " << NodeNode[Element[i][k]][No_NodeCell[k][i]] << endl;
					}
				}
			}
		}
	});

	for (int i = 1; i <= NELEM; i++)
	{
		Area_element[i] = (x[Element[i][2]] * y[Element[i][3]] + x[Element[i][1]] * y[Element[i][2]] +
			x[Element[i][3]] * y[Element[i][1]] - x[Element[i][2]] * y[Element[i][1]] - //@@@@@@@@@@@@@@@@@@@@@@@
			x[Element[i][3]] * y[Element[i][2]] - x[Element[i][1]] * y[Element[i][3]]) / 2;
		Area[i] = Area_element[i];

		Length_Face[i][1]=sqrt(pow(x[Element[i][2]]-x[Element[i][1]],2)+pow(y[Element[i][2]]-y[Element[i][1]],2));
		Length_Face[i][2]=sqrt(pow(x[Element[i][2]]-x[Element[i][3]],2)+pow(y[Element[i][2]]-y[Element[i][3]],2));
		Length_Face[i][3]=sqrt(pow(x[Element[i][1]]-x[Element[i][3]],2)+pow(y[Element[i][1]]-y[Element[i][3]],2));
		//cout<<Area_element[i]<<endl;
		//cout << i << "   " << x[Element[i][2]] << "   " << x[Element[i][1]] << "   " << y[Element[i][2]] << "   " << y[Element[i][1]] << "   " <<Length_Face[i][1]<< endl;
	}
		//system("pause");
	
	file >> A[1];
	//********************************//

	//******* ELEMENTS GROUP *******//
	for (int i = 1; i <= 13 + NELEM + 1; i++){
		file >> A[1]; //shomareye eleman
		if (A[1] == "ENDOFSECTION")
		{
			break;
		}
	}
	//********************************//



	//******* Boundary Conditions *******//
	char name [80][80];
	int *BC_data_type,**BC_element,**BC_edge;
	double **L_edge;
	//int *ND;
	//ND=i1array(NBSETS+1);
	BC_data_type=i1array(NBSETS+1);BC_element=i2array(NBSETS+1,NUMNP+1);
	BC_edge=i2array(NBSETS+1,NUMNP+1);L_edge=d2array(NBSETS+1,NUMNP+1);


	for (int j = 1; j <= NELEM; j++)
	{
		BC_Elem_ID[j] = 0; //moshakhas mikonad ke aya eleman rooye marz gharar darad ya kheyr
		BC_Edge_Counter[j] = 0; //tedade edge haii az eleman ke rooye marz gharar darand ra mishomarad
	}


	for (int i = 1; i <= NBSETS; i++){ // boundary condition set number
		file>>A[1]>>A[2]>>A[3];
		file>>name[i];
		//cout << name[i] << endl;
		file>>BC_data_type[i]; //(0=node , 1=element/cell)
		//file>>ND[i];//number of data ya tedade kole dadehaii ke dar in majmoee az B.C. ha vojoud darad
		file>>BC_NumElem[i];//number of data ya tedade kole dadehaii ke dar in majmoee az B.C. ha vojoud darad

		//BC_Npoint[i]=ND[i]; //@@@@@@@@@@@@@

		file>>A[1];  //"kole tedade node ya element haii ke dar in majmoe eraee shode ast "   "tedade maghadir baraye har dadeh zakhire shode (0 yani hich meghdari be marz ha nesbat dadeh nashodeh ast)"
		file>>A[1];  //"code boundary condition (6 yani element_side)"
		//for (j=1;j<=ND[i];j++){ //tedade kole edge haye marzi az dasteye i om
		for (int j = 1; j <= BC_NumElem[i]; j++){ //tedade kole edge haye marzi az dasteye i om
			file>>BC_element[i][j]; //BC element number
			file>>A[1]; //in sotoon ezafe ast chon dae ELEMENTS/CELLS noe eleman ra fahmidim
			file>>BC_edge[i][j]; //nth edge of boundary element
			BC_Type_Edge[j][i]=BC_edge[i][j]; //moshakhas mikonad ke face/edge shomare chand az "j"-omin eleman az "i"-omin dasteye sharte marzi rooye marz gharar darad

			BC_Type_Elem[j][i]=BC_element[i][j]; //moshakhas mikonad ke elemane shomare chand az "j"-omin eleman az "i"-omin dasteye sharte marzi rooye marz gharar darad

			BC_Elem_ID[BC_Type_Elem[j][i]] = 1;
			BC_Edge_Counter[BC_Type_Elem[j][i]] = BC_Edge_Counter[BC_Type_Elem[j][i]] + 1;
			BC_Edge[BC_Type_Elem[j][i]][BC_Edge_Counter[BC_Type_Elem[j][i]]] = BC_edge[i][j]; //BC_Edge[j][i]: elemane "j"-om, "i"-omin edge maezii an edge chandom ast
			BC_Edge_type[BC_Type_Elem[j][i]][BC_Edge_Counter[BC_Type_Elem[j][i]]] = i; //BC_Edge[j][i]: elemane "j"-om, "i"-omin edge maezii an sharte marzi chandom ast

			//cout<<" i = "<<i<< " j = "<<j<<"  BC_elem = "<<BC_Type_Elem[i][j]<<"  BC_edge = "<<BC_edge[i][j]<<endl;
		}
		//system("pause");
		file>>A[1];
	}
	//****************************//
	//system("pause");

	file.close();
}