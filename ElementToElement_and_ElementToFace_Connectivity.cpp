#include "Define_All_Includes_Header.h"
#include "Define_Arrays_Header.h"

void ElementToElement_and_ElementToFace_Connectivity(int const K, int const Nv, int const Nfaces, int ** EToV, int **EToE, int **EToF)
{

	// Purpose: triangle face connect algorithm (original by Toby Isaac). 
	// This function builds the EToE (Element to Element) and EToF (Element
	// to Face) connection information.  In other words, EToE[e][f] returns
	// the element number connected to e on face f, and EToF[e][f] returns the
	// face number for that neighboring element. 

	int **EToV2, **EToE2, **EToF2;

	EToV2 = i2array(Nfaces, K); EToE2 = i2array(Nfaces, K); EToF2 = i2array(Nfaces, K);

	//Mahya: modified part
	for (int i = 0; i < Nfaces; i++)
	{
		parallel_for (0, K , [&](int j)
		{
			EToV2[i][j] = EToV[j + 1][i + 1];
		});
	}
	//

	int **matched, *id;
	id = new int[(Nfaces)* (K)];
	matched = new int *[Nfaces];
	for (int i = 0; i<Nfaces; i++){
		matched[i] = new int[K];
		parallel_for(0, K, [&](int j)
		{
			EToE2[i][j] = j;
			EToF2[i][j] = i;
			matched[i][j] = 0;
		});
	}

	parallel_for(0, K, [&](int e)
	{
		double node0, node1;
		if (EToV2[0][e] < EToV2[1][e]){
			node0 = EToV2[0][e];
			node1 = EToV2[1][e];
		}
		else {
			node0 = EToV2[1][e];
			node1 = EToV2[0][e];
		}
		id[e] = node0 * Nv + node1 + 1;
		if (EToV2[1][e] < EToV2[2][e]){
			node0 = EToV2[1][e];
			node1 = EToV2[2][e];
		}
		else {
			node0 = EToV2[2][e];
			node1 = EToV2[1][e];
		}
		id[e + K] = node0 * Nv + node1 + 1;
		if (EToV2[2][e] < EToV2[0][e]){
			node0 = EToV2[2][e];
			node1 = EToV2[0][e];
		}
		else {
			node0 = EToV2[0][e];
			node1 = EToV2[2][e];
		}
		id[e + 2 * K] = node0 * Nv + node1 + 1;
	});
	// Scan id list for matching edge 
	parallel_for(0, K, [&](int e)
	{
		int myid;
		for (int f = 0; f<Nfaces; f++){
			myid = id[f*K + e];
			if (matched[f][e] == 0){
				for (int e2 = 0; e2<K; e2++){
					for (int f2 = 0; f2<Nfaces; f2++){
						if ((myid == id[f2*K + e2]) && (e != e2)){
							EToE2[f][e] = e2;
							EToF2[f][e] = f2;
							EToE2[f2][e2] = e;
							EToF2[f2][e2] = f;
							matched[f2][e2] = 1;
							f2 = Nfaces;
							e2 = K;
						}
					}
				}
			}
		}
	});

	//Mahya: modified part
	for (int f = 0; f < Nfaces; f++)
	{
		parallel_for(0, K, [&](int e)
		{
			EToE[e + 1][f + 1] = EToE2[f][e] + 1;
			EToF[e + 1][f + 1] = EToF2[f][e] + 1;
		});
	}
	//

	for (int i = 0; i <= Nfaces - 1; i++)
	{
		delete[] matched[i];
		//Mahya: modified part
		delete[] EToV2[i];
		delete[] EToE2[i];
		delete[] EToF2[i];
		//
	}
	delete[] matched;
	//Mahya: modified part
	delete[] EToV2;
	delete[] EToE2;
	delete[] EToF2;
	//
	delete[] id;
}
