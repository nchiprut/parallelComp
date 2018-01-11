#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MAT_SIZE 5

int main(int argc, char **argv) {
	int rank, size;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// make sure the number of processes has root
	float s = sqrt(size);
	if ((ceil(s) -s) != 0) {
		if (rank == 0) {
			printf("%d has no root\n", size);
		}
		MPI_Finalize();
		return 0;
	}
 
	// generate matrix
	srand(rank+1);
	double mat[MAT_SIZE * MAT_SIZE];
	for (int i = 0; i < MAT_SIZE * MAT_SIZE; ++i) {
		mat[i] = (rand() % 100) -50;
	}

	// calc own cordinates
	int x= rank / s;
	int y = ((int)rank % (int)s);

	// open input file and save the matrix
	char szFileName[255] = {0};
	sprintf(szFileName, "input_%d_%d.txt", x,y);
	FILE *f = fopen(szFileName, "w");
	for(int i = 1; i<= MAT_SIZE * MAT_SIZE; ++i){
		fprintf(f, "%lf ", mat[i-1]);
		if(i%(int)MAT_SIZE==0){
			fprintf(f, "\n");
		}
	}
	
	// open output file and save the matrix
	sprintf(szFileName, "output_%d_%d.txt", y,x);
	f = fopen(szFileName, "w");
	for(int i = 1; i<= MAT_SIZE; ++i){
		for(int j = 1; j<= MAT_SIZE; ++j){
			fprintf(f, "%lf ", mat[(j-1)*MAT_SIZE + i-1]);
			if(j%(int)MAT_SIZE==0){
				fprintf(f, "\n");
			}
		}
	}
	MPI_Finalize();
	return 0;
}
	

