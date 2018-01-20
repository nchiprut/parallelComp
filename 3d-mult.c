#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
#include <sys/stat.h>

#define MAT_SIZE 5
#define EPS 0.00000001
#define DIRERCTORY "mult-3d"

void saveMat(char* dir, char* type, int row, int col, double* out) { 

	// open input file and save the matrix
	char szFileName[255] = {0};
	sprintf(szFileName, "%s/mult-3d_%s_%d_%d.txt", dir, type, row, col);
	FILE *f = fopen(szFileName, "w");
	for (int i = 0; i<MAT_SIZE; ++i) {
		for (int j = 0; j<MAT_SIZE; ++j) {
			fprintf(f, "%f ", out[i*MAT_SIZE+j]);
		}
		fprintf(f, "\n ");
	}
	fclose(f);
	
}

void mult(double *a, double *b, double *c) {

	for (int i = 0; i < MAT_SIZE; ++i) {
		for (int j = 0; j < MAT_SIZE; ++j) {
			c[MAT_SIZE*i + j] = 0;

			for (int k = 0; k < MAT_SIZE; ++k) {
				c[MAT_SIZE*i + j] += a[MAT_SIZE*i + k] * b[MAT_SIZE*k + j];
			}
		}
	}
}

void printMat(double *matA, double *matB, double *matC){

	//if ((coords3D[0] == 2)&&(coords3D[1] == 0)&&(coords3D[2] == 0)){
	for (int i = 0; i<MAT_SIZE; ++i) {
		for (int j = 0; j<MAT_SIZE; ++j) {
			printf("%f ", matA[i*MAT_SIZE+j]);
		}
		printf("\n ");
	}
	printf("\n ");
	for (int i = 0; i<MAT_SIZE; ++i) {
		for (int j = 0; j<MAT_SIZE; ++j) {
			printf("%f ", matB[i*MAT_SIZE+j]);
		}
		printf("\n ");
	}
	printf("\n ");
	for (int i = 0; i<MAT_SIZE; ++i) {
		for (int j = 0; j<MAT_SIZE; ++j) {
			printf("%f ", matC[i*MAT_SIZE+j]);
		}
		printf("\n ");
	}
}

void printMatA(double *matA){

	//if ((coords3D[0] == 2)&&(coords3D[1] == 0)&&(coords3D[2] == 0)){
	for (int i = 0; i<MAT_SIZE; ++i) {
		for (int j = 0; j<MAT_SIZE; ++j) {
			printf("%f ", matA[i*MAT_SIZE+j]);
		}
		printf("\n ");
	}
	printf("\n ");
}

int main(int argc, char **argv) {
	int rank, size;
	int ndim = 3;
	int period[3] = {0};
	int dims[3];
	double matC[MAT_SIZE * MAT_SIZE] = {0};
	double matA[MAT_SIZE * MAT_SIZE] = {0};
	double matB[MAT_SIZE * MAT_SIZE] = {0};
	double output[MAT_SIZE * MAT_SIZE] = {0};
	int root = 0;

	
	MPI_Comm comm3D, commRow, commCol, commHeight;
	
	int coords3D[3];
	int id3D;
	int belongs[3];


	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// make sure the number of processes has root
	double s = pow((double)size, (double)1/3);
	if ((ceil(s) -s) >= EPS) {
		if (rank == 0) {
			printf("%d has no root\n", size);
		}
		MPI_Finalize();
		return 0;
	}
	dims[0] = round(s);
	dims[1] = round(s);
	dims[2] = round(s);
	mkdir(DIRERCTORY, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 


	/* Create 3D Cartesian topology for processes */
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, period, 0, &comm3D);
	MPI_Comm_rank(comm3D, &id3D);
	MPI_Cart_coords(comm3D, id3D, ndim, coords3D);


	/* Create 1D row subgrids */
	belongs[0] = 1; belongs[1] = 0; belongs[2] = 0;
	MPI_Cart_sub(comm3D, belongs, &commRow);
	
	/* Create 1D column subgrids */
	belongs[0] = 0; belongs[1] = 1; belongs[2] = 0;
	MPI_Cart_sub(comm3D, belongs, &commCol);

	/* Create 1D height subgrids */
	belongs[0] = 0; belongs[1] = 0; belongs[2] = 1;
	MPI_Cart_sub(comm3D, belongs, &commHeight);


	MPI_Request requests[2];
	MPI_Status status[2];

	int originRank;

	srand(rank+2);
	if (coords3D[1] == 0) {
		// generate matrix
		for (int i = 0; i < MAT_SIZE * MAT_SIZE; ++i) {
			matB[i] = (rand() % 100) -50;
		}
		saveMat(DIRERCTORY, "input_B" , coords3D[2], coords3D[0], matB);

		MPI_Cart_rank(commCol, coords3D + 1, &originRank);
		MPI_Ibcast(matB, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commCol, &requests[1]);
		
	} else {
		MPI_Cart_rank(commCol, &root, &originRank);
		MPI_Ibcast(matB, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commCol, &requests[1]);

	}
	if (coords3D[0] == 0) {
		// generate matrix
		for (int i = 0; i < MAT_SIZE * MAT_SIZE; ++i) {
			matA[i] = (rand() % 100) -50;
		}
		saveMat(DIRERCTORY, "input_A" , coords3D[1], coords3D[2], matA);

		MPI_Cart_rank(commRow, coords3D, &originRank);
		MPI_Ibcast(matA, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commRow, &requests[0]);

	} else {
		MPI_Cart_rank(commRow, &root, &originRank);
		MPI_Ibcast(matA, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commRow, &requests[0]);

	}

	MPI_Wait(&requests[1], &status[1]);
	MPI_Wait(&requests[0], &status[0]);

	mult(matA, matB, matC);

	/*
	if ((coords3D[0] == 1) && (coords3D[1] == 0) && (coords3D[2] == 0)) {
		printMat(matA, matB, matC);
	}
	*/

	/********************* REDUCE ********************/


	// the one with 0 on height axis will reduce all
	MPI_Cart_rank(commHeight, &root, &originRank);
	MPI_Reduce(matC, output, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, MPI_SUM, originRank, commHeight);


	/*
	if ((coords3D[0] == 1) && (coords3D[1] == 1) && (coords3D[2] == 1)) {
		printMatA(matC);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if ((coords3D[0] == 1) && (coords3D[1] == 1) && (coords3D[2] == 0)) {
		printMatA(matC);
		printMatA(output);
	}
	*/


	mkdir(DIRERCTORY, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (coords3D[2] == 0) {
		saveMat(DIRERCTORY, "output" , coords3D[0], coords3D[1], output);
	}

	MPI_Finalize();
	return 0;
}
	
