#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define MAT_SIZE 5
#define EPS 0.00000001

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
}

int main(int argc, char **argv) {
	int rank, size;
	int ndim = 3;
	int period[3] = {0};
	int dims[3];
	double matC[MAT_SIZE * MAT_SIZE] = {0};
	double matA[MAT_SIZE * MAT_SIZE] = {0};
	double matB[MAT_SIZE * MAT_SIZE] = {0};

	
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
 


	/* Create 2D Cartesian topology for processes */
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

	int origin[3];
	int originRank;
	origin[0] = coords3D[0];
	origin[1] = coords3D[1];
	origin[2] = coords3D[2];

	srand(rank+2);
	if (coords3D[1] == 0) {
		// generate matrix
		for (int i = 0; i < MAT_SIZE * MAT_SIZE; ++i) {
			matB[i] = (rand() % 100) -50;
		}

		MPI_Cart_rank(commCol, coords3D + 1, &originRank);
		MPI_Ibcast(matB, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commCol, &requests[1]);
		
	} else {
		origin[1] = 0;
		MPI_Cart_rank(commCol, origin + 1, &originRank);
		MPI_Ibcast(matB, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commCol, &requests[1]);
		origin[1] = coords3D[1];

	}
	if (coords3D[0] == 0) {
		// generate matrix
		for (int i = 0; i < MAT_SIZE * MAT_SIZE; ++i) {
			matA[i] = (rand() % 100) -50;
		}
		MPI_Cart_rank(commRow, coords3D, &originRank);
		MPI_Ibcast(matA, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commRow, &requests[0]);

	} else {
		origin[0] = 0;
		MPI_Cart_rank(commRow, origin, &originRank);
		MPI_Ibcast(matA, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, originRank, commRow, &requests[0]);
		origin[0] = coords3D[0];

	}

	MPI_Wait(&requests[1], &status[1]);
	MPI_Wait(&requests[0], &status[0]);

	mult(matA, matB, matC);

	/*
	int origin[3] = {0};
	int originRank;
	int rslt;
	MPI_Cart_rank(comm3D, origin, &originRank);

	if (originRank == id3D) {
		//printf("my rank: %d,%d,%d", coords3D[0], coords3D[1], coords3D[2]);
		MPI_Bcast(mat, MAT_SIZE*MAT_SIZE, MPI_DOUBLE, id3D, commRow);

	} else {
		rslt = MPI_Bcast(buf, MAT_SIZE * MAT_SIZE, MPI_DOUBLE, originRank, commRow);
		//printf("my rank: %d,%d,%d, success: %d\n", coords3D[0], coords3D[1], coords3D[2], rslt);


	}
	printf("my rank: %d,%d,%d, success: %d\n", coords3D[0], coords3D[1], coords3D[2], rslt);
	printf("%d", commRow);
	for (int i = 0; i<MAT_SIZE; ++i) {
		for (int j = 0; j<MAT_SIZE; ++j) {
			printf("%f ", buf[i*MAT_SIZE+j]);
		}
		printf("\n ");
	}
	*/

	if ((coords3D[0] == 1) && (coords3D[1] == 0) && (coords3D[2] == 0)) {
		printMat(matA, matB, matC);
	}


	MPI_Finalize();
	return 0;
}
	
