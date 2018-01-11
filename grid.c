#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

#define VEC_SIZE 5

int getRank(int x, int y, int root) {
	x += root;
	y += root;
	return (x % root)*root + (y%root);
}


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
 
	// generate vector
	srand(rank+1);
	double vec[VEC_SIZE];
	for (int i = 0; i < VEC_SIZE; ++i) {
		vec[i] = (rand() % 100) -50;
	}

	// calc own cordinates
	int x= rank / s;
	int y = ((int)rank % (int)s);

	// open input file and save the matrix
	char szFileName[255] = {0};
	sprintf(szFileName, "grid_input_%d_%d.txt", x,y);
	FILE *f = fopen(szFileName, "w");
	for(int i = 1; i<= VEC_SIZE; ++i){
		fprintf(f, "%lf ", vec[i-1]);
	}
	
	double buf[VEC_SIZE*4];

	int neighb[4] = {getRank(x,y+1,s),getRank(x+1,y,s),getRank(x,y-1,s),getRank(x-1,y,s)};
	MPI_Request requests[4];
	MPI_Request send_requests[4];
	MPI_Status status[4];

	for (int i = 0; i < 4; ++i) {

		printf("my rank: %d, sent to: %d\n", rank, neighb[i]);
		fflush(stdout);
		MPI_Isend(vec, VEC_SIZE, MPI_DOUBLE, neighb[i], 0, MPI_COMM_WORLD, &send_requests[i]);
		MPI_Irecv(buf+(i*VEC_SIZE), VEC_SIZE, MPI_DOUBLE, neighb[i], 0, MPI_COMM_WORLD, &requests[i]);
	}
	for (int i = 0; i < 4; ++i) {

		MPI_Wait(&requests[i], &status[i]);
	}
	if (rank == 0) {
		for(int i = 1; i<= VEC_SIZE*4; ++i){
			printf("%lf ", buf[i-1]);
			if ((i%5) ==0){
				printf("\n");
			}

		}
		
	}
	
	/*
	// open output file and save the matrix
	sprintf(szFileName, "grid_output_%d_%d.txt", y,x);
	f = fopen(szFileName, "w");
	for(int i = 1; i<= VEC_SIZE; ++i){
		fprintf(f, "%lf ", vec[i-1]);
	}
	*/

	MPI_Finalize();
	return 0;
}
	

