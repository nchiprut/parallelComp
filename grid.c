#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"
#include <sys/stat.h>

#define VEC_SIZE 5
#define DIRERCTORY "grid"
#define EPS 0.00000001

void saveVec(char* dir, char* type, int row, int col, double* out) { 

	// open input file and save the matrix
	char szFileName[255] = {0};
	sprintf(szFileName, "%s/mult-3d_%s_%d_%d.txt", dir, type, row, col);
	FILE *f = fopen(szFileName, "w");
	for(int i = 0; i< VEC_SIZE; ++i){
		fprintf(f, "%lf ", out[i]);
	}
	fclose(f);
	
}

void conv(double* input, double* output) {
	for (int i=0; i<VEC_SIZE; ++i) {
		for (int j=1; j<=4; ++j) {
			output[i] += j*input[i+(j-1)*VEC_SIZE];

		}
	}
}

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
	if ((ceil(s) -s) >= EPS) {
		if (rank == 0) {
			printf("%d has no root\n", size);
		}
		MPI_Finalize();
		return 0;
	}
	mkdir(DIRERCTORY, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
 
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
	saveVec(DIRERCTORY, "input", x, y, vec);
	
	double buf[VEC_SIZE*4];

	int neighb[4] = {getRank(x-1,y,s),getRank(x+1,y,s),getRank(x,y-1,s),getRank(x,y+1,s)};
	MPI_Request requests[4];
	MPI_Request send_requests[4];
	MPI_Status status[4];

	for (int i = 0; i < 4; ++i) {

		MPI_Isend(vec, VEC_SIZE, MPI_DOUBLE, neighb[i], 0, MPI_COMM_WORLD, &send_requests[i]);
		MPI_Irecv(buf+(i*VEC_SIZE), VEC_SIZE, MPI_DOUBLE, neighb[i], 0, MPI_COMM_WORLD, &requests[i]);
	}
	for (int i = 0; i < 4; ++i) {

		MPI_Wait(&requests[i], &status[i]);
	}
	conv(buf, vec);
	saveVec(DIRERCTORY, "output", x, y, vec);
	

	MPI_Finalize();
	return 0;
}
	

