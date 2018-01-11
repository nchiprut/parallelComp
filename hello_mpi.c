#include <stdio.h>
#include "mpi.h"

int main(int argc, char **argv) {
	int rank, size;
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	printf("my rank: %d, total: %d\n", rank, size);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		printf("Hello World!\n");

	}

	MPI_Finalize();
	return 0;
}
	

